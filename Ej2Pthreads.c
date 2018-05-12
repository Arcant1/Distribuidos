#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

#define COMPARAR_SECUENCIAL

// Tipo de los elementos en los vectores
// Compilar con -D_INT_ para vectores de tipo entero
// Compilar con -D_DOUBLE_ para vectores de tipo double
// Predeterminado vectores de tipo float

#ifdef _INT_
typedef int basetype;     // Tipo para elementos: int
#define tipo_dato   "ints"
#elif _DOUBLE_
typedef double basetype;  // Tipo para elementos: double
#define tipo_dato    "doubles"
#else
typedef float basetype;   // Tipo para elementos: float     PREDETERMINADO
#define tipo_dato    "floats"
#endif

// Tipo de dato que se pasa de parámetro a los threads
typedef struct param
{
	int id;
} param;


// -- Variables globales --

// Matrices
basetype *A;
basetype *B;
basetype *C;	// Matriz resultado

int cant_filas_restantes;	// Cantidad de filas que restan ordenar
int cant_filas_x_thread;	// Cantidad de filas que procesa cada thread en cada iteración (va disminuyendo)
int pos_max_global;

int N;
int CANT_THREADS;

double speedup;
double eficiencia;
double tiempo_paral;
double tiempo_sec;

basetype factor_global;

pthread_barrier_t barrera; 		// Barrerra

// -- Definición de funciones --

void *funcion_threads(void *arg);
void multiplicacion(param *parametro);
void promedioB(param* parametro);
void prodPromLU(param* parametro);
void imprimir_matriz (basetype * matriz,int N);
double dwalltime();
double tiempo_copia_total=0;

#ifdef COMPARAR_SECUENCIAL
void multiplicacion_secuencial(basetype *A,basetype *B,basetype *C_secuencial,int N);

void verificar_resultado(basetype *C,basetype *C_secuencial,int N);

basetype *C_secuencial;	// Matriz resultado
#endif


int main(int argc,char *argv[])
{
	int i;
	int j;

	if (argc != 3) {
		printf("Error en la cantidad de parámetros. Parametro 1 -> long matriz Parametro 2 -> n° threads");
		printf("Ingresar la dimension de la matriz!\n");
		return 0;
	}


	N = atoi(argv[1]);	// Dimensión de la matriz: N*N
	CANT_THREADS = atoi(argv[2]);


	printf("Dimensión de la matriz: %d*%d \n",N,N);

	//Vector usado para guardar las sumas parciales de la matriz B para luego hacer el promedio
	double * sumaParcialB = (basetype*)malloc(sizeof(basetype)*CANT_THREADS);
	double promB;

	// Reserva de memoria para las matrices
	A=(basetype*)malloc(sizeof(basetype)*N*N);			// Reserva memoria para A
	B=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para B
	C=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para C
	D=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para D
	E=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para E	
	F=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para F
	L=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para L
	U=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para U	

	#ifdef COMPARAR_SECUENCIAL
	C_secuencial=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para C_SECUENCIAL
	#endif


	// -- Inicialización ALEATORIA de las matrices --

	// Inicializa matrices A y B
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j]=rand()%5;  	// Inicializa matriz A con random
			B[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			C[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			D[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			E[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			F[i*N+j]=rand()%5; 	// Inicializa matriz B con random
		}
	}

	//inicializo las matrices triangulares L y U
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			if(i>j)
			{
				U[i*N+j]=rand%5;
				L[i*N+j]=0;	
			}
			else
			{
				L[i*N+j]=rand%5;
				U[i*N+j]=0;	
			}
		}
	}


	param parametros[CANT_THREADS];	// Arreglo de param (struct que contiene los datos para pasar a los threads)

	// -- Inicialización de threads --
	pthread_t threads[CANT_THREADS];	// Arreglo de threads

	if (pthread_barrier_init(&barrera,NULL,CANT_THREADS)!=0) {
		printf("Error creacion de barrera\n");
		return 0;
	}

	// Inicialización de parámetros
	for (i=0;i<CANT_THREADS;i++){
		parametros[i].id=i;
	}
	double tiempo_inicial=dwalltime();
	// Creacion de los threads
	for (i=0;i<CANT_THREADS;i++){
		if (pthread_create( &threads[i],NULL,funcion_threads,(void*)&parametros[i])!=0) {
			printf("Error creacion de thread\n");
			return 0;
		}
	}

	// Join de los threads
	for(i = 0; i < CANT_THREADS; i++) 
	{
		pthread_join(threads[i], NULL);

	}

	tiempo_paral = dwalltime()-tiempo_inicial;
	printf("\nTiempo Total (pthreads) : %f\n\n",dwalltime()-tiempo_inicial);

	#ifdef COMPARAR_SECUENCIAL
		tiempo_inicial=dwalltime();
		multiplicacion_secuencial(A,B,C_secuencial,N);	// C_secuencial = A*B
		tiempo_sec = dwalltime()-tiempo_inicial;
		printf("-- Fin de multiplicacion (secuencial) -->> \t Tiempo: %f \n", tiempo_sec);

		speedup = tiempo_sec / tiempo_paral;
		printf("-- Speedup conseguido: %f \n", speedup);
		eficiencia = speedup / CANT_THREADS;
		printf("-- Eficiencia: %f \n", eficiencia);

	free(C_secuencial);
	#endif


	// Libera memoria
	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	free(F);
	free(L);
	free(U);

	return(0);

}

// ------------------------
// -- FUNCIONES PTHREADS //
// ------------------------


void *funcion_threads(void *arg) {
	param* parametro = (param*)arg;
	double tiempo_inicial2;
	//printf("Mi ID es: %d \n",(*parametro).id);
	if ((*parametro).id==0){
		tiempo_inicial2=dwalltime();
	}
	multiplicacion(parametro);
	pthread_barrier_wait(&barrera); //espero a que todos los hilos finalicen
	if ((*parametro).id==0){
		printf("-- Fin de multiplicacion (pthread) -->> \t Tiempo: %f \n",dwalltime()-tiempo_inicial2);
		//printf("Matriz A:\n");
	}
	pthread_exit(NULL);
}

void multiplicacion(param* parametro)
{
	//printf("Comienzo etapa 1\n");
	int id = (*parametro).id;
	basetype total;
	int i,j,k;


	// Filas que multiplica el thread
	int cant_filas = N/CANT_THREADS;	// Cant de filas que multiplica cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	//printf("ID: %d \t fila_inicial=%d \t fila_final=%d \n",id,fila_inicial,fila_final);

	// Multiplica A*B=C
	for(i=fila_inicial;i<=fila_final;i++)
	{	// Recorre solo algunas filas
		for(j=0;j<N;j++)
		{	// Recorre todas las columnas
			total=0;
			for(k=0;k<N;k++)
			{
				total+=A[i*N+k]*B[k*N+j];	// total=A*B
			}
			C[i*N+j] = total;
		}
	}
}



void prodPromLU(param* parametro)
{
	int id = (*parametro).id;
	basetype total;
	int i,j,k;


	// Filas que suma el thread
	int cant_filas = N/CANT_THREADS;	// Cant de filas que suma cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	total=0;

	for(i=fila_inicial;i<=fila_final;i++)
	{	// Recorre solo algunas filas
		for(j=0;j<N;j++)
		{	// Recorre todas las columnas
			total=0;
			for(k=0;k<N;k++)
			{
				total+=A[i*N+k]*B[k*N+j];	// total=A*B
			}
			C[i*N+j] = total;
		}
	}
}


void promedioB(param* parametro)
{
	int id = (*parametro).id;
	basetype total;
	int i,j,k;


	// Filas que multiplica el thread
	int cant_filas = N/CANT_THREADS;	// Cant de filas que multiplica cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	//printf("ID: %d \t fila_inicial=%d \t fila_final=%d \n",id,fila_inicial,fila_final);
	total=0;
	// Multiplica A*B=C
	for(i=fila_inicial;i<=fila_final;i++)
	{										// Recorre solo algunas filas
		for(j=0;j<N;j++)					// Recorre todas las columnas
		{	
			total+=B[i*N + j];
		}
	}
	sumaParcialB[id]=total;
    pthread_barrier_wait(&barrera); //Espera a que todas terminen
    if(id == 0)						//si es el hilo principal
    {
    	for (int i = 0; i < CANT_THREADS; ++i)
    	{
    		promB+=sumaParcialB[i];
    	}
    	promB /=CANT_THREADS;
    }
}

// --------------------------
// -- FUNCIONES AUXILIARES //
// --------------------------

void imprimir_matriz (basetype * matriz,int N){
	int i;
	int j;
	for (i=0;i<N;i++){
		for (j = 0 ; j < N ; j++){
			printf ("%.1f\t",matriz [ i * N + j ]);
		}
		printf("\n");
	}
}

double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

// ----------------------------
// -- FUNCIONES SECUENCIALES //
// ----------------------------

#ifdef COMPARAR_SECUENCIAL
void multiplicacion_secuencial(basetype *A,basetype *B,basetype *C,int N){
	//printf("Comienzo etapa 1\n");

	basetype total;
	int i,j,k;

	// Multiplica A*B*D=C
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			total=0;
			for(k=0;k<N;k++){
					total+=A[i*N+k]*B[k*N+j];	// total=A*B
				}
				C[i*N+j] = total;		// C=total
			}
		}

	}

	void verificar_resultado(basetype *C,basetype *C_secuencial,int N){
		int i,j;
		int check = 1;
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				check=check&&(C[i*N+j]==C_secuencial[i*N+j]);
			}
		}

		if(check){
			printf("Resultado correcto\n");
		}
		else{
			printf("Resultado erroneo \n");
		}
	}
#endif
