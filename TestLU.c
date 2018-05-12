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

double dwalltime();
void *funcion_threads(void *arg);
void imprimir_matriz (basetype * matriz,int N);
void imprimir_superior(basetype * matriz,int N);
void sumaLU(param * parametros);



basetype *L;
basetype *U;
basetype *C;
int N;
int CANT_THREADS;


pthread_barrier_t barrera; 		// Barrerra

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
	NT = (N*(N+1))>>1;

	CANT_THREADS = atoi(argv[2]);
	L=(basetype*)malloc(NT*sizeof(basetype));			// Reserva memoria para L
	U=(basetype*)malloc(NT*sizeof(basetype));			// Reserva memoria para U	
	C=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para C


	printf("Dimensión de la matriz: %d*%d \n",N,N);

	for (int i = 0; i < NT; ++i)
	{
		for (int j = 0; j < NT; ++j)
		{
			U[i+j*(j+1)/2]=rand()%5;
			L[i+N*j - i*(i+1)/2]=rand()%5;
		}
	}

	//imprimir_matriz(L,N);
	//printf("\n");
	imprimir_matriz(U,N);
	printf("\n");
	imprimir_superior(U,N);



	
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

	//imprimir_matriz(C,N);
	free(L);
	free(U);
	free(C);


}

double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

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

void imprimir_superior (basetype * matriz,int N){
	int i;
	int j;
	for (i=0;i<N;i++){
		for (j = 0 ; j < N; j++){
			if(i>=j)
				printf ("%.1f\t",matriz [ i + (j*(j+1)/2) ]);
			else
				printf ("0\t");
		}
		printf("\n");
	}
}

void *funcion_threads(void *arg) {
	param* parametro = (param*)arg;
	double tiempo_inicial2;
	//printf("Mi ID es: %d \n",(*parametro).id);
	if ((*parametro).id==0){
		tiempo_inicial2=dwalltime();
	}
	sumaLU(parametro);
	pthread_barrier_wait(&barrera); //espero a que todos los hilos finalicen
	if ((*parametro).id==0){
		printf("-- Fin de suma (pthread) -->> \t Tiempo: %f \n",dwalltime()-tiempo_inicial2);
		//printf("Matriz A:\n");
	}
	pthread_exit(NULL);
}

void sumaLU(param* parametro)
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
			C[i*N+j] = L[i*N+j]+U[i*N+j];
		}


	}
}