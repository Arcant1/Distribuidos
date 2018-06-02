#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

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


//Variables globales

// Matrices
basetype *A;
basetype *B;
basetype *C;
basetype *D;
basetype *E;
basetype *F;
basetype *L;
basetype *U;
basetype *LT;
basetype *UT;

basetype promB 	=0;
basetype promL 	=0;
basetype promU 	=0;
basetype prodLU =0;

//Matrices de resultados intermedios
basetype * sumaParcialB;
basetype * sumaL;
basetype * sumaU;
basetype * ULA;
basetype * AC;
basetype * bL;
basetype * BE;
basetype * bD;
basetype * UF;
basetype * ULLACbLBE;
basetype * bDUF;
basetype * bLBE;
basetype * ulAAC;
basetype * resultado;


int pos_max_global;

int N;
int NT;
int CANT_THREADS;

double speedup;
double eficiencia;
double tiempo_paral;
double tiempo_sec;

basetype factor_global;

pthread_barrier_t barrera; 		// Barrera

//Definición de funciones concurrentes

void funcion_threads (void *arg);
void promedioBOMP ();
void prod_escalarOMP (basetype * m1, basetype a, basetype * m2);
void suma_matrizOMP ( basetype * m1, basetype * m2, basetype * res);
void multiplicacion_omp ( basetype * m1, basetype * m2, basetype * m3, int dim);
void multiplicacionXTriangularUOMP ( basetype * m1, basetype * m2, basetype * m3, int dim);
void multiplicacionXTriangularLOMP ( basetype * m1, basetype * m2, basetype * m3, int dim);

void prodPromLUOMP ();
void imprimir_matriz (basetype * matriz,int N);
double dwalltime ();


//Funciones secuenciales

double tiempo_copia_total;

#ifdef COMPARAR_SECUENCIAL
	void multiplicacion_secuencial (basetype *A,basetype *B,basetype *C,int N);
	void multiplicacionXTriangularLSECUENCIAL (basetype * m1, basetype * m2, basetype *m3, int dim);
	void multiplicacionXTriangularUSECUENCIAL (basetype * m1, basetype * m2, basetype *m3, int dim);
	void promedioBSECUENCIAL ();
	void prod_escalarSECUENCIAL (basetype * m1, basetype a, basetype * m2);
	void prodPromLUSECUENCIAL ();

	void suma_matrizSECUENCIAL (basetype * m1, basetype * m2, basetype * res);
	void verificar_resultado (basetype *C,basetype *C_secuencial,int N);

	basetype *C_secuencial;	// Matriz resultado
#endif

int main(int argc,char *argv[])
{
	int i;
	int j;
	register unsigned long indice;
	double tiempo_inicial;
	if (argc != 3) {
		printf("Error en la cantidad de parámetros. Parametro 1 -> long matriz Parametro 2 -> n° threads");
		printf("Ingresar la dimension de la matriz!\n");
		return 0;
	}

	N = atoi(argv[1]);	//Dimensión de la matriz: N*N
	NT = N*N-(N*(N+1))/2;
	omp_set_num_threads(atoi(argv[2]));

	printf("Dimensión de la matriz: %d*%d \n",N,N);

	//Vector usado para guardar las sumas parciales de la matriz B para luego hacer el promedio
	sumaParcialB = (basetype*)malloc(sizeof(basetype)*CANT_THREADS);
	basetype promB;

	//Arreglos usados para calcular promedios
	sumaL = (basetype*)malloc(sizeof(basetype)*CANT_THREADS);
	sumaU = (basetype*)malloc(sizeof(basetype)*CANT_THREADS);

	// Reserva de memoria para las matrices
	A=(basetype*)malloc(sizeof(basetype)*N*N); // Reserva memoria para A
	B=(basetype*)malloc(N*N*sizeof(basetype)); // Reserva memoria para B
	C=(basetype*)malloc(N*N*sizeof(basetype)); // Reserva memoria para C
	D=(basetype*)malloc(N*N*sizeof(basetype)); // Reserva memoria para D
	E=(basetype*)malloc(N*N*sizeof(basetype)); // Reserva memoria para E
	F=(basetype*)malloc(N*N*sizeof(basetype)); // Reserva memoria para F

	//Caso especial de matrices triangulares
	L=(basetype*)malloc(N*N*sizeof(basetype)); // Reserva memoria para L
	U=(basetype*)malloc(N*N*sizeof(basetype)); // Reserva memoria para U

	LT=(basetype*)malloc(NT*sizeof(basetype)); // Reserva memoria para L transformada para ahorrar espacio
	UT=(basetype*)malloc(NT*sizeof(basetype)); // Reserva memoria para U transformada para ahorrar espacio

	AC 		=	(basetype*)malloc(sizeof(basetype)*N*N);
	ULA 	=	(basetype*)malloc(sizeof(basetype)*N*N);
	ulAAC 	=	(basetype*)malloc(sizeof(basetype)*N*N);
	BE 		=	(basetype*)malloc(sizeof(basetype)*N*N);
	bLBE	=	(basetype*)malloc(sizeof(basetype)*N*N);
	bD 		=	(basetype*)malloc(sizeof(basetype)*N*N);
	UF		=	(basetype*)malloc(sizeof(basetype)*N*N);
	bDUF	=	(basetype*)malloc(sizeof(basetype)*N*N);
	ULLACbLBE=	(basetype*)malloc(sizeof(basetype)*N*N);

	printf("Matrices creadas \n");



	//Inicialización ALEATORIA de las matrices

	//Inicializa matrices A y B
	# pragma omp parallel for private(i,j)
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j]=rand()%5;  // Inicializa matriz A con random
			B[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			C[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			D[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			E[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			F[i*N+j]=rand()%5; 	// Inicializa matriz B con random
		}
	}
	printf("Matrices inicializadas 1\n");

	//Inicializacion de matrices triangulares superior por columnas e inferior por filas
# pragma omp parallel for private(i,j)
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			if(i==j)
			{
				L[i*N+j]=rand()%5;
				U[i*N+j]=rand()%5;
			}
			else if(i>j)
			{
				U[i*N+j]=rand()%5;
				L[i*N+j]=0;
			}
			else
			{
				L[i*N+j]=rand()%5;
				U[i*N+j]=0;
			}
		}
	}
	printf("Matrices LU inicializadas \n");

	//Transformo las matrices triangulares en arreglos para ahorrar espacio y libero el espacio ocupado por las triangulares
# pragma omp parallel for private(i,j)

	for ( i = 0; i < N; ++i)
	{
		for ( j = 0; j < N; ++j)
		{
			indice = N * i + j - ((i *(i+1))/2);
			UT[indice]	= U[i*N + j];
			indice = N * j + i - ((j *(j+1))/2);
			//printf("%lu \n",indice);
			LT[indice] = L[i*N+j];
			//UT[i*NT + j - i*(i+1)/2]	=	U[i*N+j];
			//LT[j*N + i - j*(j+1)/2]	=	L[j*N+i];
		}
	}
	printf("Matrices inicializadas \n");


	tiempo_paral = dwalltime()-tiempo_inicial;
	printf("\nTiempo Total (OMP) : %f\n\n",dwalltime()-tiempo_inicial);

	/*free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	free(F);
	//free(L);
	//free(U);
	//free(LT);
	//free(UT);
	free(AC);
	free(ULA);
	free(ulAAC);
	free(BE);
	free(bLBE);
	free(bD);
	free(UF);
	free(bDUF);
	free(ULLACbLBE);*/

	A=(basetype*)   malloc(N*N*sizeof(basetype)); // Reserva memoria para A
	B=(basetype*)   malloc(N*N*sizeof(basetype)); // Reserva memoria para B
	C=(basetype*)   malloc(N*N*sizeof(basetype)); // Reserva memoria para C
	D=(basetype*)   malloc(N*N*sizeof(basetype)); // Reserva memoria para D
	E=(basetype*)   malloc(N*N*sizeof(basetype)); // Reserva memoria para E
	F=(basetype*)   malloc(N*N*sizeof(basetype)); // Reserva memoria para F
																		//Caso especial de matrices triangulares
	L=(basetype*)   malloc(N*N*sizeof(basetype)); // Reserva memoria para L
	U=(basetype*)   malloc(N*N*sizeof(basetype)); // Reserva memoria para U

	LT=(basetype*)  malloc(NT*sizeof(basetype)); // Reserva memoria para L transformada para ahorrar espacio
	UT=(basetype*)  malloc(NT*sizeof(basetype)); // Reserva memoria para U transformada para ahorrar espacio
	printf("Matrices creadas\n");




	AC 		=	(basetype*)malloc(sizeof(basetype)*N*N);
	ULA 	=	(basetype*)malloc(sizeof(basetype)*N*N);
	ulAAC 	=	(basetype*)malloc(sizeof(basetype)*N*N);
	BE 		=	(basetype*)malloc(sizeof(basetype)*N*N);
	bLBE	=	(basetype*)malloc(sizeof(basetype)*N*N);
	bD 		=	(basetype*)malloc(sizeof(basetype)*N*N);
	UF		=	(basetype*)malloc(sizeof(basetype)*N*N);
	bDUF	=	(basetype*)malloc(sizeof(basetype)*N*N);
	ULLACbLBE=	(basetype*)malloc(sizeof(basetype)*N*N);

	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			A[i*N+j]=rand()%5;  // Inicializa matriz A con random
			B[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			C[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			D[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			E[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			F[i*N+j]=rand()%5; 	// Inicializa matriz B con random
		}
	}
	printf("Matrices inicializadas 1\n");

	//Inicializacion de matrices triangulares superior por columnas e inferior por filas

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			if(i==j)
			{
				L[i*N+j]=rand()%5;
				U[i*N+j]=rand()%5;
			}
			else if(i>j)
			{
				U[i*N+j]=rand()%5;
				L[i*N+j]=0;
			}
			else
			{
				L[i*N+j]=rand()%5;
				U[i*N+j]=0;
			}
		}
	}
	printf("Matrices LU inicializadas \n");

	//Transformo las matrices triangulares en arreglos para ahorrar espacio y libero el espacio ocupado por las triangulares
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			indice = (N*i)+j-(i*(i+1)/2);
			UT[indice] = U[i*N+j];
			indice = (N*j)+i-(j*(j+1)/2);
			LT[indice] = L[i*N+j];

		}
	}
	printf("Matrices LT y UT inicializadas \n");

	printf("Matrices inicializadas\n");
	tiempo_inicial=dwalltime();

	printf("etapa 0\n");
	multiplicacion_secuencial(A,C,AC,N);

	printf("etapa 1\n");
	prodPromLUSECUENCIAL();

	printf("etapa 2\n");
	prod_escalarSECUENCIAL(A,prodLU,ULA);

	printf("etapa 3\n");
	multiplicacion_secuencial(ULA,AC,ulAAC,N);


	printf("etapa 4\n");
	promedioBSECUENCIAL();

	printf("etapa 5\n");
	multiplicacion_secuencial(B,E,BE,N);

	printf("etapa 6\n");
	multiplicacionXTriangularLSECUENCIAL(BE,LT,bLBE,N);


	printf("etapa 7\n");
	prod_escalarSECUENCIAL(bLBE,promB,bLBE);

	printf("etapa 8\n");
	prod_escalarSECUENCIAL(D,promB,bD);

	printf("etapa 9\n");
	multiplicacionXTriangularUSECUENCIAL(F,UT,UF,N);


	printf("etapa 10\n");
	multiplicacion_secuencial(bD,UF,bDUF,N);


	printf("etapa 11\n");
	suma_matrizSECUENCIAL(ulAAC,bLBE,ULLACbLBE);


	printf("etapa 12\n");
	resultado= (basetype*)malloc(sizeof(basetype)*N*N);
	suma_matrizSECUENCIAL(ULLACbLBE,bDUF,resultado);


	tiempo_sec = dwalltime()-tiempo_inicial;
	printf("-- Fin de multiplicacion (secuencial) -->> \t Tiempo: %f \n", tiempo_sec);

	speedup = tiempo_sec / tiempo_paral;
	printf("-- Speedup conseguido: %f \n", speedup);
	eficiencia = speedup / CANT_THREADS;
	printf("-- Eficiencia: %f \n", eficiencia);


	//Libera memoria

	free(resultado);
	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	free(F);
	free(L);
	free(U);
	free(LT);
	free(UT);
	free(AC);
	free(ULA);
	free(ulAAC);
	free(BE);
	free(bLBE);
	free(bD);
	free(UF);
	free(bDUF);
	free(ULLACbLBE);
	return(0);
}





// ------------------------
// -- FUNCIONES PTHREADS //
// ------------------------

	void funcion_OMP()
	{
		double tiempo_inicial2;
		
		tiempo_inicial2=dwalltime();
		
		multiplicacion_omp(A,C,AC,N);

		printf("etapa 0\n");
		

		prodPromLUOMP();

		printf("etapa 1\n");

		prod_escalarOMP(A,prodLU,ULA);

		printf("etapa 2\n");

		multiplicacion_omp(ULA,AC,ulAAC,N);

		printf("etapa 3\n");
		

		promedioBOMP();

		printf("etapa 4\n");

		multiplicacion_omp(B,E,BE,N);

		printf("etapa 5\n");

		multiplicacionXTriangularLOMP(BE,LT,bLBE,N);

		printf("etapa 6\n");
		

		prod_escalarOMP(bLBE,promB,bLBE);

		printf("etapa 7\n");
		

		prod_escalarOMP(D,promB,bD);


		printf("etapa 8\n");
		
		multiplicacionXTriangularUOMP(F,UT,UF,N);


		printf("etapa 9\n");

		
		multiplicacion_omp(bD,UF,bDUF,N);

		printf("etapa 10\n");
		

		suma_matrizOMP(ulAAC,bLBE,ULLACbLBE);

		printf("etapa 11\n");
		resultado= (basetype*)malloc(sizeof(basetype)*N*N);

		
		suma_matrizOMP(ULLACbLBE,bDUF,resultado);

		printf("etapa 12\n");
	
	
		printf("-- Fin de operacion (OMP) -->> \t Tiempo: %f \n",dwalltime()-tiempo_inicial2);
		//printf("Matriz A:\n");
}

/*
*	m2 es triangular superior
*/
void multiplicacionXTriangularUOMP(basetype * m1, basetype * m2, basetype *m3, int dim)
{
	basetype total;
	basetype aux;
	int i,j,k;
	# pragma omp parallel for private(i,j,k)
	for( i=0;i<=dim;i++)
	{
		//Recorre solo algunas filas
		for( j=0;j<dim;j++)
		{
			//Recorre todas las columnas
			total=0;
			for( k=0;k<dim;k++)
			{
				if(i>=j)
				{
					aux=m2[	j + k*(k+1)/2];
				}
				else aux = 0;
				total+=m1[i*dim+k]*aux;
			}
			m3[i*dim+j] = total;
		}
	}
	
}

void multiplicacionXTriangularOMP(basetype * m1, basetype * m2, basetype *m3, int dim)
{
	basetype total;
	basetype aux;
	int i,j,k;
# pragma omp parallel for private(i,j,k)
	for( i=0;i<dim;i++)
	{
		for( j=0;j<dim;j++)
		{
			total=0;
			for( k=0;k<dim;k++)
			{
				if(i<=j)
				{
					aux=m2[	j + k*(k+1)/2];
				}
				else aux = 0;
				total+=m1[i*dim+k]*aux;
			}
			m3[i*dim+j] = total;
		}
	}
}

void multiplicacionXTriangularLOMP( basetype * m1, basetype * m2, basetype *m3, int dim)
{
	basetype total;
	basetype aux;
	int i,j,k;
# pragma omp parallel for private(i,j,k)
	for( i=0;i<dim;i++)
	{
		for( j=0;j<dim;j++)
		{
			total=0;
			for( k=0;k<dim;k++)
			{
				if(i<=j)
				{
					aux=m2[	k + j*(j+1)/2];
				}
				else aux = 0;
				total+=m1[i*dim+k]*aux;
			}
			m3[i*dim+j] = total;
		}
	}
}


void multiplicacion_omp(basetype *A,basetype *B,basetype *C,int N)
{
	int i,j,k;
# pragma omp parallel for private(i,j,k)
	for(i=0;i<N;i++){ 
		for(j=0;j<N;j++){
			C[i*N+j]=0;
			for(k=0;k<N;k++){
				C[i*N+j]= C[i*N+j] + A[i*N+k]*B[k+j*N];
			}
		}
	}   
}

void prodPromLUOMP()
{
	int i;

	sumaL[omp_get_thread_num()]=0;
	sumaU[omp_get_thread_num()]=0;

# pragma omp parallel for private(i)
	for(i=0;i<=NT;i++)
	{
		//Recorre solo algunas filas
		sumaL[omp_get_thread_num()]+=UT[i];
		sumaU[omp_get_thread_num()]+=LT[i];

	}
	
	for (int i = 0; i < omp_get_num_thread(); ++i)
	{
		promL+=sumaL[i];
		promU+=sumaU[i];
	}
	promL/=omp_get_num_thread();
	promU/=omp_get_num_thread();
	prodLU=promU*promL;
    
}

basetype prodLUSEC;
basetype sumaLSEC,sumaUSEC;

void prodPromLUSECUENCIAL()
{
	sumaLSEC=0;
	sumaUSEC=0;
	for(int i=0;i<=NT;i++)
	{
			sumaLSEC+=UT[i];
			sumaUSEC+=LT[i];
	}

	prodLUSEC=promU*promL;

}

void promedioBOMP()
{
	basetype total;
	int i,j;
	//printf("ID: %d \t fila_inicial=%d \t fila_final=%d \n",id,fila_inicial,fila_final);
	total=0;
	// Multiplica A*B=C
# pragma omp parallel for private(i,j)
	for(i=0;i<=N;i++)
	{
		// Recorre solo algunas filas
		for(j=0;j<N;j++)
		// Recorre todas las columnas
		{
			total+=B[i*N + j];
		}
	}
	sumaParcialB[omp_get_thread_num()]=total;
   
	for (int i = 0; i < omp_get_num_thread(); ++i)
	{
		promB+=sumaParcialB[i];
	}
   	promB /=(N*N);
    
}

void promedioBSECUENCIAL()
{

	basetype total=0;
	// Multiplica A*B=C
	for(int i=0;i<=N;i++)
	{
		// Recorre solo algunas filas
		for(int j=0;j<N;j++)
		// Recorre todas las columnas
		{
			total+=B[i*N + j];
		}
	}

	promB =total/(N*N);

}

// --------------------------
// -- FUNCIONES AUXILIARES //
// --------------------------

void imprimir_matriz (basetype * matriz,int N){

	for (int i=0;i<N;i++){
		for (int j = 0 ; j < N ; j++){
			printf ("%.1f\t",matriz [ i * N + j ]);
		}
		printf("\n");
	}
}

void prod_escalarOMP (basetype * m1, basetype a, basetype * m2, int N)
{
	int i,j;
# pragma omp parallel for private(i,j)
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			m2[i*N+j]=m1[i*N+j]*a;
		}
	}
}

void prod_escalarSECUENCIAL ( basetype * m1, basetype a, basetype * m2)
{

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			m2[i*N+j]=m1[i*N+j]*a;
		}
	}
}

void suma_matrizOMP (basetype * m1, basetype * m2, basetype * res,int N)
{
	int i,j;
# pragma omp parallel for private(i,j)

	for ( i = 0; i < N; ++i)
	{
		for ( j = 0; j < N; ++j)
		{
			res[i*N+j]	=	m1[i*N+j]	+	m2[i*N+j];
		}
	}
}

void suma_matrizSECUENCIAL ( basetype * m1, basetype * m2, basetype * res)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			res[i*N+j]	=	m1[i*N+j]	+	m2[i*N+j];
		}
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

void multiplicacion_secuencial(basetype *A,basetype *B,basetype *C,int N){
	//printf("Comienzo etapa 1\n");

	basetype total;
	int i,j,k;

	// Multiplica A*B*D=C
	for(i=0;i<N;i++){
		for(j=0;j<N;j++)
		{
			total=0;
			for(k=0;k<N;k++)
			{
				total+=A[i*N+k]*B[k*N+j];	// total=A*B
			}
			C[i*N+j] = total;		// C=total
		}
	}
}

void verificar_resultado(basetype *C,basetype *C_secuencial,int N){
	int i,j;
	int check = 1;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			check=check&&(C[i*N+j]==C_secuencial[i*N+j]);
		}
	}
	if(check)
	{
		printf("Resultado correcto\n");
	}
	else
	{
		printf("Resultado erroneo \n");
	}
}
