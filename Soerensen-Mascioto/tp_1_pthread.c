#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

#define CANT_THREADS 4
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

// Tipo de dato para que los threads calculen sus resultados parciales
typedef struct result_parciales
{
	basetype maxA;
	basetype minA;
	basetype maxB;
	basetype minB;
	basetype totalA;
	basetype totalB;
}result_parciales;

typedef struct max_parciales
{
	basetype valor;
	int posicion;
}max_parciales;

// -- Variables globales -- 

// Matrices
basetype *A;	
basetype *B;
basetype *C;	// Matriz resultado
basetype *D;	// Matriz diagonal

// Arreglo donde los pthreads guardan los valores calculados por cada uno de ellos
result_parciales parciales [CANT_THREADS];

max_parciales maximos [CANT_THREADS];	// Los threads guardan el máx encontrado por cada uno de ellos (ordenación --> etapa 3)

int cant_filas_restantes;	// Cantidad de filas que restan ordenar
int cant_filas_x_thread;	// Cantidad de filas que procesa cada thread en cada iteración (va disminuyendo)
int pos_max_global;

int N;

basetype factor_global;

pthread_barrier_t barrera; 		// Barrerra


// -- Definición de funciones --

void *funcion_threads(void *arg);
void etapa1(param *parametro);
void etapa2(param *parametro);
void etapa3(param* parametro);
void imprimir_matriz (basetype * matriz,int N);
double dwalltime();
double tiempo_copia_total=0;

#ifdef COMPARAR_SECUENCIAL
void etapa1_secuencial(basetype *A,basetype *B,basetype *C_secuencial,basetype *D,int N);
void etapa2_secuencial(basetype *A,basetype *B,basetype *C_secuencial,int N);
void etapa3_secuencial(basetype *C_secuencial,int N);

void verificar_resultado(basetype *C,basetype *C_secuencial,int N);

basetype *C_secuencial;	// Matriz resultado
#endif

int main(int argc,char *argv[]){
	int i;
	int j;
	
	

	if (argc != 2) {
		printf("Error en la cantidad de parámetros.");
		printf("1er parámetro: N --> dimensión de la matriz\n");
		return 0;
	}
	
	
	N = atoi(argv[1]);	// Dimensión de la matriz: N*N
		
	printf("Dimensión de la matriz: %d*%d \n",N,N);
		
	
	// Reserva de memoria para las matrices 
	A=(basetype*)malloc(sizeof(basetype)*N*N);			// Reserva memoria para A
	B=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para B
	C=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para C
	D=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para D
	
	#ifdef COMPARAR_SECUENCIAL
	C_secuencial=(basetype*)malloc(N*N*sizeof(basetype));			// Reserva memoria para C
	#endif

	
	// -- Inicialización ALEATORIA de las matrices --
	
	// Inicializa matrices A y B
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j]=rand()%5;  	// Inicializa matriz A con random
			//A[i*N+j]=1;  	// Inicializa matriz A con random
			B[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			//B[i*N+j]=1; 	// Inicializa matriz B con random		
		}	
	}


	// Inicializa matriz D (diagonal)
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			if(i==j){
				D[i*N+j]=rand()%5;  	// Inicializa matriz D	
			}
			else{
				D[i*N+j]=0;
			}
		}	
	}
	

	// -- Inicialización con valores CONOCIDOS de las matrices --
	/*
	N=3;	

	// Inicialización de A
	A[0]=1;	A[1]=2;	A[2]=6;
	A[3]=5;	A[4]=9;	A[5]=3;
	A[6]=4;	A[7]=8;	A[8]=7;

	// Inicialización de B
	B[0]=3;	B[1]=2;	B[2]=1;
	B[3]=3;	B[4]=5;	B[5]=3;
	B[6]=1;	B[7]=7;	B[8]=2;   
	
	// Inicialización de C
	D[0]=4;	D[1]=0;	D[2]=0;
	D[3]=0;	D[4]=5;	D[5]=0;
	D[6]=0;	D[7]=0;	D[8]=3;   
	*/

	/*
	printf("Matriz A:\n");
	imprimir_matriz(A,N);
	printf("Matriz B:\n");
	imprimir_matriz(B,N);
	printf("Matriz D:\n");
	imprimir_matriz(D,N);
	*/

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
	for(i = 0; i < CANT_THREADS; i++) {
        	pthread_join(threads[i], NULL);
        	
    	}
    	
    	printf("\nTiempo Total (pthreads) : %f\n\n",dwalltime()-tiempo_inicial);
    	
    	#ifdef COMPARAR_SECUENCIAL
	tiempo_inicial=dwalltime();
	etapa1_secuencial(A,B,C_secuencial,D,N);	// C = A*B*D	
	printf("-- Fin de etapa 1 (secuencial) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial);

	double tiempo_inicial2=dwalltime();	
	etapa2_secuencial(A,B,C_secuencial,N);		// C*factor
	printf("-- Fin de etapa 2 (secuencial) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial2);

	double tiempo_inicial3=dwalltime();
	etapa3_secuencial(C_secuencial,N);		// Ordenación de columnas
	printf("-- Fin de etapa 3 (secuencial) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial3);
	printf("\nTiempo Total (secuencial) : %f\n\n",dwalltime()-tiempo_inicial);
	verificar_resultado(C,C_secuencial,N);
	
	free(C_secuencial);
	#endif
	
	
	// Libera memoria
	free(A);
 	free(B);
 	free(C);
 	free(D);

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
	etapa1(parametro);
	pthread_barrier_wait(&barrera);
	if ((*parametro).id==0){
		printf("-- Fin de etapa 1 (pthread) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial2);
	}
	//	printf("Matriz C:\n");
	
	
	
	if ((*parametro).id==0){
		tiempo_inicial2=dwalltime();
	}
	etapa2(parametro);
	pthread_barrier_wait(&barrera);
	if ((*parametro).id==0){
		printf("-- Fin de etapa 2 (pthread) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial2);
	} 
	
	
	
	if ((*parametro).id==0){
		tiempo_inicial2=dwalltime();
	}
	etapa3(parametro);
	
	if ((*parametro).id==0){
		printf("-- Fin de etapa 3 (pthread) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial2);
	} 	
}

void etapa1(param* parametro){
	//printf("Comienzo etapa 1\n");
	int id = (*parametro).id;
	basetype total;
	int i,j,k;


	// Filas que multiplica el thread
	int cant_filas = N/CANT_THREADS;	// Cant de filas que multiplica cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	//printf("ID: %d \t fila_inicial=%d \t fila_final=%d \n",id,fila_inicial,fila_final);

	// Multiplica A*B*D=C
	for(i=fila_inicial;i<=fila_final;i++){	// Recorre solo algunas filas
			for(j=0;j<N;j++){	// Recorre todas las columnas
				total=0;
				for(k=0;k<N;k++){
					total+=A[i*N+k]*B[k*N+j];	// total=A*B
				}
				// D tiene una única fila llena por c/ columna
				C[i*N+j] = total * D[j*N+j];		// C=total*D  
			}
     } 

}




void etapa2(param* parametro){
	//printf("Comienzo etapa 2\n");
	int id = (*parametro).id;
	basetype maxA = -1, minA = 99999999;  
	basetype maxB = -1, minB = 99999999;
	basetype totalA = 0;
	basetype totalB = 0;
	basetype avgA, avgB;
	basetype actA, actB;

	// Filas con las que trabaja el thread
	int cant_filas = N/CANT_THREADS;	
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;
	float factor;

	int i,j;  

	
	// Calcula valores para A y B
	for(i=fila_inicial;i<=fila_final;i++){	// Recorre solo algunas filas
		for( j=0;j<N;j++){		// Recorre todas las columnas
			// Cálculos para A
			actA = A [ i * N + j];
			if(actA<minA){
				minA=actA;		// Actualiza mínimo
			}
			if(actA>maxA){
				maxA=actA;	// Actualiza máximo		
			}			

			totalA+=actA;	// Incrementa total

			// Cálculos para B
			actB = B [ i * N + j];
			if(actB<minB){
				minB=actB;		// Actualiza mínimo
			}
			if(actB>maxB){
					maxB=actB;	// Actualiza máximo		
			}			

			
			totalB+=actB;	// Incrementa total
		}	
	}
	
	// Actualiza variables globales con resultados parciales
	parciales[id].maxA=maxA;
	parciales[id].maxB=maxB;
	parciales[id].minA=minA;
	parciales[id].minB=minB;
	parciales[id].totalA=totalA;
	parciales[id].totalB=totalB;
	pthread_barrier_wait(&barrera);	// Barrera

	// El thread de ID=0 calcula los valores totales
	if (id==0){
		// Resetea variables locales para calcular los totales
		maxA = -1; 
		minA = 99999999;  
		maxB = -1;
		minB = 99999999;
		totalA = 0;
		totalB = 0;
		// Calcula los valores totales
		for (i=0;i<CANT_THREADS;i++){
			if (parciales[i].maxA>maxA){
				maxA = 	parciales[i].maxA;
			}
			if (parciales[i].maxB>maxB){
				maxB = 	parciales[i].maxB;
			}
			if (parciales[i].minA<minA){
				minA = 	parciales[i].minA;
			}
			if (parciales[i].minB<minB){
				minB = 	parciales[i].minB;
			}
			totalA+=parciales[i].totalA;
			totalB+=parciales[i].totalB;
		

		}

		// Promedios
		avgA=totalA/(N*N);
		avgB=totalB/(N*N);
	

	
	
		factor=(((maxA-minA)*(maxA-minA))/avgA) * (((maxB-minB)*(maxB-minB))/avgB);	// Calcula factor
		factor_global = factor;
		/*
		printf("maxA = %.4f \t maxB = %.4f \n",maxA,maxB);
		printf("minA = %.4f \t minB = %.4f \n",minA,minB);
		printf("avgA = %.4f \t avgB = %.4f \n",avgA,avgB);
		printf("factor: %.4f\n", factor);
		*/	
 	}

	
	pthread_barrier_wait(&barrera);	// Barrera


	// C*factor
	for(i=fila_inicial;i<=fila_final;i++){	// Recorre solo algunas filas
		for(j=0;j<N;j++){
			C[i*N+j]=C[i*N+j]*factor_global;
		}
	}

}


void etapa3(param* parametro){
	//printf("Comienzo etapa 3\n");
	int i,j,k,m;
	int id = (*parametro).id;

	basetype fila_actual[N];
	basetype fila_max [N];
	basetype valor_max = -1;
	basetype actual;

	cant_filas_restantes= N;		// Cantidad de filas que restan ordenar
	cant_filas_x_thread = N/CANT_THREADS;	// Cantidad de filas que procesa cada thread en cada iteración (va disminuyendo)

	int cant_columnas_restantes= N;		// Cantidad de filas que restan ordenar
	int cant_columnas_x_thread = N/CANT_THREADS;	// Cantidad de filas que procesa cada thread en cada iteración (va disminuyendo)

	int columna_inicial, columna_final;
	int fila_inicial, fila_final;

	for(j=0;j<N;j++){	// Recorre por columna
	
		basetype actual;
		cant_filas_restantes = N;
		cant_filas_x_thread = N/CANT_THREADS;
		
		cant_columnas_restantes= N-j;
		cant_columnas_x_thread = cant_columnas_restantes/CANT_THREADS;

		columna_inicial = j+id*cant_columnas_x_thread;

		if(id==CANT_THREADS-1){
			columna_final = N-1;		
		}
		else{
			columna_final = columna_inicial + cant_columnas_x_thread -1;
		}

		for(i=0;i<N;i++){	// Recorre por filas

			fila_inicial = i+id*cant_filas_x_thread;

			if(id==CANT_THREADS-1){
				fila_final = N-1;	
			}
			else{
				fila_final = fila_inicial + cant_filas_x_thread -1;
			}

			if (fila_inicial<=fila_final){

			//printf("Columna %d \t ID: %d \t fila_inicial: %d \t fila_final: %d\t Filas restantes %d \t Filas x thread: %d\n",j,id,fila_inicial,fila_final,cant_filas_restantes,cant_filas_x_thread);

			actual = C[fila_inicial*N+j];	// En c/ iteración "actual" es un nuevo elemento de la col (cambia de fila)
	
			maximos[id].posicion = fila_inicial;
			maximos[id].valor = actual;
			//printf("Columna %d \t Actual: %f\n",j,maximos[id].valor);

			for (k=fila_inicial+1;k<=fila_final;k++){	// Recorre las pos restantes de la col para buscar el máx
				if (C[k*N+j]>maximos[id].valor){
					maximos[id].posicion = k;	// Actualiza posición del máx
					maximos[id].valor = C[k*N+j];	// Actualiza valor del máx			
				}
			}
		}
		else{
			maximos[id].valor = -1;		// Para que no queden valores viejos
		}
			
		pthread_barrier_wait(&barrera);	// Barrera (espera a que cada thread encuentre su máximo local)

		if (id==0){	
			valor_max = -1;
			for (m=0;m<CANT_THREADS;m++){
				if (maximos[m].valor>valor_max){
					valor_max = maximos[m].valor;
					pos_max_global = maximos[m].posicion;
				}
			}
			//printf("Columna %d \t Maximo: %f\n",j,valor_max);
			cant_filas_restantes--;		// Decrementa la cantidad de filas restantes
			cant_filas_x_thread = cant_filas_restantes/CANT_THREADS;
		
		}
	
		pthread_barrier_wait(&barrera);	// Barrera (espera a que el thread 0 calcule el máximo global)

		//printf("Columna %d \t ID: %d \t fila_inicial: %d \t fila_final: %d\t Filas restantes %d \t Filas x thread: %d\n",j,id,fila_inicial,fila_final,cant_filas_restantes,cant_filas_x_thread);
		
		
		int pos1,pos2;
		int pos_max_global_local = pos_max_global;
		basetype aux;

		// Copia filas que luego swappea
		// Cada thread copia un conj de columnas (colu inicial=fila_inicial y col final=fila_final)
		for (k=columna_inicial;k<=columna_final;k++){	
			//printf("Columna: %d \t ID: %d \t K:%d \t Col inicial:%d \t Col final: %d\n",j,id,k,columna_inicial,columna_final);
			pos1=i*N+k;
			pos2=pos_max_global*N+k;
			fila_actual [k] = C [pos1];
			fila_max [k] = C [pos2]; 
		}	
		
		
		// Swappea filas
		for (k=columna_inicial;k<=columna_final;k++){
			pos1=i*N+k;
			pos2=pos_max_global*N+k;
			C [pos1] = fila_max[k];
			C [pos2]= fila_actual[k]; 
		}
		

		/*
		for (k=columna_inicial;k<=columna_final;k++){
			pos1=i*N+k;
			pos2=pos_max_global*N+k;
			aux = C [pos1];
			C [pos1] = C [pos2]; 
			C [pos2] = aux;
		}
		*/
		

		pthread_barrier_wait(&barrera);	// Barrera (espera a que el thread 0 calcule el máximo global)
		
		}	// Fin recorrido por fila


		

	}	// Fin recorrido por columna
			


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
void etapa1_secuencial(basetype *A,basetype *B,basetype *C,basetype *D,int N){
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
				// D tiene una única fila llena por c/ columna
				C[i*N+j] = total * D[j*N+j];		// C=total*D  
			}
     	} 

}




void etapa2_secuencial(basetype *A,basetype *B,basetype *C,int N){
	//printf("Comienzo etapa 2\n");
	basetype maxA = -1, minA = 99999999;  
	basetype maxB = -1, minB = 99999999;
	basetype totalA = 0;
	basetype totalB = 0;
	basetype avgA, avgB;
	basetype actA, actB;

	int i,j;   
	
	// Calcula valores para A y B
	for(i=0;i<N;i++){
		for( j=0;j<N;j++){
			// Cálculos para A
			actA = A [ i * N + j];
			if(actA<minA){
				minA=actA;		// Actualiza mínimo
			}
			if(actA>maxA){
				maxA=actA;	// Actualiza máximo	
			}			
			totalA+=actA;	// Incrementa total

			// Cálculos para B
			actB = B [ i * N + j];
			if(actB<minB){
				minB=actB;		// Actualiza mínimo
			}
			if(actB>maxB){
				maxB=actB;	// Actualiza máximo		
			}			

			totalB+=actB;	// Incrementa total
		}	
	}
	
	
	// Promedios
	avgA=totalA/(N*N);
	avgB=totalB/(N*N);

	float factor=(((maxA-minA)*(maxA-minA))/avgA) * (((maxB-minB)*(maxB-minB))/avgB);	// Calcula factor

	/*
	printf("maxA = %.4f \t maxB = %.4f \n",maxA,maxB);
	printf("minA = %.4f \t minB = %.4f \n",minA,minB);
	printf("avgA = %.4f \t avgB = %.4f \n",avgA,avgB);
	printf("factor: %.4f\n", factor);
	*/

	// C*factor
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			C[i*N+j]=C[i*N+j]*factor;
		}
	}
}


void etapa3_secuencial(basetype *C,int N){
	//printf("Comienzo etapa 3\n");

	int i,j,k;


	// -- Ordenación por selección --
	int pos_max;
	basetype valor_max;
	basetype actual;


	basetype fila_actual[N];
	basetype fila_max [N];

	for(j=0;j<N;j++){	// Recorre por columna

		for(i=0;i<N;i++){	// Recorre por filas
			actual = C[i*N+j];	// En cada iteración "actual" toma el valor de un nuevo elemento de la columna (cambia de fila)
			pos_max = i;
			valor_max = actual;

			for (k=i+1;k<N;k++){	// Recorre todas las posiciones restantes de la columna para buscar el máximo
				if (C[k*N+j]>valor_max){
					pos_max = k;	// Actualiza posición del máximo
					valor_max = C[k*N+j];	// Actualiza valor del máximo			
				}
			}
			
			// Copia filas que luego swappea
			for (k=j;k<N;k++){
				fila_actual [k] = C [i*N+k];
				fila_max [k] = C [pos_max*N+k]; 
			}

			// Swappea filas
			for (k=j;k<N;k++){
				C [i*N+k] = fila_max[k];
				C [pos_max*N+k]= fila_actual[k]; 
			}

		}
		//printf("Matriz C - Ordenación columna %d \n",j);
		//imprimir_matriz(C,N);
		
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


