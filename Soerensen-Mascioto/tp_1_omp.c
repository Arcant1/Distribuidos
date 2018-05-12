#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

// Tipo de los elementos en los vectores
// Compilar con -D_INT_ para vectores de tipo entero
// Compilar con -D_DOUBLE_ para vectores de tipo double
// Predeterminado vectores de tipo float

#define CANT_THREADS 3
#define COMPARAR_SECUENCIAL

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

int N;


result_parciales parciales [CANT_THREADS];	// Los threads guardan los valores calculados por cada uno de ellos (etapa 2)
max_parciales maximos [CANT_THREADS];		// Los threads guardan el máx encontrado por cada uno de ellos (ordenación --> etapa 3)


// -- Definición de funciones --

void etapa1();
void etapa2();
void etapa3();
void imprimir_matriz (basetype * matriz,int N);
double dwalltime();

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
	
	printf ("\033c");	// Limpia la consola
	if (argc != 2) {
		printf("Error en la cantidad de parámetros.");
		printf("1er parámetro: N --> dimensión de la matriz\n");
		return 0;
	}
	
	
	N = atoi(argv[1]);	// Dimensión de la matriz: N*N
		
	
		

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
	
	printf("Comienzo de las etapas. Dimensión de la matriz: %d*%d \n",N,N);

	// -- Comienzo de las etapas --
	omp_set_num_threads(CANT_THREADS);	// Setea el número de threads
	
	double tiempo_inicial=dwalltime();
	etapa1();	// C = A*B*D	
	printf("-- Fin de etapa 1 (omp) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial);
	//	printf("Matriz C:\n");
	//imprimir_matriz(C,N);

	double tiempo_inicial2=dwalltime();
	etapa2();	// C*factor
	printf("-- Fin de etapa 2 (omp) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial2);
	//printf("Matriz C:\n");
	//imprimir_matriz(C,N);
	
	double tiempo_inicial3=dwalltime();
	etapa3();	// Ordenación de columnas
	printf("-- Fin de etapa 3 (omp) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial3);
	//printf("Matriz C:\n");
	//imprimir_matriz(C,N);
	
	printf("\nTiempo Total (omp) : %f\n\n",dwalltime()-tiempo_inicial);

	#ifdef COMPARAR_SECUENCIAL
	tiempo_inicial=dwalltime();
	etapa1_secuencial(A,B,C_secuencial,D,N);	// C = A*B*D	
	printf("-- Fin de etapa 1 (secuencial) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial);

	tiempo_inicial2=dwalltime();	
	etapa2_secuencial(A,B,C_secuencial,N);		// C*factor
	printf("-- Fin de etapa 2 (secuencial) -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial2);

	tiempo_inicial3=dwalltime();
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


// ----------------------
// -- FUNCIONES OPENMP //
// ----------------------

void etapa1(){
	//printf("Comienzo etapa 1\n");

	// -- Multiplica A*B*D=C --

	// -- FORMA 1 --
	int i,j,k;
	basetype total;
	#pragma omp parallel for private(i,j,k,total)
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




void etapa2(){
	//printf("Comienzo etapa 2\n");
	int cant_filas = N/CANT_THREADS;
	int i,j,k;
	
	basetype total;
	basetype maxA = -1, minA = 999;  
	basetype maxB = -1, minB = 999;
	basetype totalA = 0;
	basetype totalB = 0;
	basetype actA, actB;

	#pragma omp parallel for  firstprivate(maxA,minA,maxB,minB,totalA,totalB) private(actA,actB,i,j)
	for(i=0;i<N;i++){ 	// Calcula valores para A y B
		int id = omp_get_thread_num();
		for( j=0;j<N;j++){
			// Cálculos para A
			actA = A [ i * N + j];
			if(actA<minA){
				parciales[id].minA=actA;		// Actualiza mínimo
				minA=actA;
			}
			if(actA>maxA){
				parciales[id].maxA=actA;	// Actualiza máximo	
				maxA=actA;	
			}			

			parciales[id].totalA+=actA;	// Incrementa total

			// Cálculos para B
			actB = B [ i * N + j];
			if(actB<minB){
				parciales[id].minB=actB;		// Actualiza mínimo
				minB=actB;
			}
			if(actB>maxB){
				parciales[id].maxB=actB;	// Actualiza máximo	
				maxB=actB;	
			}			
			parciales[id].totalB+=actB;	// Incrementa total
		}	
	}
	
	/*
	printf("maxA = %.4f \t maxB = %.4f \n",parciales[0].maxA,parciales[0].maxB);
	printf("minA = %.4f \t minB = %.4f \n",parciales[0].minA,parciales[0].minB);
	printf("totalA = %.4f \t totalB = %.4f \n",parciales[0].totalA,parciales[0].totalB);
	*/


	basetype maxA_total = -1, minA_total = 99999999;  
	basetype maxB_total = -1, minB_total = 99999999;
	basetype totalA_total = 0;
	basetype totalB_total = 0;
	basetype avgA, avgB;
	basetype factor;
	int m;

	// Calcula los valores totales
	for (m=0;m<CANT_THREADS;m++){
		if (parciales[m].maxA>maxA_total){
			maxA_total = 	parciales[m].maxA;
		}
		if (parciales[m].maxB>maxB_total){
			maxB_total = 	parciales[m].maxB;
		}
		if (parciales[m].minA<minA_total){
			minA_total = 	parciales[m].minA;
		}
		if (parciales[m].minB<minB_total){
			minB_total = 	parciales[m].minB;
		}
		totalA_total+=parciales[m].totalA;
		totalB_total+=parciales[m].totalB;
	

	}

	// Promedios
	avgA=totalA_total/(N*N);
	avgB=totalB_total/(N*N);

	factor=(((maxA_total-minA_total)*(maxA_total-minA_total))/avgA) * (((maxB_total-minB_total)*(maxB_total-minB_total))/avgB);	// Calcula factor

	/*
	printf("maxA = %.4f \t maxB = %.4f \n",maxA,maxB);
	printf("minA = %.4f \t minB = %.4f \n",minA,minB);
	printf("avgA = %.4f \t avgB = %.4f \n",avgA,avgB);
	printf("factor: %.4f\n", factor);
	*/

	
	// -- C*factor --
	#pragma omp parallel for private (i,j)
	for(i=0;i<=N;i++){
		for(j=0;j<N;j++){
			C[i*N+j]=C[i*N+j]*factor;
		}
	}
	
}


void etapa3(){
	//printf("Comienzo etapa 3\n");

	int j;
	

	// -- Ordenación por selección --
	
	basetype fila_actual[N];
	basetype fila_max [N];
	basetype valor_max = -1;
	int pos_max;

	int cant_filas_restantes= N;		// Cantidad de filas que restan ordenar
	int cant_filas_x_thread = N/CANT_THREADS;	// Cantidad de filas que procesa cada thread en cada iteración (va disminuyendo)

	int cant_columnas_restantes= N;		// Cantidad de filas que restan ordenar
	int cant_columnas_x_thread = N/CANT_THREADS;	// Cantidad de filas que procesa cada thread en cada iteración (va disminuyendo)

	for(j=0;j<N;j++){	// Recorre por columna
		#pragma omp parallel 
		{
			int i,k,m;
			int id = omp_get_thread_num();
			
		
			basetype actual;
			cant_filas_restantes = N;
			cant_filas_x_thread = N/CANT_THREADS;
			
			cant_columnas_restantes= N-j;
			cant_columnas_x_thread = cant_columnas_restantes/CANT_THREADS;

			int columna_inicial = j+id*cant_columnas_x_thread;
			int columna_final;

			if(id==CANT_THREADS-1){
				columna_final = N-1;		
			}
			else{
				columna_final = columna_inicial + cant_columnas_x_thread -1;
			}
			

			for(i=0;i<N;i++){	// Recorre por filas

				int fila_inicial = i+id*cant_filas_x_thread;
				int fila_final;
				
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
				#pragma omp barrier 	// Barrera (espera a que cada thread encuentre su máximo local)

				#pragma omp master	// Master calcula el maximo global (valor y posición)
				{	
					valor_max = -1;
					for (m=0;m<CANT_THREADS;m++){
						if (maximos[m].valor>valor_max){
							valor_max = maximos[m].valor;
							pos_max = maximos[m].posicion;
						}
					}
					//printf("Columna %d \t Maximo: %f\n",j,valor_max);
					cant_filas_restantes--;		// Decrementa la cantidad de filas restantes
					cant_filas_x_thread = cant_filas_restantes/CANT_THREADS;
				
				}
				
				#pragma omp barrier 	// Barrera (espera a que el thread 0 calcule el máximo global)

				//printf("Columna %d \t ID: %d \t fila_inicial: %d \t fila_final: %d\t Filas restantes %d \t Filas x thread: %d\n",j,id,fila_inicial,fila_final,cant_filas_restantes,cant_filas_x_thread);

				// Copia filas que luego swappea
				// Cada thread copia un conj de columnas (colu inicial=fila_inicial y col final=fila_final)
				for (k=columna_inicial;k<=columna_final;k++){	
					fila_actual [k] = C [i*N+k];
					fila_max [k] = C [pos_max*N+k]; 
				}

				// Swappea filas
				for (k=columna_inicial;k<=columna_final;k++){
					C [i*N+k] = fila_max[k];
					C [pos_max*N+k]= fila_actual[k]; 
				}

				#pragma omp barrier 
				#pragma omp master
				{
					//printf("Matriz C - Ordenación columna %d \n",j);
					//imprimir_matriz(C,N);
				}			
			}
			
		}

		//printf("Matriz C - Ordenación columna %d \n",j);
		//imprimir_matriz(C,N);
		
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

