#include<stdio.h>
#include<stdlib.h>
#include <sys/time.h>


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


// -- Definición de funciones --

void etapa1(basetype *A,basetype *B,basetype *C,basetype *D,int N);
void etapa2(basetype *A,basetype *B,basetype *C,int N);
void etapa3(basetype *C,int N);
void imprimir_matriz (basetype * matriz,int N);
double dwalltime();

int main(int argc,char *argv[]){
	
	int N;
	int i;
	int j;
	
	// Matrices
	basetype *A;	
	basetype *B;
	basetype *C;	// Matriz resultado
	basetype *D;	// Matriz diagonal

	printf ("\033c");	// Limpia la consola
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


	// -- Inicialización ALEATORIA de las matrices --
	
	// Inicializa matrices A y B
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			//A[i*N+j]=rand()%5;  	// Inicializa matriz A con random
			A[i*N+j]=1;  	// Inicializa matriz A con 1
			//B[i*N+j]=rand()%5; 	// Inicializa matriz B con random
			B[i*N+j]=1; 	// Inicializa matriz B con 1		
		}	
	}
	A[0]=2; B[0]=2;


	// Inicializa matriz D (diagonal)
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			if(i==j){
				// D[i*N+j]=rand()%5;  	// Inicializa matriz D	
				D[i*N+j]=1;  	// Inicializa matriz D	
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
	
	
	// -- Comienzo de las etapas --

	double tiempo_inicial=dwalltime();
	etapa1(A,B,C,D,N);	// C = A*B*D	
	printf("-- Fin de etapa 1 -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial);
	//printf("Matriz C:\n");
	//imprimir_matriz(C,N);

	double tiempo_inicial2=dwalltime();
	etapa2(A,B,C,N);	// C*factor
	printf("-- Fin de etapa 2 -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial2);
	//printf("Matriz C:\n");
	//imprimir_matriz(C,N);
	
	double tiempo_inicial3=dwalltime();
	etapa3(C,N);	// Ordenación de columnas
	printf("-- Fin de etapa 3 -- \t Tiempo: %f \n",dwalltime()-tiempo_inicial3);
	//printf("Matriz C:\n");
	//imprimir_matriz(C,N);
	
	
	printf("\nTiempo Total: %f\n\n",dwalltime()-tiempo_inicial);
	
	// Libera memoria
	free(A);
 	free(B);
 	free(C);
 	free(D);

	return(0);

}




void etapa1(basetype *A,basetype *B,basetype *C,basetype *D,int N){
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




void etapa2(basetype *A,basetype *B,basetype *C,int N){
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


void etapa3(basetype *C,int N){
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

