#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#ifdef _INT_
typedef int basetype;     // Tipo para elementos: int
#define tipo_dato   "ints"
#elif _DOUBLE_
typedef float basetype;  // Tipo para elementos: basetype
#define tipo_dato    "basetypes"
#else
typedef double basetype;   // Tipo para elementos: float     PREDETERMINADO
#define tipo_dato    "floats"
#endif

//Para calcular tiempo
double dwalltime(){
        basetype sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

void multiplicacionXTriangularUSECUENCIAL(basetype * m1, basetype * m2, basetype *m3, int N, int NT);

void multiplicacionXTriangularLSECUENCIAL( basetype * m1, basetype * m2, basetype *m3, int N,int NT);

void imprimir_matriz (basetype * matriz,int N);

basetype * prod_escalarSECUENCIAL ( basetype * m1, basetype a,int N);

basetype * multiplicacion_secuencial(basetype *A,basetype *B,int N);

basetype* suma_matrizSECUENCIAL ( basetype * m1, basetype * m2,int N);


int main(int argc, char** argv){

	int cantidadDeProcesos;
	int N; // Dimension de la matriz
	int sizeMatrix; // Cantidad total de datos matriz 	
	int sizePart,sizePartLU; // Cantidad de elementos por proceso
	//basetype *A_buf, *D_buf,*ab_temp,*de_temp,*abc_temp,*def_temp;
	//basetype *A,*B,*C,*D,*E,*F,*M;
	basetype *A_buf, *D_buf,*ab_temp,*LT_temp,*UT_temp,*LC_temp,*DU_temp;
	basetype *A,*B,*C,*D,*L,*U, *LT,*UT,*M;
	basetype uu=1.0,ll=1.0,escalar;
	basetype timetick;

	
	/***** PROGRAMA *****/
	if (argc < 2){
		printf("\n Falta un parametro ");
		printf("\n 1. Dimension de la matriz ");
		return 0;
	}

	N = atoi(argv[1]);
 	basetype NT = (N*(N+1))/2;
 	sizeMatrix=N*N; // Cantidad total de datos matriz
 	sizePartLU = NT;
 	
	
	ab_temp=(basetype*)malloc(sizeof(basetype)*sizeMatrix);
	LC_temp=(basetype*)malloc(sizeof(basetype)*sizeMatrix);
	DU_temp=(basetype*)malloc(sizeof(basetype)*sizeMatrix);
	
	 // El proceso con ID=0 inicializa y distribuye los datos
	
	A=(basetype*)malloc(sizeof(basetype)*sizeMatrix);
	B=(basetype*)malloc(sizeof(basetype)*sizeMatrix);
	C=(basetype*)malloc(sizeof(basetype)*sizeMatrix);
	L=(basetype*)malloc(sizeof(basetype)*sizeMatrix);
	U=(basetype*)malloc(sizeof(basetype)*sizeMatrix);
	LT=(basetype*)malloc(sizeof(basetype)*NT);
	UT=(basetype*)malloc(sizeof(basetype)*NT);
	D=(basetype*)malloc(sizeof(basetype)*sizeMatrix);
	M=(basetype*)malloc(sizeof(basetype)*sizeMatrix);

	for(int i=0 ;i<sizeMatrix;i++) { // Inicializa las matrices A y B en 1, el resultado sera una matriz con todos sus valores en N
		A[i]=1.0;
		B[i]=1.0;
		C[i]=1.0;
		D[i]=1.0;
		L[i]=1.0;
		U[i]=1.0;
  	} //end i

  	int indice;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			indice = N * i + j - ((i *(i+1))/2);
			UT[indice]	= U[i*N + j];

			indice = N * j + i - ((j *(j+1))/2);
			LT[indice] = L[i*N+j];			
		}
	}
	printf("Fin Inicializacion\n");

	

	timetick = dwalltime();

	basetype t1,t2,t3,t4,t5,t6;

	//*********Promedio L**********

	// Buffer para la Matriz LT
  	// Distribuye los elementos de la matriz L
	// Promedio de B
	basetype temp=0;

	for(int i=0;i<sizePartLU;i++)
		temp+=LT[i];

	temp/=sizeMatrix;


	
	//////////////////////////////////////////////////////////////////////

	//*********Promedio U**********

	// Buffer para la Matriz UT
  	// Distribuye los elementos de la matriz L
	// Promedio de B
	temp=0;

	for(int i=0;i<sizePartLU;i++)
		temp+=UT[i];

	temp/=sizeMatrix;


	
	///////////////////////////////////////////////////////////////////////////////



	escalar = uu*ll;

	

	t1=dwalltime();

	ab_temp= multiplicacion_secuencial(A,B,N);				//AB
	printf("Multiplicacion 1 tarda = %f  \n", dwalltime() - t1);


	//multiplicacionXTriangularUSECUENCIAL(D,UT,DU_temp,N,NT);	//DU
	t2=dwalltime();
	DU_temp= multiplicacion_secuencial(D,U,N);				
	printf("Multiplicacion 2 tarda = %f  \n", dwalltime() - t2);

	//multiplicacionXTriangularLSECUENCIAL(C,LT,LC_temp,N,NT);	//CL
	t3=dwalltime();
	LC_temp=multiplicacion_secuencial(C,L,N);				
	printf("Multiplicacion 3 tarda = %f  \n", dwalltime() - t3);

	t4=dwalltime();
	ab_temp = suma_matrizSECUENCIAL(ab_temp,LC_temp,N);
	printf("sum1 1 tarda = %f  \n", dwalltime() - t4);

	t5=dwalltime();
	ab_temp = suma_matrizSECUENCIAL(ab_temp,DU_temp,N);
	printf("Suma 2 tarda = %f  \n", dwalltime() - t5);

	t6=dwalltime();
	ab_temp= prod_escalarSECUENCIAL(ab_temp,escalar,N);
	printf("Prod escalar tarda = %f  \n", dwalltime() - t6);

											//

	printf("Tiempo en segundos %f \n", dwalltime() - timetick);

	//***** FIN PROGRAMA ****
	//imprimir_matriz(LC_temp,N);
	
	free(M);
	free(A);
	free(B);
	free(C);
	free(L);
	free(U);
	free(LT);
	free(UT);
	free(D);
	
	return(0);
}




void multiplicacionXTriangularUSECUENCIAL(basetype * m1, basetype * m2, basetype *m3, int N,int NT)
{
	basetype aux;
	basetype total;
	for(int i=0;i<N;i++)
	{
		//Recorre solo algunas filas
		for(int j=0;j<N;j++)
		{
			//Recorre todas las columnas
			total=0;
			for(int k=0;k<NT;k++)
			{
				if(i>=j)
				{
					aux=m2[	j + k*(k+1)/2];

				}
				else aux = 0;
				total+=m1[i*N+k]*aux;
			}
			m3[i*N+j] = total;
		}
	}
}

void multiplicacionXTriangularLSECUENCIAL( basetype * m1, basetype * m2, basetype *m3, int N,int NT)
{
	basetype total;
	basetype aux;
	for(int i=0;i<N;i++)
	{
		// Recorre solo algunas filas
		for(int j=0;j<N;j++)
		{
			// Recorre todas las columnas
			total=0;
			for(int k=0;k<NT;k++)
			{
				if(i<=j)
				{
					aux=m2[	k + j*(j+1)/2];
				}
				else aux = 0;
				total+=m1[i*N+k]*aux;
			}
			m3[i*N+j] = total;
		}
	}
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


basetype * prod_escalarSECUENCIAL ( basetype * m1, basetype a,int N)
{
	basetype * C=(basetype*)malloc(sizeof(basetype)*N*N);

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			C[i*N+j]=m1[i*N+j]*a;
		}
	}
	return C;
}



basetype * multiplicacion_secuencial(basetype *A,basetype *B,int N)
{
	//printf("Comienzo etapa 1\n");
	basetype * C=(basetype*)malloc(sizeof(basetype)*N*N);
	basetype total;
	int i,j,k;

	// Multiplica A*B*D=C
	for(i=0;i<N;i++){
		for(j=0;j<N;j++)
		{
			total=0;
			for(k=0;k<N;k++)
			{
				total+=A[i*N+k]*B[k+N*j];	// total=A*B
			}
			C[i*N+j] = total;		// C=total
		}
	}
	return C;
}

basetype * suma_matrizSECUENCIAL ( basetype * m1, basetype * m2,int N)
{
	basetype * res = (basetype*)malloc(sizeof(basetype)*N*N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			res[i*N+j]	=	m1[i*N+j]	+	m2[i*N+j];
		}
	}
	return res;
}
