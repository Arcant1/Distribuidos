#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

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

//Para calcular tiempo
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

int main(int argc, char** argv){

	int miID; 
	int cantidadDeProcesos;
	int N; // Dimension de la matriz
	int sizeMatrix; // Cantidad total de datos matriz 	
	int sizePart; // Cantidad de elementos por proceso
	//double *A_buf, *D_buf,*ab_temp,*de_temp,*abc_temp,*def_temp;
	//double *A,*B,*C,*D,*E,*F,*M;
	double *A_buf, *D_buf,*AB_temp,*LC_temp,*DU_temp;
	double *A,*B,*C,*D,*L,*U;
	double uu=0.0,ll=0.0;
	double timetick;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&miID);
	MPI_Comm_size(MPI_COMM_WORLD,&cantidadDeProcesos);

	/***** PROGRAMA *****/
	if (argc < 2){
		printf("\n Falta un parametro ");
		printf("\n 1. Dimension de la matriz ");
		return 0;
	}

	N = atoi(argv[1]);
 	sizeMatrix=N*N; // Cantidad total de datos matriz
 	sizePart=sizeMatrix/cantidadDeProcesos; // Cantidad de elementos por bloque x Cantidad de bloques por proceso 	
 	NT = (N*(N+1))/2;
	
	A_buf=(double*)malloc(sizeof(double)*sizePart);
	B=(double*)malloc(sizeof(double)*sizeMatrix);
	C=(double*)malloc(sizeof(double)*sizeMatrix);
	D=(double*)malloc(sizeof(double)*sizePart);
	L=(double*)malloc(sizeof(double)*sizeMatrix);
	U=(double*)malloc(sizeof(double)*sizeMatrix);
	LT=(double*)malloc(sizeof(double)*sizeMatrix);
	UT=(double*)malloc(sizeof(double)*sizeMatrix);

	AB_temp=(double*)malloc(sizeof(double)*sizePart);
	LC_temp=(double*)malloc(sizeof(double)*sizePart);
	DU_temp=(double*)malloc(sizeof(double)*sizePart);
	
	
	B=(double*)malloc(sizeof(double)*sizeMatrix);
	C=(double*)malloc(sizeof(double)*sizeMatrix);
	D=(double*)malloc(sizeof(double)*sizePart);
	E=(double*)malloc(sizeof(double)*sizeMatrix);
	F=(double*)malloc(sizeof(double)*sizeMatrix);

	if(miID==0) { // El proceso con ID=0 inicializa y distribuye los datos
	
		A=(double*)malloc(sizeof(double)*sizeMatrix);
		D=(double*)malloc(sizeof(double)*sizeMatrix);
		M=(double*)malloc(sizeof(double)*sizeMatrix);

		for(int i=0 ;i<sizeMatrix;i++) { // Inicializa las matrices A y B en 1, el resultado sera una matriz con todos sus valores en N
			A[i]=1.0;
			B[i]=1.0;
			C[i]=1.0;
			D[i]=1.0;
			L[i]=1.0;
			U[i]=1.0;
	  	} //end i

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
	}

	timetick = dwalltime();

	//*********Promedio L**********

	// Buffer para la Matriz LT
  	// Distribuye los elementos de la matriz L
  	MPI_Scatter(B,sizePart,MPI_DOUBLE,A_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);
	// Promedio de B
	double temp=0;

	for(int i=0;i<sizePart;i++)
		temp+=A_buf[i];

	temp/=sizeMatrix;

	MPI_Allreduce(&temp, &b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 		//MPI_SUM supongo que suma todo lo que recibe en temp

	printf("Proceso %d, Promedio b = %lf \n",miID,b);

	//**********Promedio D**********

  	// Distribuye los elementos de la matriz D
  	MPI_Scatter(D,sizePart,MPI_DOUBLE,A_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);
	// Promedio de D
	temp=0;

	for(int i=0;i<sizePart;i++)
		temp+=A_buf[i];

	temp/=sizeMatrix;

	MPI_Allreduce(&temp, &d, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	printf("Proceso %d, Promedio d = %lf \n",miID,d);

	MPI_Scatter(A,sizePart,MPI_DOUBLE,A_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(B, sizeMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(C, sizeMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(D,sizePart,MPI_DOUBLE,D_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(E, sizeMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(F, sizeMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
												//A*B
	for(int i=0;i<sizePart/N;i++){
		for(int j=0;j<N;j++){
			ab_temp[i*N+j]=0;
			for(int k=0;k<N;k++){
				ab_temp[i*N+j]= ab_temp[i*N+j] + A_buf[i*N+k]*B[k+j*N];
			}
		}
	}
												//(A*B)*C
	for(int i=0;i<sizePart/N;i++){
		for(int j=0;j<N;j++){
			abc_temp[i*N+j]=0;
			for(int k=0;k<N;k++){
				abc_temp[i*N+j]= abc_temp[i*N+j] + ab_temp[i*N+k]*C[k+j*N];
			}
		}
	}
												//dABC
	for (int i=0;i<sizePart;i++){
		abc_temp[i]*=d;
	}

												//DE
	for(int i=0;i<sizePart/N;i++){
		for(int j=0;j<N;j++){
			de_temp[i*N+j]=0;
			for(int k=0;k<N;k++){
				de_temp[i*N+j]= de_temp[i*N+j] + D_buf[i*N+k]*E[k+j*N];
			}
		}
	}
												//(DE)F
	for(int i=0;i<sizePart/N;i++){
		for(int j=0;j<N;j++){
			def_temp[i*N+j]=0;
			for(int k=0;k<N;k++){
				def_temp[i*N+j]= def_temp[i*N+j] + de_temp[i*N+k]*F[k+j*N];
			}
		}
	}

												//b*DEF => DEF
	for (int i=0;i<sizePart;i++){
		def_temp[i]*=b;
	}

												//dABC+bDEF =>dABC
	for (int i=0;i<sizePart;i++){
		abc_temp[i]+=def_temp[i];
	}

												//
	MPI_Gather(abc_temp, sizePart, MPI_DOUBLE, M, sizePart, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	printf("Tiempo en segundos %f \n", dwalltime() - timetick);

	/*if(miID==0){
		for (int j=0; j<N; j++){
			for (int i=0; i<N; i++){
				printf ("%lf \t",M[j*N+i]);
			}
			printf ("\n");
		}
	}*/

	//***** FIN PROGRAMA ****
	
	free(B);
	free(C);
	free(E);
	free(F);
	free(ab_temp);
	free(abc_temp);
	free(de_temp);
	free(def_temp);
	free(A_buf);
	free(D_buf);
	if(miID==0){
		free(M);
		free(A);
		free(D);
	}

	MPI_Finalize();

	return(0);
}



void multiplicacionXTriangularU(param* parametro, basetype * m1, basetype * m2, basetype *m3, int dim)
{
	int id = (*parametro).id;
	basetype total;
	basetype aux;

	//Filas que multiplica el thread
	int cant_filas = dim/CANT_THREADS;	//Cant de filas que multiplica cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	for(int i=fila_inicial;i<=fila_final;i++)
	{
		//Recorre solo algunas filas
		for(int j=0;j<dim;j++)
		{
			//Recorre todas las columnas
			total=0;
			for(int k=0;k<dim;k++)
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

void multiplicacionXTriangularUSECUENCIAL(basetype * m1, basetype * m2, basetype *m3, int dim)
{
	basetype aux;
	basetype total;
	for(int i=0;i<dim;i++)
	{
		//Recorre solo algunas filas
		for(int j=0;j<dim;j++)
		{
			//Recorre todas las columnas
			total=0;
			for(int k=0;k<dim;k++)
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

void multiplicacionXTriangularLSECUENCIAL( basetype * m1, basetype * m2, basetype *m3, int dim)
{
	basetype total;
	basetype aux;
	for(int i=0;i<dim;i++)
	{
		// Recorre solo algunas filas
		for(int j=0;j<dim;j++)
		{
			// Recorre todas las columnas
			total=0;
			for(int k=0;k<dim;k++)
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

/*
*	m2 es triangular superior
*/
void multiplicacionXTriangularL(param* parametro, basetype * m1, basetype * m2, basetype *m3, int dim)
{
	int id = (*parametro).id;
	basetype total;
	basetype aux;
	//Filas que multiplica el thread
	int cant_filas = dim/CANT_THREADS;	//Cant de filas que multiplica cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	for(int i=fila_inicial;i<=fila_final;i++)
	{
		//Recorre solo algunas filas
		for(int j=0;j<dim;j++)
		{
			//Recorre todas las columnas
			total=0;
			for(int k=0;k<N;k++)
			{
				if(i<=j)
				{
					//p=j + j*(j+1)/2;
					//if (p > 1500000)printf("%d %d\n",j,k);
					//printf("%lu \n",p);
					aux=m2[	k + j*(j+1)/2];
				}
				else
					aux = 0;
				total 	+=	m1[i*dim + k]*aux;
			}
			m3[i*dim+j] = total;
		}
	}
}

void multiplicacion(param* parametro, basetype * m1, basetype * m2, basetype * m3, int dim)
{
	int id = (*parametro).id;
	basetype total;
	int i,j,k;

	//Filas que multiplica el thread
	int cant_filas = dim/CANT_THREADS;	// Cant de filas que multiplica cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	//printf("ID: %d \t fila_inicial=%d \t fila_final=%d \n",id,fila_inicial,fila_final);

	// Multiplica A*B=C
	for(i=fila_inicial;i<=fila_final;i++)
	{
		//Recorre solo algunas filas
		for(j=0;j<dim;j++)
		{
			//Recorre todas las columnas
			total=0;
			for(k=0;k<dim;k++)
			{
				total+=m1[i*dim+k]*m2[k*dim+j];	//total=A*B
			}
			m3[i*dim+j] = total;
		}
	}
}

void prodPromLU(param* parametro)
{
	int id = (*parametro).id;
	int i;

	// Filas que suma el thread
	int cant_filas = NT/CANT_THREADS;	// Cant de filas que suma cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	sumaL[id]=0;
	sumaU[id]=0;


	for(i=fila_inicial;i<=fila_final;i++)
	{
		//Recorre solo algunas filas
		sumaL[id]+=UT[i];
		sumaU[id]+=LT[i];

	}
	pthread_barrier_wait(&barrera); //Espera a que todas terminen
    if(id == 0)						//si es el hilo principal
    {
    	for (int i = 0; i < CANT_THREADS; ++i)
    	{
    		promL+=sumaL[i];
    		promU+=sumaU[i];
    	}
    	promL/=CANT_THREADS;
    	promU/=CANT_THREADS;
    	prodLU=promU*promL;
    }
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

void promedioB(param* parametro)
{
	int id = (*parametro).id;
	basetype total;
	int i,j;

	// Filas que multiplica el thread
	int cant_filas = N/CANT_THREADS;	// Cant de filas que multiplica cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	//printf("ID: %d \t fila_inicial=%d \t fila_final=%d \n",id,fila_inicial,fila_final);
	total=0;
	// Multiplica A*B=C

	for(i=fila_inicial;i<=fila_final;i++)
	{
		// Recorre solo algunas filas
		for(j=0;j<N;j++)
		// Recorre todas las columnas
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
    	promB /=(N*N);
    }
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

void prod_escalar (param * parametro, basetype * m1, basetype a, basetype * m2)
{
	int id = (*parametro).id;
	int cant_filas = N/CANT_THREADS;	// Cant de filas que multiplica cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	for (int i = fila_inicial; i < fila_final; ++i)
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

void suma_matriz (param * parametro, basetype * m1, basetype * m2, basetype * res)
{
	int id = (*parametro).id;
	int cant_filas = N/CANT_THREADS;	// Cant de filas que multiplica cada thread
	int fila_inicial = id*cant_filas;
	int fila_final = fila_inicial + cant_filas -1;

	for (int i = fila_inicial; i < fila_final; ++i)
	{
		for (int j = 0; j < N; ++j)
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