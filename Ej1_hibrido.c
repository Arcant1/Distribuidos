#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#include <omp.h>

//Para calcular tiempo
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

int main(int argc, char** argv){

	int miID; int cantProcesos;
	int N; // Dimension de la matriz
	int dim; // Cantidad total de datos matriz 	
	int buf; // Cantidad de elementos por proceso
	int dimTriangular;
	double *A_buf, *L_buf, *D_buf,*ab_temp,*lc_temp,*du_temp;
	double *A,*B,*C,*D,*L,*U,*M;
	double u=0.0,l=0.0;
	double tiempo_paral, tiempo_balance;
	double t2;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&miID);
	MPI_Comm_size(MPI_COMM_WORLD,&cantProcesos);

	/***** PROGRAMA *****/
	if (argc < 2){
		printf("\n Falta un parametro ");
		printf("\n 1. Dimension de la matriz ");
		return 0;
	}

	N = atoi(argv[1]);
 	dim=N*N; // Cantidad total de datos matriz
 	buf=dim/cantProcesos; // Cantidad de elementos por bloque x Cantidad de bloques por proceso 	
	dimTriangular = (N+1)*N/2;
	int posicion = (N/cantProcesos)*miID;

	int CANT_THREADS = atoi(argv[2]);
    omp_set_num_threads(CANT_THREADS);

 	A_buf=(double*)malloc(sizeof(double)*buf);
	L_buf=(double*)malloc(sizeof(double)*buf);
	D_buf=(double*)malloc(sizeof(double)*buf);

	ab_temp=(double*)malloc(sizeof(double)*buf);
	lc_temp=(double*)malloc(sizeof(double)*buf);
	du_temp=(double*)malloc(sizeof(double)*buf);
	
	B=(double*)malloc(sizeof(double)*dim);
	C=(double*)malloc(sizeof(double)*dim);
	U=(double*)malloc(sizeof(double)*dimTriangular);      

	if(miID==0) { // El proceso con ID=0 inicializa y distribuye los datos
	
		A=(double*)malloc(sizeof(double)*dim);
		L=(double*)malloc(sizeof(double)*dim);
		D=(double*)malloc(sizeof(double)*dim);
		M=(double*)malloc(sizeof(double)*dim);

	  	for(int i=0;i<N;i++){
     		for(int j=0;j<N;j++){
       			A[i*N+j]=1;
       			D[i*N+j]=1;
	   			if(i>=j){
	   				L[i*N+j]=1;
				}else{
					L[i*N+j]=0;	
				}
       			B[i+N*j]=1;
       			C[i+N*j]=1;
		   		U[i+j*(j+1)/2]=1;
     		}
   		}
	}

 	tiempo_paral = dwalltime();

	MPI_Scatter(A,buf,MPI_DOUBLE,A_buf,buf,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(B, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(L,buf,MPI_DOUBLE,L_buf,buf,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(C, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(D,buf,MPI_DOUBLE,D_buf,buf,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(U, dimTriangular, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	tiempo_balance = dwalltime();

	// Promedio de U
	double temp=0;
	t2=dwalltime();
	#pragma omp parallel for ordered reduction(+ : temp) schedule(static)
	for(int i=0;i<dimTriangular;i++)
		temp+=U[i];

	temp/=dim;	

	MPI_Allreduce(&temp, &u, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


	// Promedio de L
	temp=0;
	#pragma omp parallel for ordered reduction(+ : temp) schedule(static)
	for(int i=0;i<buf/N;i++){
		for (int j = 0; j < posicion+i+1; j++){
			temp+=L_buf[i*N+j];
		}
	}

	temp/=dim;

	MPI_Allreduce(&temp, &l, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


	#pragma omp parallel
	{
	#pragma omp parallel for
	for(int i=0;i<buf/N;i++){
		for(int j=0;j<N;j++){
			ab_temp[i*N+j]=0;
			for(int k=0;k<N;k++){
				ab_temp[i*N+j]= ab_temp[i*N+j] + A_buf[i*N+k]*B[k+j*N];
			}
		}
	}

	#pragma omp parallel for
	for (int i=0;i<buf;i++){
		ab_temp[i]*=u*l;
	}

	#pragma omp parallel for
	for(int i=0;i<buf/N;i++){
		for(int j=0;j<N;j++){
			lc_temp[i*N+j]=0;
			for(int k=0;k<posicion+i+1;k++){
				lc_temp[i*N+j]= lc_temp[i*N+j] + L_buf[i*N+k]*C[k+j*N];
			}
		}
	}

	#pragma omp parallel for
	for (int i=0;i<buf;i++){
		lc_temp[i]*=u*l;
	}

	#pragma omp parallel for
	for(int i=0;i<buf/N;i++){
		for(int j=0;j<N;j++){
			du_temp[i*N+j]=0;
			for(int k=0;k<j+1;k++){
				du_temp[i*N+j]= du_temp[i*N+j] + D_buf[i*N+k]*U[k+j*(j+1)/2];
			}
		}
	}

	#pragma omp parallel for
	for (int i=0;i<buf;i++){
		du_temp[i]*=u*l;
	}

	#pragma omp parallel for
	for (int i=0;i<buf;i++){
		ab_temp[i]+=lc_temp[i]+du_temp[i];
	}
	}
	t2=dwalltime()-t2;

	printf("Tiempo de procesamiento en cada host %f \n", t2);
	MPI_Gather(ab_temp, buf, MPI_DOUBLE, M, buf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	
	tiempo_paral = dwalltime() - tiempo_paral;
	tiempo_balance = dwalltime() - tiempo_balance;
	printf("Tiempo en segundos paralelo %f \n", tiempo_paral);
	printf("Tiempo en segundos balance %f \n", tiempo_balance);
	if(miID==0)printf("Tiempo de overhead %f \n", tiempo_paral - t2);

	//***** FIN PROGRAMA ****
	
	free(B);
	free(C);
	free(U);
	free(ab_temp);
	free(lc_temp);
	free(du_temp);
	free(A_buf);
	free(L_buf);
	free(D_buf);
	if(miID==0){
		free(A);
		free(L);
		free(D);
		free(M);		
	}
	MPI_Finalize();
	return(0);
}


