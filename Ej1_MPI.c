#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

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
	int sizeTrian;
	double *A_buf, *L_buf, *D_buf,*ab_temp,*lc_temp,*du_temp;
	double *A,*B,*C,*D,*L,*U,*M;
	double u=0.0,l=0.0;
	double t_p;
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
	sizeTrian = (N+1)*N/2;
	int posicion = (N/cantProcesos)*miID;

	A_buf=(double*)malloc(sizeof(double)*buf);
	L_buf=(double*)malloc(sizeof(double)*buf);
	D_buf=(double*)malloc(sizeof(double)*buf);

	ab_temp=(double*)malloc(sizeof(double)*buf);
	lc_temp=(double*)malloc(sizeof(double)*buf);
	du_temp=(double*)malloc(sizeof(double)*buf);
	
	B=(double*)malloc(sizeof(double)*dim);
	C=(double*)malloc(sizeof(double)*dim);
	U=(double*)malloc(sizeof(double)*sizeTrian);     
	
	if(miID==0) { // El proceso con ID=0 inicializa y distribuye los datos
	
		A=(double*)malloc(sizeof(double)*dim);
		L=(double*)malloc(sizeof(double)*dim);
		D=(double*)malloc(sizeof(double)*dim);
		M=(double*)malloc(sizeof(double)*dim);
		for(int i=0;i<sizeTrian;i++){U[i]=0;}
  	for(int i=0;i<N;i++){
   		for(int j=0;j<N;j++){
   			//Ordenadas por "fila"
     			A[i*N+j]=1;
     			D[i*N+j]=1;
     			//inicializa una matriz triangular inferior
   			if(i>=j){
   				L[i*N+j]=1;
			}else{
				L[i*N+j]=0;	
			}
     			//Ordenadas por "columna"
     			B[i+N*j]=1;
     			C[i+N*j]=1;
   			//inicializa una matriz triangular superior
	   		U[i+j*(j+1)/2]=1;
   		}
 		}
	}

 	t_p = dwalltime();

	MPI_Scatter(A,buf,MPI_DOUBLE,A_buf,buf,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(B, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(L,buf,MPI_DOUBLE,L_buf,buf,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(C, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(D,buf,MPI_DOUBLE,D_buf,buf,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(U, sizeTrian, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Promedio de U

	double temp=0;

	for(int i=0;i<sizeTrian;i++)
		temp+=U[i];

	temp/=dim;	

	MPI_Allreduce(&temp, &u, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	printf("Proceso %d, Promedio u = %lf \n",miID,u);

	// Promedio de L
	temp=0;

	for(int i=0;i<buf/N;i++){
		for (int j = 0; j < i+posicion+1; j++){
			temp+=L_buf[i*N+j];
		}
	}

	temp/=dim;

	MPI_Allreduce(&temp, &l, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	printf("Proceso %d, Promedio l = %lf \n",miID,l);

//--------------------------------------------------------
	t2 = dwalltime();
	for(int i=0;i<buf/N;i++){
		for(int j=0;j<N;j++){
			ab_temp[i*N+j]=0;
			for(int k=0;k<N;k++){
				ab_temp[i*N+j]= ab_temp[i*N+j] + A_buf[i*N+k]*B[k+j*N];
			}
		}
	}

	for (int i=0;i<buf;i++){
		ab_temp[i]*=u*l;
	}

	for(int i=0;i<buf/N;i++){
		for(int j=0;j<N;j++){
			lc_temp[i*N+j]=0;
			for(int k=0;k<posicion+i+1;k++){
				lc_temp[i*N+j]= lc_temp[i*N+j] + L_buf[i*N+k]*C[k+j*N];
			}
		}
	}

	for (int i=0;i<buf;i++){
		lc_temp[i]*=u*l;
	}

	for(int i=0;i<buf/N;i++){
		for(int j=0;j<N;j++){
			du_temp[i*N+j]=0;
			for(int k=0;k<j+1;k++){
				du_temp[i*N+j]= du_temp[i*N+j] + D_buf[i*N+k]*U[k+j*(j+1)/2];
			}
		}
	}

	for (int i=0;i<buf;i++){
		du_temp[i]*=u*l;
	}

	for (int i=0;i<buf;i++){
		ab_temp[i]+=lc_temp[i]+du_temp[i];
	}
	t2=dwalltime() - t2;
	printf("Tiempo de procesamiento en cada host %f \n", t2);

	MPI_Gather(ab_temp, buf, MPI_DOUBLE, M, buf, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	t_p = dwalltime() - t_p;
	printf("Tiempo de overhead de comunicacion %f \n", t_p - t2);
	printf("Tiempo en segundos %f \n", t_p);

	
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

