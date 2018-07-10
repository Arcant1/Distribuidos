#include<stdio.h>
#include<stdlib.h>
#include <sys/time.h>
double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

int main(int argc, char* argv[]){
	
	unsigned long long N = atol(argv[1]);
	unsigned long m = atol(argv[2]);
	
	unsigned long* resultModulo;
	unsigned long* resultOpt;
	
	resultModulo = (unsigned long*)malloc(sizeof(unsigned long)*N);
	resultOpt = (unsigned long*)malloc(sizeof(unsigned long)*N);
	unsigned long i;
	double t1,t2;	
	

	//Calculo usando %
	double timetick = dwalltime();
	for(i=0;i<N;i++){
		resultModulo[i] = i%m;
	}
	t1=dwalltime() - timetick;
	printf("Tiempo usando %%: %f \n", t1);
	 
	//Calculo usando la equivalencia
	timetick = dwalltime();
	for(i=0;i<N;i++){
		resultOpt[i] = i&(m-1);
	}
	t2=dwalltime() - timetick;
	printf("Tiempo usando equivalencia: %f \n", t2);
	
	//Chequeo de resultados
	for(i=0;i<N;i++){
		if(resultModulo[i]!=resultOpt[i]) 
			printf("Error\n");
	}
	printf("mejora aprox= %f \n",t1/t2 );
	
	return 0;
}