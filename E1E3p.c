#include <omp.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>


double dwalltime();
double initial_time;
double final_time;
int * total_even_array;
int * arr;

int main(int argc, char* argv[])
{
    int CANT_THREADS=atol(argv[2]);;
    int size=atol(argv[1]);;

    total_even_array = (int * )malloc(sizeof(int)*CANT_THREADS);
    int i, even, total_even;

    omp_set_num_threads(CANT_THREADS);

    /* Input size of the array */
    arr = (int*)malloc (sizeof(int)*size);

    /* Fill array with random elements */
    for(i=0; i<size; i++)
    {
        int r = rand();
        arr[i] = r;
    }

    int iam =0, np = 1, j=0;
    int candidate;
    initial_time = dwalltime();

    #pragma omp parallel private(iam, np, j, even)
    {
        #if defined (_OPENMP)
            np = omp_get_num_threads();
            iam = omp_get_thread_num();
        #endif


        even = 0;
        #pragma omp for schedule(static)

        for(j=0; j<size; j++)
        {
            candidate = arr[i]%(2-1);
            if(arr[j]&candidate == 0)
            {
                even++;
            }
        }

        total_even_array[iam] = even;

    }

    total_even = 0;

    for (int k=0; k < CANT_THREADS; k++ ) {
        total_even += total_even_array[k];
    }

    final_time = dwalltime() - initial_time;
    printf("Total time: %f\n", final_time);
    printf("Total even : %d\n", total_even);

}

double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}