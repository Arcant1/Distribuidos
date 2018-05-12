#include <omp.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>

#define CANT_THREADS 2

double dwalltime();
double initial_time;
double final_time;

int main()
{
    
    int total_even_array[CANT_THREADS];
    int i, size, even, total_even;

    omp_set_num_threads(CANT_THREADS);

    /* Input size of the array */
    printf("Enter size of the array: ");
    scanf("%d", &size);
    int arr[size];

    /* Fill array with random elements */
    for(i=0; i<size; i++)
    {
        int r = rand();
        arr[i] = r;
    }

    int iam =0, np = 1, j=0;
    initial_time = dwalltime();

    #pragma omp parallel private(iam, np, j, even)
    {
        #if defined (_OPENMP)
            np = omp_get_num_threads();
            iam = omp_get_thread_num();
        #endif

        printf("Hello from thread %d out of %d\n",iam,np);

        even = 0;
        #pragma omp for schedule(static)

        for(j=0; j<size; j++)
        {
            /* If the current element of array is even then increment even count */
            if(arr[j]%2 == 0)
            {
                even++;
            }
        }

        total_even_array[iam] = even;

    }

    total_even = 0;

    for (int k=0; k < CANT_THREADS; k++ ) {
        total_even += total_even_array[k];
        printf("Position %d and value %d\n",k,total_even_array[k]);
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