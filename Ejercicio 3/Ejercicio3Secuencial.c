/**
 * C program to count total number of even and odd elements in an array
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double dwalltime();
double initial_time;
double final_time;

int main()
{
    int i, size, even;

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

    even = 0;

    initial_time = dwalltime();
    for(i=0; i<size; i++)
    {
        /* If the current element of array is even then increment even count */
        if(arr[i]%2 == 0)
        {
            even++;
        }
    }
    final_time = dwalltime() - initial_time;

    printf("Total even elements: %d\n", even);
    printf("Total time: %f\n", final_time);

    return 0;
}

double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}