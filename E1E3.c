

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double dwalltime();
double initial_time;
double final_time;
int * arr;
int main(int argc, char* argv[])
{
    int i, size, even,candidate;
    size = atol(argv[1]);
    arr = (int*)malloc(size*sizeof(int));

    for(i=0; i<size; i++)
    {
        int r = rand();
        arr[i] = r;
    }

    even = 0;

    initial_time = dwalltime();
    for(i=0; i<size; i++)
    {
        candidate = arr[i]%(2-1);
        if(arr[i]&candidate == 0)
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