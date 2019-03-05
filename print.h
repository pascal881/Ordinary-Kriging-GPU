void printmatrix(float *a,int m,int n,int lda)
{
	int i=0,j=0;
	for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
				printf("%f\t",a[(i*lda)+j]);
			printf("\n");
		}

	printf("\n");
}
