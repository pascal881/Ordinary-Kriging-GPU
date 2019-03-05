#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cusolverDn.h>
#include <cublas_v2.h>
#include <time.h>

#include "mylib.h"
#include "print.h"
#include "dist2.h"

#define MAXBUFF 4096

__global__ void computeMatrix(float* matrice,int m,point2d *kp)
{
	int i = threadIdx.x + (blockIdx.x * blockDim.x);
	int lda = m;
	
	
	if(i>=1 && i<m)
	{
	float w,s,d;
	
	int j = 1;
	
	for(j = 1; j < m ; j++)
	{
		int tempi,tempj,mm;
		tempi = i-1;
		tempj = j-1;

		mm = lda-1;

	

		d = distGPU(kp[i],kp[j]);
		s = 10;
		w = s * exp(-3*(d/10));
		
		
		//Computing the index for access to matrix (vector' format)
		matrice[(tempi*lda)+tempj] = w;
		
		matrice[(tempi*lda)+m-1] = 1;
		matrice[(mm*lda)+tempj] = 1;
		
	}
	
	}
}

__global__ void computeMatrixD(float *d_D,point2d* d_kp,int m)
{
	int i = threadIdx.x + (blockIdx.x * blockDim.x);

	if ( i >= 1 && i < m)
		{

		float d,w,s;
		s = 10;
		d = distGPU(d_kp[i],d_kp[0]);
		w = s * exp(-3 * (d/10));
	
		d_D[i-1] = w;

		__syncthreads();
		}

	// Only master thread
	if( i == 0 )
		{
			d_D[m-1] = 1; //Lagrange factor

		}
	if(i<m)

	__syncthreads();

}

__global__ void swapValue(point2d* d_kp,point2d* d_qp,int k)
{
	int i = threadIdx.x + (blockIdx.x * blockDim.x);


	//Master thread swap the #k point in queryPoint set to first location of knowPoints	
	if( i == 0 )
		{

			d_kp[i] = d_qp[k];

		}

		__syncthreads();		
}

int main(int argc,char* argv[])
{
	char path1[MAXBUFF],path2[MAXBUFF];
	
	int flag_save = 0; // Flag to choose if save the results

	//Data for time
	cudaEvent_t start1,stop;

	point2d *kp,*qp;
	point2d *kpd,*qpd;

	float *z,*zfinal;
	float *zd,*zfinald;


	int i = 0;


	dim3 nb,nt;

	int n,m;

	float *matricegpu;
	//float *matricecpu;

	float *dd;

	/* CUSOLVER DATA */
	cusolverDnHandle_t start;
	cusolverStatus_t status;

	int *Lwork;
	float *Workspace;
	int *devIpiv; //Pivot LU
	int *devInfo;

	start = NULL;
	
	/* CUSOLVER DATA */
	int lda,ldb,nrhs;
	nrhs = 1;

	/* CUBLAS DATA */
	cublasHandle_t handlet;
        cublasStatus_t statusCub;

	if(argc>1)
		{
			flag_save = atoi(argv[1]);
			strcpy(path1,argv[2]);
			strcpy(path2,argv[3]);
			
			if(argc>4)
				nt.x = atoi(argv[4]);
		}
	else
		{
			printf("Error usage : ./exec <flagSave> <knownPointsfile> <queryPointsDataset>\n");
			exit(-1);
		}


	m = getLines(path1);
	n = getLines(path2);

//	m--;
	n--;		
	lda = m;
	ldb = m;


	if(argc<3)
	nt.x = m;

	kp = (point2d*)malloc(m*sizeof(point2d));
	qp = (point2d*)malloc(n*sizeof(point2d));
	z = (float*)malloc(m*sizeof(point2d));
	zfinal = (float*)malloc(n*sizeof(float));

	//Calcolo del numero di blocchi
	nb.x = 	m/nt.x + ( (m%nt.x) == 0 ? 0 : 1);

	printf("Flag save results : ");
	if(flag_save == 1 )
		printf("yes\n");
	else printf("no\n");

	printf("KnownPoints : %d\n",m);
	printf("QueryPoints : %d\n",n);
	printf("#Threads : %d\n",nt.x);
	printf("#Blocks : %d\n",nb.x);	
	
	/* APERTURA E LETTURA DA FILE */
	FILE *f1,*f2;
	f1 = fopen(path1,"r");
	
	for( i=0 ; i<m; i++)
		fscanf(f1,"%f;%f;%f\n",&kp[i+1].x,&kp[i+1].y,&z[i]);	
	fclose(f1);

	
	f2 = fopen(path2,"r");
	
	for(i=0;i<n;i++)
		fscanf(f2,"%f;%f;\n",&qp[i].x,&qp[i].y);
	
	fclose(f2);
/*
	printf("KnwonPoints :\nX\t\tY\t\tZ\n");
	for(i=0;i<m;i++)
		printf("%f\t%f\t%f\n",kp[i].x,kp[i].y,z[i]);
	
	printf("\nQueryPoints:\nX\t\tY\n");
	for(i=0;i<n;i++)
		printf("%f\t%f\n",qp[i].x,qp[i].y);
*/	
	cudaEventCreate(&start1);
	cudaEventCreate(&stop);

	//Allocazione dei vettori dei punti su GPU
	cudaMalloc((void**)&kpd,m*sizeof(point2d));
	cudaMalloc((void**)&qpd,n*sizeof(point2d));
	cudaMalloc((void**)&zd,m*sizeof(float));
	cudaMalloc((void**)&zfinald,n*sizeof(float));

	//Allocazione matrice 
	cudaMalloc((void**)&matricegpu,m*m*sizeof(float));
	cudaMalloc((void**)&dd,m*sizeof(float));

//	matricecpu = (float*)malloc(m*m*sizeof(float));
	//d = (float*)malloc(m*sizeof(float));

	//Copia dei vettori da CPU a GPU
	cudaMemcpy(kpd,kp,m*sizeof(point2d),cudaMemcpyHostToDevice);
	cudaMemcpy(qpd,qp,n*sizeof(point2d),cudaMemcpyHostToDevice);
	cudaMemcpy(zd,z,m*sizeof(float),cudaMemcpyHostToDevice);

	Lwork = (int*)malloc(sizeof(int));
	
	cudaEventRecord(start1);

	computeMatrix<<<nb.x,nt.x>>>(matricegpu,m,kpd);

	//cudaMemcpy(matricecpu,matricegpu,m*m*sizeof(float),cudaMemcpyDeviceToHost);
	
	 

	//Creazione contensto CUSOLVER
	if( cusolverDnCreate(&start) == CUSOLVER_STATUS_SUCCESS)
		printf("CUSOLVER OK\n");	

	status = cusolverDnSgetrf_bufferSize(start,m,m,matricegpu,lda,Lwork);
	
	//Calcolo buffersize
/*	if( status == CUSOLVER_STATUS_SUCCESS)
		printf("BUFFERSIZE OK\n");	
	printf("Workspace : %d\n",(*Lwork));
*/
	//Allocazione su GPU del workspace
	cudaMalloc((void**)&Workspace,sizeof(float)*(*Lwork));
	cudaMalloc((void**)&devIpiv,sizeof(int)*m);
	cudaMalloc((void**)&devInfo,sizeof(int));

	//Fattorizzazione
	status = cusolverDnSgetrf(start,m,m,matricegpu,lda,Workspace,devIpiv,devInfo);
	if( status == CUSOLVER_STATUS_SUCCESS)
		printf("LU Factorization ok\n");
	
	

	int k = 0;
	

	if(cublasCreate(&handlet) == CUBLAS_STATUS_SUCCESS)
		printf("CUBLAS OK\n");
	

	float o = 0;
	int incx = 1;
	int incy = 1;
	
	for(k=0;k<n;k++)
	{
	
	swapValue<<<1,1>>>(kpd,qpd,k);
	computeMatrixD<<<nb.x,nt.x>>>(dd,kpd,m);

		
	status = cusolverDnSgetrs(start,CUBLAS_OP_C,m,nrhs,matricegpu,lda,devIpiv,dd,ldb,devInfo);
	if(status != CUSOLVER_STATUS_SUCCESS)
	{
		printf("Solver cusolver problem\n");
		exit(-1);
	}


	
//	float o;
//	int incx = 1;
//	int incy = 1;

	// solve the dot product between 2 vector and save the result in "o" variable
	statusCub = cublasSdot(handlet,m-1,dd,incx,zd,incy,&o);
	if(statusCub != CUBLAS_STATUS_SUCCESS)
		{
			printf("CUBLAS problem\n");
			exit(-1);
		}

	zfinal[k] = o;
	
	}

	cudaEventRecord(stop);

	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds,start1,stop);
	printf("GPU Time elapsed : %f s\n",milliseconds*0.001);

	//Copia della matrice creata , viene copiata dalla GPU alla CPU per debug
	//cudaMemcpy(matricecpu,matricegpu,m*m*sizeof(float),cudaMemcpyDeviceToHost);
	
/*	printf("Valori finali: \nX\t\tY\t\tZ\n");
	for(i=0;i<n;i++)
		printf("%f\t%f\t%f\n",qp[i].x,qp[i].y,zfinal[i]);
	printf("\n");
*/
	
	/* Writing time elapsed */

	if(flag_save == 1)
	{	
	FILE *f6;
   	time_t t;
    	struct tm *tm;
    	 char date[30], *myDir;

    	t = time(NULL);
    	tm = localtime(&t);
    	strftime(date, sizeof(date)-1, "%d-%m-%Y_%H:%M:%S", tm);
    	static const char *directory = "ResultsKrigingTimeGPU-";

    	myDir = (char *)malloc(strlen(directory)+strlen(date)+1);
    	strcpy(myDir, directory);
    	strcat(myDir, date);
		
	f6 = fopen(myDir,"w");
	fprintf(f6,"Total time : %f s\n",milliseconds*0.001);

	fclose(f6);
	
	saveData(qp,zfinal,n);

	free(myDir);

	}

	// Free heap data
	cusolverDnDestroy(start);
	cublasDestroy(handlet);

	cudaFree(kpd);
	cudaFree(qpd);
	cudaFree(zd);	
	cudaFree(matricegpu);
	cudaFree(zfinald);
	cudaFree(Workspace);
	cudaFree(devIpiv);
	cudaFree(devInfo);
	
	free(kp);
	free(qp);
	free(Lwork);
	free(z);
	free(zfinal);	
	
	exit(0);	
}
