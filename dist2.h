#include <math.h>

struct p
{
	float x,y;
};

typedef struct p point2d;


__device__ float distGPU(point2d p1,point2d p2)
{
	float diffx = p2.x - p1.x;
	float diffy = p2.y - p1.y;
	float dxs = pow(diffx,2);
	float dys = pow(diffy,2);
	
	float distance = sqrt(dxs + dys);

	return distance;
}

float dist(point2d p1,point2d p2)
{

	float diffx = p2.x - p1.x;
	float diffy = p2.y - p1.y;

	float dxs = pow(diffx,2);
	float dys = pow(diffy,2);

	float distance = sqrt(dxs + dys);

	return distance;


}

float dist2(point2d p1,point2d p2)
{

        float diffx = p2.x - p1.x;
        float diffy = p2.y - p1.y;

        float dxs = pow(diffx,2);
        float dys = pow(diffy,2);

        float distance = dxs + dys;

        return distance;
}
void saveData(point2d *p, float *res,int n)
{
        int i=0;

        FILE *myfile;
	
    time_t t;
    struct tm *tm;
    char date[30], *myDir;

    t = time(NULL);
    tm = localtime(&t);
    strftime(date, sizeof(date)-1, "%d-%m-%Y_%H:%M:%S", tm);
    static const char *directory = "ResultsKrigingGPU-";

    myDir = (char *)malloc(strlen(directory)+strlen(date)+1);
    strcpy(myDir, directory);
    strcat(myDir, date);
        myfile = fopen(myDir,"w");

        fprintf(myfile,"X;Y;Z\n");

        for(i=0;i<n;i++)
                fprintf(myfile,"%f;%f;%f\n",p[i].x,p[i].y,res[i]);

        fclose(myfile);


}

/*
void saveData(point2d *p, float *res,int n)
{
        int i=0;

        FILE *myfile;
        myfile = fopen("webello.txt","w");

        fprintf(myfile,"X\t\tY\t\tZ\n");

        for(i=0;i<n;i++)
                fprintf(myfile,"%f;%f;%f\n",p[i].x,p[i].y,res[i]);

        fclose(myfile);


}
*/                                                                                                                                                                            
