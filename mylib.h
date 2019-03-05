/*
void print_matrix_rowmajor(char *desc, lapack_int m, lapack_int n, float* mat, lapack_int ldm)
{
         lapack_int i, j;
         printf( "\n %s\n", desc );
 
         for( i = 0; i < m; i++ ) {
                 for( j = 0; j < n; j++ ) printf( " %6.2f", mat[j+i*ldm] );
                 printf( "\n" );
        }

}
*/

int getLines(char *filename)
{
    FILE *fp = fopen(filename,"r");
    int ch = 0, cont = 1;

    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == '\n')
        {
            cont++;
        }
    }

    fclose(fp);

    return cont;
}
/*
void saveData(point2d *p, float *res,int n)
{
	int i=0;
	
	FILE *myfile;
	myfile = fopen("webello.txt","w");

	fprintf(myfile,"X\tY\tZ\n");

	for(i=0;i<n;i++)
		fprintf(myfile,"%f;%f;%f\n",p[i].x,p[i].y,res[i]);

	fclose(myfile);
	

}*/
