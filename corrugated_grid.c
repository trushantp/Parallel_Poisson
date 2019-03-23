#include<stdio.h>
#include<math.h>
main()
{
	int i,j,k;
	int p,q,C,N,n,S,count,col,N1;
	double len,wid;
	
	count=1;
	printf("\nEnter the no. of rows:");
	scanf("%d",&q);
	printf("\nEnter the no. of columns:");
	scanf("%d",&p);
	printf("\nEnter the no. of rough mountains:");
	scanf("%d",&col);
        printf("\nLength:");
        scanf("%lf",&len);
        printf("\nWidth:");
        scanf("%lf",&wid);
	
	double xc[p+1][q+1],yc[p+1][q+1],xc1[p+1][q+1],yc1[p+1][q+1];
	double dx,dy;

	dx = len/(2*p*col);
	dy = wid/q;
	
	n=q;
	
	
	for(i=0;i<=p;i++)
	{
		for(j=0;j<=q;j++)
		{ 
			xc[i][j]=i;
			yc[i][j]=j;
			
			xc1[i][j]=i;
			yc1[i][j]=j;
			
			if((xc[i][j]-yc[i][j])>0)
			{
				xc[i][j]=0;
				yc[i][j]=0;
			}
			
			if((xc1[i][j]+yc1[i][j]<p))
			{
				xc1[i][j]=0;
				yc1[i][j]=0;
			}
			
//			printf("%lf\t%lf",xc[i][j],yc[i][j]);
			if (xc[i][j] != 0 || yc[i][j] != 0 )
			 count++;
			if (xc1[i][j] != 0 )
			 count++;
			
		}
//		printf("\n");
	}
	C= 2*col*(((q-p)*p) + ((p)*(p+1)/2));
	N=col*count-(col-1)*(q+1);
	//Basic no. of points
	N1=count-q-1;
	S=5*C - 2*p*col;
	count = 0;
	//int flag[]
	
	FILE *fp;
	FILE *pf;

	fp = fopen("grid_triangular.vtk","w");
	pf = fopen("corrugated.grid","w");

	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"Triangular grid - Uniform\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp,"POINTS %d double\n",N);
	fprintf(pf,"%d\n",N);

	double l=0;
	fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",xc[0][0],yc[0][0],l);
	fprintf(pf,"%24.16E\t%24.16E\t%24.16E\n",xc[0][0],yc[0][0],l);
	i=0;
	for (j=1;j<=q;j++)
	{
 		if( xc[i][j]==0)
		{
 		//	fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",dx*i,dy*j,l);
		fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",l,yc[i][j]*dx,l);
                fprintf(pf,"%24.16E\t%24.16E\t%24.16E\n",l,yc[i][j]*dx,l);
		}
	}
	
	for(k=0;k<=col-1;k++)
	{
		for (i=0;i<=p;i++)
		{
	 		for (j=0;j<=q;j++)
	 		{
	 			if((xc[i][j]-yc[i][j])<=0 && yc[i][j]>0 && xc[i][j]>0)
				{
	  		//	fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",dx*i,dy*j,l);
	  				fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",(xc[i][j]+k*2*p)*dx,yc[i][j]*dx,l);
	                                fprintf(pf,"%24.16E\t%24.16E\t%24.16E\n",(xc[i][j]+k*2*p)*dx,yc[i][j]*dx,l);
				}
			}
	  			
		}
		for (i=0;i<=p;i++)
		{
		 	for (j=0;j<=q;j++)
		 	{
		 		if((xc1[i][j]+yc1[i][j])>=p && xc1[i][j]>0)
				{
		  		//	fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",dx*i,dy*j,l);
		  			fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",(xc1[i][j]+p+k*2*p)*dx,yc1[i][j]*dx,l);
                                        fprintf(pf,"%24.16E\t%24.16E\t%24.16E\n",(xc1[i][j]+p+k*2*p)*dx,yc1[i][j]*dx,l);
				}
		  	}
		  		
		  		
		}
	}

	fprintf(fp,"CELLS %d %d\n",C,S);
	fprintf(pf,"%d\n",C);

	int u=0,a,b,c,d,t;
	n=q;
	a=q-p;
	t=q+1;
	b=p*t-p*(p-1)/2;
	for(k=0;k<=col-1;k++)
	{
		for (i=1;i<=p;i++)
		{
		
	 		for (j=1;j<=q;j++)
	 		{
	 			if (j==1)
	 			{
	                                //b=i-a;
	                c=(i*(i-1))/2;
	                d=(i-2)*(i-1)/2;
	 				fprintf(fp,"3\t%d\t%d\t%d\n",((i-1)*t-d)+k*N1,((i-1)*t-d+1)+k*N1,i*t-c+k*N1);
                                        fprintf(pf,"3\t%d\t%d\t%d\n",((i-1)*t-d)+k*N1+1,((i-1)*t-d+1)+k*N1+1,i*t-c+k*N1+1);
	 			}
				else if (j > 1)
				{
					c=(i*(i-1))/2;
	                d=(i-2)*(i-1)/2;
	 				fprintf(fp,"4\t%d\t%d\t%d\t%d\n",((i-1)*t-d+j-1)+k*N1,((i-1)*t-d+j)+k*N1,i*t-c+j-1+k*N1,i*t-c+j-2+k*N1);
                                        fprintf(pf,"4\t%d\t%d\t%d\t%d\n",((i-1)*t-d+j-1)+k*N1+1,((i-1)*t-d+j)+k*N1+1,i*t-c+j-1+k*N1+1,i*t-c+j-2+k*N1+1);
				}
	 		}
			q--;
		}
		
		q=n;
		
		for (i=1;i<=p;i++)
		{
		
	 		for (j=1;j<=(q-p+1);j++)
	 		{
	 			if (j==1)
	 			{
	                                //b=i-a;
	                c=(i*(i-1))/2;
	                d=(i)*(i+1)/2;
	 				fprintf(fp,"3\t%d\t%d\t%d\n",(b+(i-1)*a+c+k*N1),b+(i*a+1+d)+k*N1,b+(i*a+d)+k*N1);
                                        fprintf(pf,"3\t%d\t%d\t%d\n",(b+(i-1)*a+c+k*N1+1),b+(i*a+1+d)+k*N1+1,b+(i*a+d)+k*N1+1);
	 			}
				else if (j > 1)
				{
					c=(i*(i-1))/2;
	                d=(i)*(i+1)/2;
	                fprintf(fp,"4\t%d\t%d\t%d\t%d\n",(b+(i-1)*a+c)+j-2+k*N1,b+(i-1)*a+c+j-1+k*N1,b+(i*a+d)+j+k*N1,b+(i*a+d)+j-1+k*N1);
                        fprintf(pf,"4\t%d\t%d\t%d\t%d\n",(b+(i-1)*a+c)+j-2+k*N1+1,b+(i-1)*a+c+j-1+k*N1+1,b+(i*a+d)+j+k*N1+1,b+(i*a+d)+j-1+k*N1+1);
				}
	 		}
			q++;
		}
		
		q=n;
	}
	fprintf(fp,"CELL_TYPES %d\n",C);
	q=n;
	for(k=0;k<=col-1;k++)
	{
		for (i=1;i<=(p);i++)
		{
			for(j=1;j<=q;j++)
			{
				if(j==1)
					fprintf(fp,"5\n");
				else
					fprintf(fp,"9\n");
			}
			q--;
		}
		
		q=n;
		
		for (i=1;i<=p;i++)
		{
		
	 		for (j=1;j<=(q-p+1);j++)
	 		{
			 	if(j==1)
					fprintf(fp,"5\n");
				else
					fprintf(fp,"9\n");
			}
			q++;
		}
		q=n;
	}
fclose(fp);

fprintf(pf,"%d\n",(((4*p*(q-p) + p + q) + 2*(p*(p+1))))*col - (col-1)*q);

int z=0,bc,y;

for (i=1;i<=q;i++)
{
 count = i;
 bc = 1;
 fprintf(pf,"%d\t%d\t%d\t%d\t%d\n",count,count+1,z,count,bc);
}

int dummy=q;
int e,f;

f = q+1;

for (i=1;i<=col;i++)
{
 dummy = q;
 for (j=1;j<=p;j++)
 {
  count++;
  for (y=1;y<=dummy-1;y++)
  {
   bc = 0;
   count++;
   if (y == 1 && j == p )
   {
    f++;
   }
   bc = 0;
   fprintf(pf,"%d\t%d\t%d\t%d\t%d\n",count,count+1,count-q-(i-1)*p,f,bc);
   e = count - q - (i-1)*p;
   f++;
  }
  dummy--;
 }

// dummy = q;

 for (j=1;j<=p;j++)
 {
  count++;
  for (y=1;y<=dummy+1;y++)
  {
   bc = 0;
   count++;
   e++;
   f++;
   if (y == 1 && j == p)
   {
    f--;
   }

   if (i == col && j == p)
   {
    bc = 3;
    fprintf(pf,"%d\t%d\t%d\t%d\t%d\n",count,count+1,z,e,bc);
   }
   else
   {
    bc = 0;
    fprintf(pf,"%d\t%d\t%d\t%d\t%d\n",count,count+1,e,f,bc);
   }
  }
  f++;
  dummy++;
 }
}

 dummy=q;
 count = 0;
 e = 1;
 f = 0;


for (i=1;i<=col;i++)
{
 for (j=1;j<=p;j++)
 {
  count++;
  for (y=1;y<=dummy;y++)
  {
   count++;
   if ( y == dummy)
   {
    bc = 2;
    fprintf(pf,"%d\t%d\t%d\t%d\t%d\n",count,count+dummy,z,++f,bc);
    ++e;
   }
   else
   {
    bc = 0;
    fprintf(pf,"%d\t%d\t%d\t%d\t%d\n",count,count+dummy,++e,++f,bc);
   }
  }
  dummy--;
 }

 dummy++;

 for (j=1;j<=p;j++)
 {
  for (y=1;y<=dummy;y++)
  {
   count++;
   if ( y == dummy)
   {
    bc = 2;
    fprintf(pf,"%d\t%d\t%d\t%d\t%d\n",count,count+dummy+1,z,++f,bc);
    ++e;
   }
   else
   {
    bc = 0;
    fprintf(pf,"%d\t%d\t%d\t%d\t%d\n",count,count+dummy+1,++e,++f,bc);
   }
  }
  dummy++;
 }
 dummy--;

}

count = 1;
dummy = q;
e = 1;

for (i=1;i<=col;i++)
{
 for (j=1;j<=p;j++)
 {
  bc = 4;
  fprintf(pf,"%d\t%d\t%d\t%d\t%d\n",count,count+dummy+1,z,e,bc);
  e = e + dummy;
  count = count + dummy + 1;
  dummy--;
 }
 dummy++;
 for (j=1;j<=p;j++)
 {
  fprintf(pf,"%d\t%d\t%d\t%d\t%d\n",count,count+dummy,z,e,bc);
  e = e + dummy;
  count = count + dummy;
  dummy++;
 }
 dummy--;
}



fclose(pf);


}

