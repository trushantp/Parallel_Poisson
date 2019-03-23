#include<stdio.h>
#include<math.h>

int main()
{

 int i,j,n,m,count=0;

 double r,r1,theta,dr,pi,dtheta;

 printf("\nEnter the radius of outer circle :");
 scanf("%lf",&r);

 printf("Enter the radius of the inner circle :");
 scanf("%lf",&r1);

 printf("Enter the number parts to be divided into :");
 scanf("%d",&n);

 printf("Enter no. of divisions in theta: ");
 scanf("%d",&m);

 pi = atan(1.0)*4;

 dtheta = 360.0/m;

 dr = (r-r1)/n;

 double xc[n+2][m+1],yc[n+2][m+1];

 for (i=1;i<=n+1;i++)
 {
  for (j=1;j<=m;j++)
  {
    theta = pi*(j-1)*dtheta/180;
    xc[i][j] = r*cos(theta);
    yc[i][j] = r*sin(theta);
    count++;
  }
  r = r - dr;
 }
	
 FILE *fp;

 fp = fopen("grid_circular_pipe.vtk","w");

 fprintf(fp,"# vtk DataFile Version 3.0\n");
 fprintf(fp,"Circular grid - Uniform\n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
 fprintf(fp,"POINTS %d double\n",count);

 double l=0;

 for (i=1;i<=n+1;i++)
 {
  for (j=1;j<=m;j++)
  {
   fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",xc[i][j],yc[i][j],l);
  }
 }

 fprintf(fp,"CELLS %d %d\n",m*n,m*n*5);

 int u=0;

 for (i=1;i<=n;i++)
 {
  for (j=1;j<=m;j++)
  {
   if (j==m)
    fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j-1,u+j-m,u+j,u+j+m-1);
   else
    fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j-1,u+j,u+j+m,u+j+m-1);
  }
  u= m + u;
 }

 fprintf(fp,"CELL_TYPES %d\n",m*n);

 for (i=1;i<=n;i++)
 {
  for(j=1;j<=m;j++)
  {
   fprintf(fp,"9\n");
  }
 }

 fclose(fp);

 fp = fopen("circular_pipe.grid","w");

 fprintf(fp,"%d\n",count);
 
 for (i=1;i<=n+1;i++)
 {
  for (j=1;j<=m;j++)
   fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",xc[i][j],yc[i][j],l);
 }

 fprintf(fp,"%d\n",m*n);

 u=0;

 for (i=1;i<=n;i++)
 {
  for (j=1;j<=m;j++)
  {
   if (j == m)
    fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j,u+j+1-m,u+j+1,u+j+m);
   else
    fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j,u+j+1,u+j+m+1,u+j+m);
  }
  u= m + u;
 }

 fprintf(fp,"%d\n",m*(2*n+1));

 printf("1 - Dirichlet condition - Series 1001 to 1999\n"); 
 printf("2 - Homogeneous Neumann condition - Series 2001 to 2999\n");
 printf("3 - Non - Homogeneous Neumann condition - Series 3001 to 3999\n");
 printf("4 - Robin's or Mixed condition - Series 4001 to 4999\n");

 int dummy1,dummy2;
 count = 0;
 u=0;

 for (i=1;i<=n+1;i++)
 {
  if (i == 1)
  {
   printf("Enter boundary condition for outer surface(number 1):");
   scanf("%d",&dummy1);
  }
  else if (i == n+1)
  {
   printf("Enter boundary condition for inner surface(number 2):");
   scanf("%d",&dummy1);
  }
  else
   dummy1 = 0;
  for (j=1;j<=m;j++)
  {
   if (i==1 || i == n+1)
    dummy2 = 0;
   else
    dummy2 = count - m + 1;
   count++;
   if (j == m && i != n+1)
    fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",u+j,u+j+1-m,dummy2,count,dummy1);
   else
   {
    if (i == n+1)
    {
     if (j == m)
      fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",u+j,u+j-m+1,dummy2,count-m,dummy1); 
     else
      fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",u+j,u+j+1,dummy2,count-m,dummy1);
    }
    else
     fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",u+j,u+j+1,dummy2,count,dummy1);
   }
  }
  u = u + m;
 }

 count = 0;
 u=0;

 for (i=1;i<=n;i++)
 {
  for (j=1;j<=m;j++)
  {
   count++;
   if (j==1)
   {
    fprintf(fp,"%d\t%d\t%d\t%d\t0\n",count,count+m,count,count+m-1);
   }
   else
    fprintf(fp,"%d\t%d\t%d\t%d\t0\n",count,count+m,count,count-1);
  }
 } 

 fclose(fp); 

 return(0);

}

