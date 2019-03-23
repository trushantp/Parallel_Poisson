#include<stdio.h>
#include<math.h>

int main()
{

 int i,j,n,m,count=0;

 double r,theta,dr,pi,dtheta;

 printf("\nEnter the radius of outer circle :");
 scanf("%lf",&r);

 printf("Enter the number parts to be divided into :");
 scanf("%d",&n);

 printf("Enter no. of divisions in theta: ");
 scanf("%d",&m);

 pi = atan(1.0)*4;

 dtheta = 90.0/m;

 dr = r/n;

 double xc[n+1][m+2],yc[n+1][m+2];

 for (i=1;i<=n;i++)
 {
  for (j=1;j<=m+1;j++)
  {
    theta = pi*(j-1)*dtheta/180;
    xc[i][j] = r*cos(theta);
    yc[i][j] = r*sin(theta);
    count++;
  }
  r = r - dr;
 }
	
 FILE *fp;

 fp = fopen("grid_circular.vtk","w");

 fprintf(fp,"# vtk DataFile Version 3.0\n");
 fprintf(fp,"Circular grid - Uniform\n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
 fprintf(fp,"POINTS %d double\n",count+1);

 double l=0;

 for (i=1;i<=n;i++)
 {
  for (j=1;j<=m+1;j++)
   fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",xc[i][j],yc[i][j],l);
 }

 fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",l,l,l);

 fprintf(fp,"CELLS %d %d\n",m*n,m*4+(m*(n-1)*5));

 int u=0;

 for (i=1;i<=n;i++)
 {
  for (j=1;j<=m;j++)
  {
   if (i == n)
    fprintf(fp,"3\t%d\t%d\t%d\n",u+j-1,u+j,count);
   else
    fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j-1,u+j,u+j+m+1,u+j+m);
  }
  u= m + u + 1;
 }

 fprintf(fp,"CELL_TYPES %d\n",m*n);

 for (i=1;i<=n;i++)
 {
  for(j=1;j<=m;j++)
  {
   if(i == n)
    fprintf(fp,"5\n");
   else
    fprintf(fp,"9\n");
  }
 }

 fclose(fp);

 fp = fopen("quarter.grid","w");

 fprintf(fp,"%d\n",count+1);
 
 for (i=1;i<=n;i++)
 {
  for (j=1;j<=m+1;j++)
   fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",xc[i][j],yc[i][j],l);
 }

 fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",l,l,l);

 fprintf(fp,"%d\n",m*n);

 u=0;

 for (i=1;i<=n;i++)
 {
  for (j=1;j<=m;j++)
  {
   if (i == n)
    fprintf(fp,"3\t%d\t%d\t%d\n",u+j,u+j+1,count+1);
   else
    fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j,u+j+1,u+j+m+2,u+j+m+1);
  }
  u= m + u + 1;
 }

 fprintf(fp,"%d\n",n*((2*m) + 1)); 

 printf("1 - Dirichlet condition - Series 1001 to 1999\n"); 
 printf("2 - Homogeneous Neumann condition - Series 2001 to 2999\n");
 printf("3 - Non - Homogeneous Neumann condition - Series 3001 to 3999\n");
 printf("4 - Robin's or Mixed condition - Series 4001 to 4999\n");

 int dummy1,dummy2;
 count = 0;
 u=0;

 for (i=1;i<=n;i++)
 {
  if (i == 1)
  {
   printf("Enter boundary condition for curved surface(number 3):");
   scanf("%d",&dummy1);
  }
  else
   dummy1 = 0;
  for (j=1;j<=m;j++)
  {
   if (i==1)
    dummy2 = 0;
   else 
    dummy2 = count - m + 1;
   count++;
   fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",u+j,u+j+1,dummy2,count,dummy1);
  }
  u = u + m + 1;
 }

 count = 0;
 u=0;

 for (i=1;i<=m+1;i++)
 {
  if (i == 1)
  {
   printf("Enter boundary condition for south side(number 2):");
   scanf("%d",&dummy1);
  }
  else if (i == m+1)
  {
   printf("Enter boundary condition for west side(number 1):");
   scanf("%d",&dummy1);
  }
  else
   dummy1 = 0;
  for (j=1;j<=n;j++)
  {
   count = (m*(j-1)) + i;
   if (i == m+1)
    count = count - 1;
   if (i == 1 || i == m+1)
    dummy2 = 0;
   else
    dummy2 = count - 1;
   if (j == n)
    fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",u+1,(m+1)*n+1,dummy2,count,dummy1);
   else
    fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",u+1,u+m+2,dummy2,count,dummy1);
   u = u+m+1;
  }
  u = u - n*m -(n-1);
 }

 fclose(fp); 

 return(0);

}

