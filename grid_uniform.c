#include<stdio.h>
#include<stdlib.h>

int main()
{

 int i,j,p,q;

 double len,wid,dx,dy;

 printf("Length :");
 scanf("%lf",&len);

 printf("Width :");
 scanf("%lf",&wid);

 printf("Enter number of rows :");
 scanf("%d",&q);

 printf("Enter number of columns :");
 scanf("%d",&p);

 dx = len/p;
 dy = wid/q;

 FILE *fp;

 fp = fopen("uniform_grid.vtk","w");

 fprintf(fp,"# vtk DataFile Version 3.0\n");
 fprintf(fp,"Non - uniform grid\n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
 fprintf(fp,"POINTS %d double\n",(p+1)*(q+1));

 double l=0;

 for (i=0;i<=p;i++)
 {
  for (j=0;j<=q;j++)
   fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",dx*i,dy*j,l);
 }

 fprintf(fp,"CELLS %d %d\n",p*q,p*q*5);

 int u=0;

 for (i=1;i<=p;i++)
 {
  for (j=1;j<=q;j++)
  {
   fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j-1,u+j,u+j+q+1,u+j+q);
  }
  u= q + u + 1;
 }

 fprintf(fp,"CELL_TYPES %d\n",p*q);

 for (i=1;i<=p*q;i++)
  fprintf(fp,"9\n");

 fclose(fp);

 fp = fopen("uniform.grid","w");

 fprintf(fp,"%d\n",(p+1)*(q+1));
 
 for (i=0;i<=p;i++)
 {
  for (j=0;j<=q;j++)
   fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",dx*i,dy*j,l);
 }
 
 fprintf(fp,"%d\n",p*q);

 u=1;

 for (i=1;i<=p;i++)
 {
  for (j=1;j<=q;j++)
  {
   fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j-1,u+j,u+j+q+1,u+j+q);
  }
  u= q + u + 1;
 }

 fprintf(fp,"%d\n",(p*(q+1))+(q*(p+1))); 
 
 int dummy=-1;
 u=1;
 int dummy1=0,dummy2;

 printf("1 - Dirichlet condition - Series 1001 to 1999\n");
 printf("2 - Homogeneous Neumann condition - Series 2001 to 2999\n");
 printf("3 - Non - Homogeneous Neumann condition - Series 3001 to 3999\n");
 printf("4 - Mixed or Robin's Condition - Series 4001 to 4999\n");

 int d=1;
 int count = 1,flag=0;

 for (i=0;i<=p;i++)
 {
  if (i == 0 || i == p)
  {
   printf("Enter the boundary condition for this wall %d:",d);
   scanf("%d",&dummy2);
   d=d+2;
  }
  else
   dummy2 = 0;
  for (j=1;j<=q;j++)
  {
   fprintf(fp,"%d\t%d\t",u,u+1);
   if (i == 0 || i == p)
    dummy1 = 0;
   else
    dummy1 = count-p;
   if (i == p && flag == 0)
   {
    count = count - q;
    flag = 1;
   }
   fprintf(fp,"%d\t%d\t%d\n",dummy1,count,dummy2);
   u++;
   count++;
  }
  u++;
 }

 d= 2;
 u = 1;
 count =1;
 flag = 0;

 for (j=0;j<=q;j++)
 {
  if (j == 0 || j == q)
  {
   printf("Enter the boundary condition for this wall %d:",d);
   scanf("%d",&dummy2);
   d = d+2;
  }
  else
   dummy2 = 0;
  for (i=1;i<=p;i++)
  {
   fprintf(fp,"%d\t%d\t",u,u+(p+1));
   if ( j == 0 || j == q)
    dummy1=0;
   else
    dummy1 = count -1;
   if ( j == q && flag == 0)
   {
    count = count - 1;
    flag = 1;
   }
   fprintf(fp,"%d\t%d\t%d\n",dummy1,count,dummy2);
   count = count + p;
   u=u+p+1;
  }
  count = count - (p*q) +1;
  u = u - (p+1)*q + 1;
 }

 fclose(fp);

 return(0);

}
