#include<stdio.h>
#include<math.h>

int main()

{

 int i,j;
 int p,q,C,N,n,S,count;
	
 count=0;
	
 printf("\nEnter the no. of rows:");
 scanf("%d",&q);
 printf("\nEnter the no. of columns:");
 scanf("%d",&p);

 double xc[p+1][q+1],yc[p+1][q+1];

 n=p;

 for(i=0;i<=q;i++)
 {
  for(j=0;j<=p;j++)
  { 
   xc[i][j]=i;
   yc[i][j]=j;
   if((xc[i][j]+yc[i][j])>p)
   {
    xc[i][j]=0;
    yc[i][j]=0;
   }
  }
 }

 C=(n)*(n+1)/2;
 N=(n+2)*(n+1)/2;
 S=5*N-n;	
	
 FILE *fp;

 fp = fopen("grid_triangular.vtk","w");

 fprintf(fp,"# vtk DataFile Version 3.0\n");
 fprintf(fp,"Triangular grid - Uniform\n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
 fprintf(fp,"POINTS %d double\n",N);

 double l=0;
 fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",xc[0][0],yc[0][0],l);
 for (i=0;i<=p;i++)
 {
  for (j=0;j<=q;j++)
  {
   if((xc[i][j]+yc[i][j])<=p && (xc[i][j]+yc[i][j])>0)
    fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",xc[i][j],yc[i][j],l);
  } 		
 }

 fprintf(fp,"CELLS %d %d\n",C,S);

 int u=0;

 for (i=1;i<=p;i++)
 {
  for (j=1;j<=q;j++)
  {
   if (j==q)
    fprintf(fp,"3\t%d\t%d\t%d\n",u+j-1,u+j,u+j+q);
   else
    fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j-1,u+j,u+j+q+1,u+j+q);
  }
  u= q + u + 1;
  q--;
 }

 fprintf(fp,"CELL_TYPES %d\n",C);
 q=p;
 for (i=0;i<=p-1;i++)
 {
  for(j=0;j<=q-1;j++)
  {
   if(j==q-1)
    fprintf(fp,"5\n");
   else
    fprintf(fp,"9\n");
  }
  q--;
 }

 fclose(fp);

 fp = fopen("triangular.grid","w");

 fprintf(fp,"%d\n",N);

 q=p;

 fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",xc[0][0],yc[0][0],l);
 for (i=0;i<=p;i++)
 {
  for (j=0;j<=q;j++)
  {
   if((xc[i][j]+yc[i][j])<=p && (xc[i][j]+yc[i][j])>0)
    fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",xc[i][j],yc[i][j],l);
  } 		
 }

 fprintf(fp,"%d\n",C);

 u=0;

 for (i=1;i<=p;i++)
 {
  for (j=1;j<=q;j++)
  {
   if (j==q)
    fprintf(fp,"3\t%d\t%d\t%d\n",u+j,u+j+1,u+j+q+1);
   else
    fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j,u+j+1,u+j+q+2,u+j+q+1);
  }
  u= q + u + 1;
  q--;
 }

 q=p;

 int dummy1 = -1;
 int dummy2;
 count=1;
 u=1;

 fprintf(fp,"%d\n",((p+1)*q + p));

 printf("1 - Dirichlet condition - Series 1001 to 1999\n"); 
 printf("2 - Homogeneous Neumann condition - Series 2001 to 2999\n");
 printf("3 - Non - Homogeneous Neumann condition - Series 3001 to 3999\n");
 printf("4 - Robin's or Mixed condition - Series 4001 to 4999\n");

 for (i=0;i<=p;i++)
 {
  if (i == 0)
  {
   printf("Enter the boundary condition for the west wall(number 1):");
   scanf("%d",&dummy2);
  }
  else
   dummy2 = 0;
  for (j=1;j<=q;j++)
  {
   fprintf(fp,"%d\t%d\t",u,u+1);
   if (i == 0)
    dummy1 = 0;
   else
    dummy1 = count-q-1;
   fprintf(fp,"%d\t%d\t%d\n",dummy1,count,dummy2);
   u++;
   count++;
  }
  u++;
  q--;
 }

 q=p;
 u=1;
 count=1;

 for (j=0;j<=q;j++)
 {
  if (j == 0)
  {
   printf("Enter the boundary condition for the east wall(number 2):");
   scanf("%d",&dummy2);
  }
  else
   dummy2 = 0;
  for (i=1;i<=p;i++)
  {
   fprintf(fp,"%d\t%d\t",u,u+q+2-i);
   if (j == 0)
    dummy1 = 0;
   else
    dummy1 = count-1;
   fprintf(fp,"%d\t%d\t%d\n",dummy1,count,dummy2);
   count = count + q + 1 - i;
   u = u+q+2-i;
  }
  p--;
  count = q - p + 1;
  u = q - p + 1;
 }

 p = q;

 printf("Enter the boundary condition for hypotenuse wall(number 3):");
 scanf("%d",&dummy2);

 u = p+1;
 count = p;
 dummy1 = 0;

 for (i=1;i<=p;i++)
 {
  fprintf(fp,"%d\t%d\t",u,u+p+1-i);
  fprintf(fp,"%d\t%d\t%d\n",dummy1,count,dummy2);
  count = count + p -i;
  u = u+p+1-i;
 }

 fclose(fp);

}

