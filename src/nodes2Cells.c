#include"global.h"

int nodes2Cells()
{

 printf("*****BUILDING NODE TO CELL CONNECTIVITY*****\n");

 int y;

 n2c = (struct n2c_data *) malloc ((n+1)*sizeof(struct n2c_data));

 for (i=1;i<=n;i++)
 {
  n2c[i].num_of_cells = 0;
  n2c[i].flag = 0;
 }

 for (i=1;i<=f;i++)
 {
  if (face[i].cell[1] != 0)
  {
   for (y=1;y<=2;y++)
   {
    if (face[i].node[y] != -1)
    {	
     if (n2c[face[i].node[y]].num_of_cells > 0)
     {
      n2c[face[i].node[y]].flag = 0;
      for (j=1;j<=n2c[face[i].node[y]].num_of_cells;j++)
      {
       if (n2c[face[i].node[y]].cell[j] == face[i].cell[1])
       {
        n2c[face[i].node[y]].flag = 1;
        break; 
       }
      }
      if (n2c[face[i].node[y]].flag == 0)
      {
       n2c[face[i].node[y]].num_of_cells++;
       n2c[face[i].node[y]].cell[n2c[face[i].node[y]].num_of_cells] = face[i].cell[1];
      }
     }
     else
     {
       n2c[face[i].node[y]].num_of_cells++;
       n2c[face[i].node[y]].cell[n2c[face[i].node[y]].num_of_cells] = face[i].cell[1];
     }
    }
   }
  }
  if (face[i].cell[2] != 0)
  {
   for (y=1;y<=2;y++)
   {
    if (face[i].node[y] != -1)
    {	
     if (n2c[face[i].node[y]].num_of_cells > 0)
     {
      n2c[face[i].node[y]].flag = 0;
      for (j=1;j<=n2c[face[i].node[y]].num_of_cells;j++)
      {
       if (n2c[face[i].node[y]].cell[j] == face[i].cell[2])
       {
        n2c[face[i].node[y]].flag = 1;
        break;
       }
      }
      if (n2c[face[i].node[y]].flag == 0)
      {
       n2c[face[i].node[y]].num_of_cells++;
       n2c[face[i].node[y]].cell[n2c[face[i].node[y]].num_of_cells] = face[i].cell[2];
      }
     }
     else
     {
       n2c[face[i].node[y]].num_of_cells++;
       n2c[face[i].node[y]].cell[n2c[face[i].node[y]].num_of_cells] = face[i].cell[2];  
     } 
    }
   }
  }
 }

 for(i=1;i<=n;i++)
 {
  if(n2c[i].num_of_cells == 0)
  {
   printf("Error in Number of nodes for cell %d\n",i);
   printf("Hence Exiting the program\n");
   return(1);
  }
 }

/* FILE *fp;
 
 fp = fopen("abc.out","w");

 for (i=1;i<=n;i++)
 {
  fprintf(fp,"%d\t",n2c[i].num_of_cells);
  for (j=1;j<=n2c[i].num_of_cells;j++)
   fprintf(fp,"%d\t",n2c[i].cell[j]);
  fprintf(fp,"\n");
 }*/

 printf("*****NODE TO CELL CONNECTIVITY BUILT*****\n\n");
 
 return(0);

}