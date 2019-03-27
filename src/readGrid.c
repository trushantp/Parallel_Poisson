#include"global.h"

int readGrid()
{

 printf("*****READING THE GRID FILE STARTED*****\n");

 FILE *fp;

 fp = fopen("uniform.grid","r");

 fscanf(fp,"%d",&n);

 node = (struct node_data *) malloc((n+1)*sizeof(struct node_data));

 for (i=1;i<=n;i++)
  fscanf(fp,"%lf %lf %lf",&node[i].xn,&node[i].yn,&node[i].zn);

 fscanf(fp,"%d",&c);

 cell = (struct cell_data *) malloc((c+1)*sizeof(struct cell_data));

 for (i=1;i<=c;i++)
 {
  fscanf(fp,"%d",&cell[i].num_nodes);
  count++;
  for (j=1;j<=cell[i].num_nodes;j++)
  {
   fscanf(fp,"%d",&cell[i].node[j]);
   count++;
  }
 }

 fscanf(fp,"%d",&f);

 face = (struct face_data *) malloc((f+1)*sizeof(struct face_data));

 for (i=1;i<=f;i++)
  fscanf(fp,"%d %d %d %d %d",&face[i].node[1],&face[i].node[2],&face[i].cell[1],&face[i].cell[2],&face[i].bc);

 fclose(fp);

 printf("*****READING THE GRID FILE COMPLETED*****\n\n");

 return(0);

}