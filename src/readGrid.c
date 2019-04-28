// This function reads the grid file for a given format and stores it in memory

#include"global.h"

int readGrid()
{

 printf("*****READING THE GRID FILE STARTED*****\n");
 
 // File pointer
 FILE *fp;
 
 // Opening the file and reading it
 fp = fopen("unstructured.grid","r");
 
 // Reading the total number of nodes
 fscanf(fp,"%d",&n);
 
 // Allocating the memory for storing the co-ordinates of all the nodes
 node = (struct node_data *) malloc((n+1)*sizeof(struct node_data));
 
 // Reading the co-ordinates of all the nodes and storing it in memory
 for (i=1;i<=n;i++)
  fscanf(fp,"%lf %lf %lf",&node[i].xn,&node[i].yn,&node[i].zn);
 
 // Reading the total number of cells
 fscanf(fp,"%d",&c);
 
 // Allocating the memory for storing all the vertice location for a given cell
 cell = (struct cell_data *) malloc((c+1)*sizeof(struct cell_data));
 
 // Reading vertices information for a given cell
 for (i=1;i<=c;i++)
 {
  // Reads how many vertices constitute the cell
  fscanf(fp,"%d",&cell[i].num_nodes);
  count++;
  // Reads the vertices location
  for (j=1;j<=cell[i].num_nodes;j++)
  {
   fscanf(fp,"%d",&cell[i].node[j]);
   count++;
  }
 }
 
 // Reading the total number of faces
 fscanf(fp,"%d",&f);
 
 // Allocating the memory to store the face information
 face = (struct face_data *) malloc((f+1)*sizeof(struct face_data));
 
 // Reading the vertices which constitute the face, the cells which are sharing the face. If the face is a boundary face, one of the cell would be 0, and the other would be the cell number
 for (i=1;i<=f;i++)
  fscanf(fp,"%d %d %d %d %d",&face[i].node[1],&face[i].node[2],&face[i].cell[1],&face[i].cell[2],&face[i].bc);
 
 // Closing the grid file
 fclose(fp);

 printf("*****READING THE GRID FILE COMPLETED*****\n\n");

 return(0);

}