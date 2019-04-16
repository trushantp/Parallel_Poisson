// This function calculates the neighbouring cell locations and faces constituing every cell
#include"global.h"

int cellNeighbours()
{

 printf("\n*****FINDING NEIGHBOURING CELLS AND FACES*****\n");
 
 // Allocating the memory to store the neighbouring cell information data
 ncell = (struct neigh_cell_data *) malloc((c+1)*sizeof(struct neigh_cell_data));
 
 // Looping over all the cells and initializing the number of neighbouring cells and neighbouring faces to 0
 for (i=1;i<=c;i++)
 {
  ncell[i].num_of_ncells = 0;
  ncell[i].num_of_faces = 0;
 }
 
 // Looping over all the faces and identifying the cells that share that face, and correspondingly find the neighbouring cells
 for (i=1;i<=f;i++)
 {
  if (face[i].bc != 0)
  {
   ncell[face[i].cell[2]].num_of_faces++;
   ncell[face[i].cell[2]].face[ncell[face[i].cell[2]].num_of_faces] = i;
   ncell[face[i].cell[2]].ncells[ncell[face[i].cell[2]].num_of_faces] = 0;
  }
  else
  {
   ncell[face[i].cell[1]].num_of_faces++;
   ncell[face[i].cell[1]].num_of_ncells++;
   ncell[face[i].cell[1]].face[ncell[face[i].cell[1]].num_of_faces] = i;
   ncell[face[i].cell[1]].ncells[ncell[face[i].cell[1]].num_of_faces] = face[i].cell[2];

   ncell[face[i].cell[2]].num_of_faces++;
   ncell[face[i].cell[2]].num_of_ncells++;
   ncell[face[i].cell[2]].face[ncell[face[i].cell[2]].num_of_faces] = i;
   ncell[face[i].cell[2]].ncells[ncell[face[i].cell[2]].num_of_faces] = face[i].cell[1];
  }
 }

 // File pointer
 FILE *fp;
 
 // Writing the file for neighbouring cells
 fp = fopen("neigh_cells.out","w");
 
 // Looping over the cells and writing the neighbouring cells and faces for every cell
 for (i=1;i<=c;i++)
 {
  fprintf(fp,"%d\t%d\t",ncell[i].num_of_ncells,ncell[i].num_of_faces);
  for (j=1;j<=ncell[i].num_of_faces;j++)
   fprintf(fp,"%d\t%d\t",ncell[i].face[j],ncell[i].ncells[j]);
  fprintf(fp,"\n");
 }
 
 // Closing the file for neighbouring cell
 fclose(fp);

 printf("*****NEIGHBOURING CELLS AND FACES FOUND*****\n\n");

 return(0);

}