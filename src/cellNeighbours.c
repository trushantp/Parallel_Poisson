#include"global.h"

int cellNeighbours()
{

 printf("\n*****FINDING NEIGHBOURING CELLS AND FACES*****\n");

 ncell = (struct neigh_cell_data *) malloc((c+1)*sizeof(struct neigh_cell_data));


 for (i=1;i<=c;i++)
 {
  ncell[i].num_of_ncells = 0;
  ncell[i].num_of_faces = 0;
 }

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

 FILE *fp;

 fp = fopen("neigh_cells.out","w");

 for (i=1;i<=c;i++)
 {
  fprintf(fp,"%d\t%d\t",ncell[i].num_of_ncells,ncell[i].num_of_faces);
  for (j=1;j<=ncell[i].num_of_faces;j++)
   fprintf(fp,"%d\t%d\t",ncell[i].face[j],ncell[i].ncells[j]);
  fprintf(fp,"\n");
 }

 fclose(fp);

 printf("*****NEIGHBOURING CELLS AND FACES FOUND*****\n\n");

 return(0);

}