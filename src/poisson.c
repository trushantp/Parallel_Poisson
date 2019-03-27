#include"global.h"

struct node_data *node;
struct cell_data *cell;
struct face_data *face;
struct neigh_cell_data *ncell;
struct n2c_data *n2c;

int i,j;
int n,c,f;
int count=0;
double t[5],hf[5],alpha[5],Tinf[5];
double *d_f;
double *T,*xc,*yc,*zc;

int choice,w_timestep;
double t_end,dt,time;
double diffusivity,beta;

double Re;
double *w;

int readGrid();
int initializationSteady();
int cellNeighbours();
int faceDistance();
int solver();
int nodes2Cells();
int writingVTK();

int main()
{

 int error;

 error = readGrid();
 if (error != 0)
 {
  printf("Error in reading\nHence exiting\n");
  return(0);
 }

 error = initializationSteady();
 if (error != 0)
 {
  printf("Error in initialization\nHence exiting\n");
  return(0);
 } 

 error = cellNeighbours();
 if (error != 0)
 {
  printf("Error in assigning neighbouring cells\nHence exiting\n");
  return(0);
 }

 error = faceDistance();
 if (error != 0)
 {
  printf("Error in calculating face distance\nHence exiting\n");
  return(0);
 }

 error = solver();
 if (error != 0)
 {
  printf("Error in solver\nHence exiting\n");
  return(0);
 }

 error = nodes2Cells();
 if (error != 0)
 {
  printf("Error in assigning node to cell pair\nHence exiting\n");
  return(0);
 }

 error = writingVTK();
 if (error != 0)
 {
  printf("Error in writing vtk file\nHence exiting\n");
  return(0);
 }

 return(0);

}