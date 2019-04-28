  
// Main file for the solver. Calls all the functions / subroutines to do the work.

#include"global.h"
#include "mpi.h"
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
double *a_c;
double *T,*Told,*xc,*yc,*zc;

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
int myid;
int size;

int main(int argc,char* argv[])
{

 MPI_Init(&argc,&argv);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
 MPI_Comm_rank(MPI_COMM_WORLD, &myid);
 // Dummy variable to check if there is any error while running the functions
 int error;
 
 // Reads the mesh and stores it in memory
 error = readGrid();
 if (error != 0)
 {
  printf("Error in reading\nHence exiting\n");
  return(0);
 }
 
 // Initializes the boundary condition for each face and internal cells too
 error = initializationSteady();
 if (error != 0)
 {
  printf("Error in initialization\nHence exiting\n");
  return(0);
 } 
 
 // Calculates which are the neighbouring cells for each cell
 error = cellNeighbours();
 if (error != 0)
 {
  printf("Error in assigning neighbouring cells\nHence exiting\n");
  return(0);
 }
 
 // Calculates the face distance (i.e. edge distance for 2-d) for each face based on the co-ordinate of the vertices
 error = faceDistance();
 if (error != 0)
 {
  printf("Error in calculating face distance\nHence exiting\n");
  return(0);
 }
 
 // Solves the Poisson equation with Gauss-seidel successive over-relaxation
 error = solver();
 if (error != 0)
 {
  printf("Error in solver\nHence exiting\n");
  return(0);
 }
 // TODO the functions below should be executed by only one process, right?
 if(myid == 0){
   // Gets the information of which cells share a node for interpolation to nodes in VTK file
   error = nodes2Cells();
   if (error != 0)
   {
    printf("Error in assigning node to cell pair\nHence exiting\n");
    return(0);
   }
   
   // Writes paraview legacy VTK file for post-processing
   error = writingVTK();
   if (error != 0)
   {
    printf("Error in writing vtk file\nHence exiting\n");
    return(0);
   }
 }
 
 MPI_Finalize();
 return(0);

}