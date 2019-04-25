// Global file. All the variables are declared here.

#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

// Structure for storing the co-ordinates of node
struct node_data
{
 double xn,yn,zn;
};
extern struct node_data *node;

// Structure for storing the cell vertex information
struct cell_data
{
 int num_nodes,node[10];
};
extern struct cell_data *cell;

// Structure for storing the face information
struct face_data
{
 int node[3],cell[3],bc;
};
extern struct face_data *face;

// Structure for storing the neighbouring cell and face information
struct neigh_cell_data
{
 int num_of_ncells,ncells[10],num_of_faces,face[10];
};
extern struct neigh_cell_data *ncell;

// Structure for storing the cells which are shared by a vertex
struct n2c_data
{
 int num_of_cells,cell[10],flag;
};
extern struct n2c_data *n2c;

// Looping variables
extern int i,j;
// Total number of nodes, cells, and faces
extern int n,c,f;
// Temporary variable for counting something
extern int count;
// Properties of some of the boundary conditions
// t - temperature for Dirichlet
// hf - heat flux for Neumann
// alpha and Tinf - convection coeff/thermal conductivity and ambient temperature for mixed or Robins condition
extern double t[5],hf[5],alpha[5],Tinf[5];
// Face distance or face length variable
extern double *d_f;
// Area of each cell
extern double *a_c;
// Temperature, x,y, and z centroid co-ordinate
extern double *T,*Told,*xc,*yc,*zc;

extern int choice,w_timestep;
extern double t_end,dt,time;
extern double diffusivity,beta;

extern double Re;
extern double *w;

extern int readGrid();
extern int initializationSteady();
extern int initialization_unsteady();
extern int cellNeighbours();
extern int faceDistance();
extern int solver();
extern int nodes2Cells();
extern int writingVTK();

#endif // GLOBAL_H_INCLUDED