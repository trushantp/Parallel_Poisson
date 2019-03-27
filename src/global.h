#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

struct node_data
{
 double xn,yn,zn;
};
extern struct node_data *node;

struct cell_data
{
 int num_nodes,node[10];
};
extern struct cell_data *cell;

struct face_data
{
 int node[3],cell[3],bc;
};
extern struct face_data *face;

struct neigh_cell_data
{
 int num_of_ncells,ncells[10],num_of_faces,face[10];
};
extern struct neigh_cell_data *ncell;

struct n2c_data
{
 int num_of_cells,cell[10],flag;
};
extern struct n2c_data *n2c;

extern int i,j;
extern int n,c,f;
extern int count;
extern double t[5],hf[5],alpha[5],Tinf[5];
extern double *d_f;
extern double *T,*xc,*yc,*zc;

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