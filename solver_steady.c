#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"init.c"

int read();
int initialization();
int cell_neighbours();
int face_dist();
int solver();
int nodes_to_cells();
int writing_vtk();

int main()
{

 int error;

 error = read();
 if (error != 0)
 {
  printf("Error in reading\nHence exiting\n");
  return(0);
 }

 error = initialization_steady();
 if (error != 0)
 {
  printf("Error in initialization\nHence exiting\n");
  return(0);
 } 

 error = cell_neighbours();
 if (error != 0)
 {
  printf("Error in assigning neighbouring cells\nHence exiting\n");
  return(0);
 }

 error = face_dist();
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

 error = nodes_to_cells();
 if (error != 0)
 {
  printf("Error in assigning node to cell pair\nHence exiting\n");
  return(0);
 }

 error = writing_vtk();
 if (error != 0)
 {
  printf("Error in writing vtk file\nHence exiting\n");
  return(0);
 }

 return(0);

}

//---------------------------------------------------------------------------------------------------------------------------------------------------

int solver()
{

 printf("*****SOLVER FOR STEADY STATE STARTED*****\n");

 double a[10],ap,S;
 double temp,old,diff,r;
 int count1=0;

 T = (double *) malloc((c+1)*sizeof(double));
 xc = (double *) malloc((c+1)*sizeof(double));
 yc = (double *) malloc((c+1)*sizeof(double));
 zc = (double *) malloc((c+1)*sizeof(double));

 for (i=1;i<=c;i++)
 {
  xc[i] = 0.0;
  yc[i] = 0.0;
  zc[i] = 0.0;
  for (j=1;j<=cell[i].num_nodes;j++)
  {
   xc[i] = xc[i] + node[cell[i].node[j]].xn;
   yc[i] = yc[i] + node[cell[i].node[j]].yn;
   zc[i] = zc[i] + node[cell[i].node[j]].zn;
  }
  xc[i] = xc[i]/cell[i].num_nodes;
  yc[i] = yc[i]/cell[i].num_nodes;
  zc[i] = zc[i]/cell[i].num_nodes; 
  T[i] = 0.0;
 }

 diff = 10;

 while ( diff > 10e-7)
 {
  r = 0;
  for (i=1;i<=c;i++)
  {
   old = T[i];
   if (ncell[i].num_of_ncells == ncell[i].num_of_faces)
   {
    ap = 0.0;
    T[i] = 0.0;
    for (j=1;j<=ncell[i].num_of_faces;j++)
    {
     temp = 0.0;
     temp = pow((xc[i] - xc[ncell[i].ncells[j]]),2);
     temp = temp + pow((yc[i] - yc[ncell[i].ncells[j]]),2);
     temp = pow(temp,0.5);
     a[j] = d_f[ncell[i].face[j]] / temp;
     ap = ap + a[j];
     T[i] = T[i] + a[j]*T[ncell[i].ncells[j]];
    }
    S = 0.0;
    T[i] = (T[i] + S) / ap;
   }
   if (ncell[i].num_of_ncells < ncell[i].num_of_faces)
   {
    ap = 0.0;
    T[i] = 0.0;
    for (j=1;j<=ncell[i].num_of_faces;j++)
    {
     if (ncell[i].ncells[j] != 0)
     {     
      temp = 0.0;
      temp = pow((xc[i] - xc[ncell[i].ncells[j]]),2);
      temp = temp + pow((yc[i] - yc[ncell[i].ncells[j]]),2);
      temp = pow(temp,0.5);
      a[j] = d_f[ncell[i].face[j]] / temp;
      ap = ap + a[j];
      T[i] = T[i] + a[j]*T[ncell[i].ncells[j]];
     }
     if (ncell[i].ncells[j] == 0)
     {
      if (face[ncell[i].face[j]].bc == 1001 || face[ncell[i].face[j]].bc == 1002 || face[ncell[i].face[j]].bc == 1003 || face[ncell[i].face[j]].bc == 1004)
      {
       temp = 0.0;
       temp = (node[face[ncell[i].face[j]].node[1]].xn + node[face[ncell[i].face[j]].node[2]].xn)/2;
       temp = pow(2*(temp - xc[i]),2);
       temp = temp + pow(2*(((node[face[ncell[i].face[j]].node[1]].yn + node[face[ncell[i].face[j]].node[2]].yn)/2) - yc[i]),2);
       temp = pow(temp,0.5);
       temp = 2*temp;
       a[j] = d_f[ncell[i].face[j]] / temp;
       ap = ap + 2*a[j];
       S = 2*a[j]*t[face[ncell[i].face[j]].bc-1000];
       T[i] = T[i] + S;
      }
      if (face[ncell[i].face[j]].bc == 2001 || face[ncell[i].face[j]].bc == 2002 || face[ncell[i].face[j]].bc == 2003 || face[ncell[i].face[j]].bc == 2004)
       ;
      if (face[ncell[i].face[j]].bc == 3001 || face[ncell[i].face[j]].bc == 3002 || face[ncell[i].face[j]].bc == 3003 || face[ncell[i].face[j]].bc == 3004)
      {
       S = - hf[face[ncell[i].face[j]].bc-3000]*d_f[ncell[i].face[j]];
       T[i] = T[i] + S;
      }
      if (face[ncell[i].face[j]].bc == 4001 || face[ncell[i].face[j]].bc == 4002 || face[ncell[i].face[j]].bc == 4003 || face[ncell[i].face[j]].bc == 4004)
      { 
       temp = 0.0;
       temp = (node[face[ncell[i].face[j]].node[1]].xn + node[face[ncell[i].face[j]].node[2]].xn)/2;
       temp = pow(2*(temp - xc[i]),2);
       temp = temp + pow(2*(((node[face[ncell[i].face[j]].node[1]].yn + node[face[ncell[i].face[j]].node[2]].yn)/2) - yc[i]),2);
       temp = pow(temp,0.5);
       temp = 2*temp;
       a[j] = d_f[ncell[i].face[j]] / temp;
       ap = ap - a[j]*((2*alpha[face[ncell[i].face[j]].bc-4000]*temp)/(2 - alpha[face[ncell[i].face[j]].bc-4000]*temp));
       S = - a[j]*(2*temp*alpha[face[ncell[i].face[j]].bc-4000]*Tinf[face[ncell[i].face[j]].bc-4000])/(2 - temp*alpha[face[ncell[i].face[j]].bc-4000]);
       T[i] = T[i] + S;
      }
     }
    }
    T[i] = T[i]/ap;
   }
   r = r + pow(T[i]-old,2);
  }
  diff = pow(r,0.5);
  count1++;
  printf("%d\t%lf\n",count1,diff);
 }

 printf("*****SOLVER FOR STEADY STATE COMPLETED*****\n\n");

 return(0);

}

//---------------------------------------------------------------------------------------------------------------------------------------------------
