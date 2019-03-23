#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"init.c"

int read();
int initialization();
int initialization_unsteady();
int cell_neighbours();
int face_dist();
int unsteady_solver();
int nodes_to_cells();
int writing_vtk();

int main()
{

 int error;

 error = read();
 if (error != 0)
 {
  printf("Error in reading\nHence exiting\n");
  return(1);
 }

 error = initialization_lid_driven();
 if (error != 0)
 {
  printf("Error in intialization - LID DRIVEN CAVITY\nHence exiting\n");
  return(1);
 }

/* error = initialization_steady();
 if (error != 0)
 {
  printf("Error in initialization\nHence exiting\n");
  return(1);
 }

 error = initialization_unsteady();
 if (error != 0)
 {
  printf("Error in initialization of unsteady term\nHence exiting\n");
  return(1);
 }*/

 error = cell_neighbours();
 if (error != 0)
 {
  printf("Error in assigning cell neighbours\nHence exiting\n");
  return(1);
 }

 error = face_dist();
 if (error != 0)
 {
  printf("Error in calculating face distances\nHence exiting\n");
  return(1);
 }

 error = nodes_to_cells();
 if (error != 0)
 {
  printf("Error in converting nodes to cell list\nHence exiting\n");
  return(1);
 }

/* error = unsteady_solver();
 if (error != 0)
 {
  printf("Error in solver for unsteady gauss seidel method\nHence exiting\n");
  return(1);
 }

 if (choice == 1)
 {
  error = writing_vtk();
  if (error != 0)
  {
   printf("Error in writing vtk file\nHence exiting\n");
   return(1);
  }
 }*/

 return(0);

}

//---------------------------------------------------------------------------------------------------------------------------------------------------

int solver_lid_driven()
{

 printf("*****SOLVER FOR LID DRIVEN CAVITY STARTED*****\n");

 double a[10],ap,ap0,S,sigma;
 double temp,old,diff,r;
 double *w_new,*w_corr,*R;
 double *vx,*vy;
 double temporary;
 int count1=0,count2=0,counter=0,error;
 double u=1;

 struct wb_data
 {
  double wbb[4];
 };
 struct wb_data *wb;

 w = (double *) malloc((c+1)*sizeof(double));
 w_new = (double *) malloc((c+1)*sizeof(double));
 w_corr = (double *) malloc((c+1)*sizeof(double));
 wb = (struct wb_data *) malloc((c+1)*sizeof(struct wb_data));
 vx = (double *) malloc((c+1)*sizeof(double));
 vy = (double *) malloc((c+1)*sizeof(double));
 R = (double *) malloc((c+1)*sizeof(double));
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
  w[i] = 0.0;
  w_new[i] = 0.0;
  R[i] = 0.0;
  w_corr[i] = 0.0;
 }

 S = 0.0;

 for (i=1;i<=f;i++)
 {
  if (face[i].bc != 0)
  {
   temp = (node[face[i].node[1]].xn + node[face[i].node[2]].xn)/2;
   temp = (node[face[i].cell[2]].xn) - temp;
   temp = pow(temp,2) + pow(((node[face[i].node[1]].yn + node[face[i].node[2]].yn)/2)-node[face[i].cell[2]].yn,2);
   if (face[i].bc == 1)
    wb[face[i].cell[2]].wbb[1] = vy[face[i].cell[2]]/temp;
   if (face[i].bc == 2)
    wb[face[i].cell[2]].wbb[2] = (vx[face[i].cell[2]] - u)/temp;
   if (face[i].bc == 3)
    wb[face[i].cell[2]].wbb[3] = -vy[face[i].cell[2]]/temp;
   if (face[i].bc == 4)
    wb[face[i].cell[2]].wbb[4] = (vy[face[i].cell[2]] - vx[face[i].cell[1]])*sqrt(2)/temp;
  }
 }

  
/* for (time=0.0;time<=t_end;time=time+dt)
 {
  counter++;
  r = 10;
  count2 = 0;
  ap0 = 1;
  while (r >= 10e-7)
  {
   r = 0;
   for (i=1;i<=c;i++)
   {
    if (count2 == 0)
     w_new[i] = w[i];

    if (ncell[i].num_of_ncells == ncell[i].num_of_faces)
    {
     ap = 0.0;
     sigma = 0.0;
     for (j=1;j<=ncell[i].num_of_faces;j++)
     {
      if (yc[i] - yc[ncell[i].ncells[j]] == 0)
      {
       temporary = xc[i] - xc[ncell[i].ncells[j]];
       if (temporary < 0)
        temporary = - temporary;
       a[j] = (((((-vx[i] - vx[ncell[i].ncells[j]])/4)*(xc[i] - xc[ncell[i].ncells[j]])))-(1/Re))*dt/(temporary*(d_f[ncell[i].face[j]]));
      }
      if (xc[i] - xc[ncell[i].ncells[j]] == 0)
      {
       temporary = yc[i] - yc[ncell[i].ncells[j]];
       if (temporary < 0)
        temporary = - temporary;
       a[j] = (((((-vy[i] - vy[ncell[i].ncells[j]])/4)*(yc[i] - yc[ncell[i].ncells[j]])))-(1/Re))*dt/(temporary*(d_f[ncell[i].face[j]])); 
      }
      temp = 0.0;
      temp = pow((xc[i] - xc[ncell[i].ncells[j]]),2);
      temp = temp + pow((yc[i] - yc[ncell[i].ncells[j]]),2);
      ap = 1 + (2*dt/(Re*d_f[ncell[i].face[j]]*temp));
      ap = ap + a[j];
      sigma = a[j]*w_new[ncell[i].ncells[j]] + sigma;
     }
     ap = ap + 1;
     sigma = sigma + ap*w_new[i];
    }

    else if (ncell[i].num_of_ncells < ncell[i].num_of_faces)
    {
     ap = 0.0;
     sigma = 0.0;
     for (j=1;j<=ncell[i].num_of_faces;j++)
     {
      if (ncell[i].ncells[j] != 0)
      {
       temp = 0.0;
       temp = pow((xc[i] - xc[ncell[i].ncells[j]]),2);
       temp = temp + pow((yc[i] - yc[ncell[i].ncells[j]]),2);
       a[j] = - diffusivity*dt/temp;
       ap = ap - a[j];
       sigma = sigma + a[j]*T_new[ncell[i].ncells[j]];
      }
      if (ncell[i].ncells[j] == 0)
      {
       if (face[ncell[i].face[j]].bc == 1001 || face[ncell[i].face[j]].bc == 1002 || face[ncell[i].face[j]].bc == 1003 || face[ncell[i].face[j]].bc == 1004)
       {
        temp = 0.0;
        temp = (node[face[ncell[i].face[j]].node[1]].xn + node[face[ncell[i].face[j]].node[2]].xn)/2;
        temp = pow(2*(temp-xc[i]),2);
        temp = temp + pow(2*(((node[face[ncell[i].face[j]].node[1]].yn + node[face[ncell[i].face[j]].node[2]].yn)/2) - yc[i]),2);
        temp = sqrt(temp);
        temp = 2*temp;
        a[j] = - diffusivity*dt/pow(temp,2);
        ap = ap - 2*a[j];
        S = 2*a[j]*t[face[ncell[i].face[j]].bc-1000];
        sigma = sigma + S;
       }
      }
     }
     ap = ap + 1;
     sigma = sigma + ap*T_new[i];
    }
    R[i] = ap0*T[i] - sigma;
    T_corr[i] = R[i]*beta/ap;
    T_new[i] = T_new[i] + T_corr[i];
    r = r + pow(R[i],2);
   }
   r = r/c;
   r = sqrt(r);
   count2++;
  }

  printf("%lf\t%d\n",time,count2);

  for (i=1;i<=c;i++)
   T[i] = T_new[i];

  if (count2 == 1 && r < 10e-7)
   break;

 }

 printf("Time : %lf\n",time);

 printf("*****SOLVER FOR LID DRIVEN CAVITY COMPLETED SUCCESSFULLY*****\n");

 return(0);
*/
}
