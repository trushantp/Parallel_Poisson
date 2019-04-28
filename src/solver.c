// This function solves the poisson equation with Gauss-Seidel successive over-relaxation method
#include "global.h"
#include "mpi.h"
int myid;
int size;

int solver()
{
  // printf("solver : id %d\n", myid);
  double a[10], ap, S;
  double temp, old, diff, r;
  int count1 = 0;

  // Allocating memory for temperature for each cell
  T = (double *)malloc((c + 1) * sizeof(double));
  Told = (double *)malloc((c + 1) * sizeof(double));
  // Allocating memory for centroid x, y and z co-ordinate for a cell
  xc = (double *)malloc((c + 1) * sizeof(double));
  yc = (double *)malloc((c + 1) * sizeof(double));
  zc = (double *)malloc((c + 1) * sizeof(double));

  // Looping over all the cells to find the centroid co-ordinate
  for (i = 1; i <= c; i++)
  {
    xc[i] = 0.0;
    yc[i] = 0.0;
    zc[i] = 0.0;
    for (j = 1; j <= cell[i].num_nodes; j++)
    {
      xc[i] = xc[i] + node[cell[i].node[j]].xn;
      yc[i] = yc[i] + node[cell[i].node[j]].yn;
      zc[i] = zc[i] + node[cell[i].node[j]].zn;
    }
    xc[i] = xc[i] / cell[i].num_nodes;
    yc[i] = yc[i] / cell[i].num_nodes;
    zc[i] = zc[i] / cell[i].num_nodes;
    T[i] = 0.0;
  }
  // Initializing the rms error to 10, so that the loop runs
  diff = 10;

  if(myid == 0){
    printf("*****SOLVER FOR STEADY STATE STARTED*****\n");
    // Solver loop runs till the rms error goes below 10e-6
    while (diff > 1e-6)
    {
      r = 0.0;

      // do calculation for r
      for (i = 1; i <= (int)(c/size); i++)
      {
        if (ncell[i].num_of_ncells == ncell[i].num_of_faces)
        {
          ap = 0.0;
          T[i] = 0.0;
          for (j = 1; j <= ncell[i].num_of_faces; j++)
          {
            temp = 0.0;
            temp = sqrt(pow((xc[i] - xc[ncell[i].ncells[j]]), 2) + pow((yc[i] - yc[ncell[i].ncells[j]]), 2));
            a[j] = d_f[ncell[i].face[j]] / temp;
            ap = ap + a[j];
            T[i] = T[i] + a[j] * Told[ncell[i].ncells[j]];
          }
          // Source term is -2x + 2x^2 - 2y + 2y^2
          // The analytical solution corresponding to this term is xy(1-x)(1-y) 
          S = (2 * xc[i] - 2 * pow(xc[i], 2) + 2 * yc[i] - 2 * pow(yc[i], 2)) * a_c[i];
          T[i] = (T[i] + S) / ap;
        }
        if (ncell[i].num_of_ncells < ncell[i].num_of_faces)
        {
          ap = 0.0;
          T[i] = 0.0;
          for (j = 1; j <= ncell[i].num_of_faces; j++)
          {
            if (ncell[i].ncells[j] != 0)
            {
              temp = 0.0;
              temp = sqrt(pow((xc[i] - xc[ncell[i].ncells[j]]), 2) + pow((yc[i] - yc[ncell[i].ncells[j]]), 2));
              a[j] = d_f[ncell[i].face[j]] / temp;
              ap = ap + a[j];
              T[i] = T[i] + a[j] * Told[ncell[i].ncells[j]];
            }
            if (ncell[i].ncells[j] == 0)
            {
              if (face[ncell[i].face[j]].bc == 1001 || face[ncell[i].face[j]].bc == 1002 || face[ncell[i].face[j]].bc == 1003 || face[ncell[i].face[j]].bc == 1004)
              {
                temp = 0.0;
                temp = (node[face[ncell[i].face[j]].node[1]].xn + node[face[ncell[i].face[j]].node[2]].xn) / 2;
                temp = pow((temp - xc[i]), 2);
                temp = temp + pow((((node[face[ncell[i].face[j]].node[1]].yn + node[face[ncell[i].face[j]].node[2]].yn) / 2) - yc[i]), 2);
                temp = pow(temp, 0.5);
                temp = 2 * temp;
                a[j] = d_f[ncell[i].face[j]] / temp;
                ap = ap + 2 * a[j];
                S = 2 * a[j] * t[face[ncell[i].face[j]].bc - 1000] + (2 * xc[i] - 2 * pow(xc[i], 2) + 2 * yc[i] - 2 * pow(yc[i], 2)) * a_c[i];
                T[i] = T[i] + S;
              }
            }
          }
          T[i] = T[i] / ap;
        }
        r = r + pow(T[i] - Told[i], 2);
      }

      // receive r_part from other processes and calculate diff
      for (i = 1; i < size; ++i)
      {
        double r_part = 0;
        MPI_Recv(&r_part, 1, MPI_DOUBLE, i, MY_MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        r += r_part;
      }
      diff = pow(r, 0.5);

      for (i = 1; i <= c; i++)
        Told[i] = T[i];

      // send diff to all others
      for (i = 1; i < size; ++i)
      {
        MPI_Send(&diff, 1, MPI_DOUBLE, i, MY_MPI_TAG, MPI_COMM_WORLD);
      }
      count1++;
      if (count1 % 100 == 0)
        printf("%d\t%lf\n", count1, diff);
    }

    // get T[i] values from other processes
    for (i = 1; i < size; i++){
      MPI_Recv(&T[1 + (int)(c*i/size)], (int)(c * (i + 1)/size) - (int)(c * i /size), MPI_DOUBLE, i, MY_MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    printf("*****SOLVER FOR STEADY STATE COMPLETED*****\n\n");
  }
  else{

    // Solver loop runs till the rms error goes below 10e-6
    while (diff > 1e-6)
    {
      r = 0.0;

      // do calculation for r
      for (i = (int)(c * myid /size) + 1; i <= (int)(c * (myid + 1)/size); i++)
      {
        if (ncell[i].num_of_ncells == ncell[i].num_of_faces)
        {
          ap = 0.0;
          T[i] = 0.0;
          for (j = 1; j <= ncell[i].num_of_faces; j++)
          {
            temp = 0.0;
            temp = sqrt(pow((xc[i] - xc[ncell[i].ncells[j]]), 2) + pow((yc[i] - yc[ncell[i].ncells[j]]), 2));
            a[j] = d_f[ncell[i].face[j]] / temp;
            ap = ap + a[j];
            T[i] = T[i] + a[j] * Told[ncell[i].ncells[j]];
          }
          // Source term is -2x + 2x^2 - 2y + 2y^2
          // The analytical solution corresponding to this term is xy(1-x)(1-y) 
          S = (2 * xc[i] - 2 * pow(xc[i], 2) + 2 * yc[i] - 2 * pow(yc[i], 2)) * a_c[i];
          T[i] = (T[i] + S) / ap;
        }
        if (ncell[i].num_of_ncells < ncell[i].num_of_faces)
        {
          ap = 0.0;
          T[i] = 0.0;
          for (j = 1; j <= ncell[i].num_of_faces; j++)
          {
            if (ncell[i].ncells[j] != 0)
            {
              temp = 0.0;
              temp = sqrt(pow((xc[i] - xc[ncell[i].ncells[j]]), 2) + pow((yc[i] - yc[ncell[i].ncells[j]]), 2));
              a[j] = d_f[ncell[i].face[j]] / temp;
              ap = ap + a[j];
              T[i] = T[i] + a[j] * Told[ncell[i].ncells[j]];
            }
            if (ncell[i].ncells[j] == 0)
            {
              if (face[ncell[i].face[j]].bc == 1001 || face[ncell[i].face[j]].bc == 1002 || face[ncell[i].face[j]].bc == 1003 || face[ncell[i].face[j]].bc == 1004)
              {
                temp = 0.0;
                temp = (node[face[ncell[i].face[j]].node[1]].xn + node[face[ncell[i].face[j]].node[2]].xn) / 2;
                temp = pow((temp - xc[i]), 2);
                temp = temp + pow((((node[face[ncell[i].face[j]].node[1]].yn + node[face[ncell[i].face[j]].node[2]].yn) / 2) - yc[i]), 2);
                temp = pow(temp, 0.5);
                temp = 2 * temp;
                a[j] = d_f[ncell[i].face[j]] / temp;
                ap = ap + 2 * a[j];
                S = 2 * a[j] * t[face[ncell[i].face[j]].bc - 1000] + (2 * xc[i] - 2 * pow(xc[i], 2) + 2 * yc[i] - 2 * pow(yc[i], 2)) * a_c[i];
                T[i] = T[i] + S;
              }
            }
          }
          T[i] = T[i] / ap;
        }
        r = r + pow(T[i] - Told[i], 2);
      }
      for (i = 1; i <= c; i++)
        Told[i] = T[i];
      // send r to process 0
      MPI_Send(&r, 1, MPI_DOUBLE, 0, MY_MPI_TAG, MPI_COMM_WORLD);

      // receive diff from process 0
      MPI_Recv(&diff, 1, MPI_DOUBLE, 0, MY_MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // send T[i] values to process 0
    MPI_Send(&T[(int)(c * myid /size) + 1], (int)(c * (myid + 1)/size) - (int)(c * myid /size), MPI_DOUBLE, 0, MY_MPI_TAG, MPI_COMM_WORLD);
  }

  return (0);
}