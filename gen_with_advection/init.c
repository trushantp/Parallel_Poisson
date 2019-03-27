#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

struct node_data
{
 double xn,yn,zn;
};
struct node_data *node;

struct cell_data
{
 int num_nodes,node[10];
};
struct cell_data *cell;

struct face_data
{
 int node[3],cell[3],bc;
};
struct face_data *face;

struct neigh_cell_data
{
 int num_of_ncells,ncells[10],num_of_faces,face[10];
};
struct neigh_cell_data *ncell;

struct n2c_data
{
 int num_of_cells,cell[10],flag;
};
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

int read();
int initialization_steady();
int initialization_unsteady();
int cell_neighbours();
int face_dist();
int nodes_to_cells();
int writing_vtk();
int writing_multiple_vtk();
int initialization_lid_driven();

//---------------------------------------------------------------------------------------------------------------------------------------------------

int read()
{

 printf("*****READING THE GRID FILE STARTED*****\n");

 FILE *fp;

 fp = fopen("corrugated.grid","r");

 fscanf(fp,"%d",&n);

 node = (struct node_data *) malloc((n+1)*sizeof(struct node_data));

 for (i=1;i<=n;i++)
  fscanf(fp,"%lf %lf %lf",&node[i].xn,&node[i].yn,&node[i].zn);

 fscanf(fp,"%d",&c);

 cell = (struct cell_data *) malloc((c+1)*sizeof(struct cell_data));

 for (i=1;i<=c;i++)
 {
  fscanf(fp,"%d",&cell[i].num_nodes);
  count++;
  for (j=1;j<=cell[i].num_nodes;j++)
  {
   fscanf(fp,"%d",&cell[i].node[j]);
   count++;
  }
 }

 fscanf(fp,"%d",&f);

 face = (struct face_data *) malloc((f+1)*sizeof(struct face_data));

 for (i=1;i<=f;i++)
  fscanf(fp,"%d %d %d %d %d",&face[i].node[1],&face[i].node[2],&face[i].cell[1],&face[i].cell[2],&face[i].bc);

 fclose(fp);

 printf("*****READING THE GRID FILE COMPLETED*****\n\n");

 return(0);

}

//---------------------------------------------------------------------------------------------------------------------------------------------------

int initialization_steady()
{

 int face_flag[5];

 for (i=1;i<=4;i++)
  face_flag[i] = 0;

 for (i=1;i<=f;i++)
 {
  if ((face[i].bc == 1001 || face[i].bc == 1002 || face[i].bc == 1003 || face[i].bc == 1004 ) && face_flag[face[i].bc-1000] == 0)
  {
   printf("Dirichlet condition for wall %d:\n",face[i].bc-1000);
   printf("Enter the temperature:");
   scanf("%lf",&t[face[i].bc-1000]);
   face_flag[face[i].bc-1000] = 1;
  }

  if ((face[i].bc == 2001 || face[i].bc == 2002 || face[i].bc == 2003 || face[i].bc == 2004 ) && face_flag[face[i].bc-2000] == 0)
  {
   printf("Homogeneous Neumann Condition for wall %d:\n",face[i].bc-2000);
   printf("Heat flux = 0\nInsulated boundary condition\n");
   face_flag[face[i].bc-2000] = 1;
  }

  if ((face[i].bc == 3001 || face[i].bc == 3002 || face[i].bc == 3003 || face[i].bc == 3004 ) && face_flag[face[i].bc-3000] == 0)
  {
   printf("Non-Homogeneous Neumann Condition for wall %d:\n",face[i].bc-3000);
   printf("Enter the value of heat flux:");
   scanf("%lf",&hf[face[i].bc-3000]);
   face_flag[face[i].bc-3000] = 1;
  }

  if ((face[i].bc == 4001 || face[i].bc == 4002 || face[i].bc == 4003 || face[i].bc == 4004 ) && face_flag[face[i].bc-4000] == 0)
  {
   printf("Robin's or Mixed condition for wall %d:\n",face[i].bc-4000);
   printf("Enter value of alpha(convection coeff./thermal conductivity):");
   scanf("%lf",&alpha[face[i].bc-4000]);
   printf("Enter ambient temperature:");
   scanf("%lf",&Tinf[face[i].bc-4000]);
   face_flag[face[i].bc-4000] = 1;
  }

 }

 return(0); 

}

//---------------------------------------------------------------------------------------------------------------------------------------------------

int cell_neighbours()

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

//---------------------------------------------------------------------------------------------------------------------------------------------------

int face_dist()
{

 printf("*****CALCULATING FACE DISTANCE*****\n");

 d_f = (double *) malloc((f+1)*sizeof(double));

 for (i=1;i<=f;i++)
 {
  d_f[i] = pow(node[face[i].node[1]].xn - node[face[i].node[2]].xn,2);
  d_f[i] = d_f[i] + pow(node[face[i].node[1]].yn - node[face[i].node[2]].yn,2);
  d_f[i] = pow(d_f[i],0.5);
 }

 printf("*****FACE DISTANCE CALCULATED*****\n\n");

 return(0);

}

//---------------------------------------------------------------------------------------------------------------------------------------------------

int nodes_to_cells()
{

 printf("*****BUILDING NODE TO CELL CONNECTIVITY*****\n");

 int y;

 n2c = (struct n2c_data *) malloc ((n+1)*sizeof(struct n2c_data));

 for (i=1;i<=n;i++)
 {
  n2c[i].num_of_cells = 0;
  n2c[i].flag = 0;
 }

 for (i=1;i<=f;i++)
 {
  if (face[i].cell[1] != 0)
  {
   for (y=1;y<=2;y++)
   {
    if (face[i].node[y] != -1)
    {	
     if (n2c[face[i].node[y]].num_of_cells > 0)
     {
      n2c[face[i].node[y]].flag = 0;
      for (j=1;j<=n2c[face[i].node[y]].num_of_cells;j++)
      {
       if (n2c[face[i].node[y]].cell[j] == face[i].cell[1])
       {
        n2c[face[i].node[y]].flag = 1;
        break; 
       }
      }
      if (n2c[face[i].node[y]].flag == 0)
      {
       n2c[face[i].node[y]].num_of_cells++;
       n2c[face[i].node[y]].cell[n2c[face[i].node[y]].num_of_cells] = face[i].cell[1];
      }
     }
     else
     {
       n2c[face[i].node[y]].num_of_cells++;
       n2c[face[i].node[y]].cell[n2c[face[i].node[y]].num_of_cells] = face[i].cell[1];
     }
    }
   }
  }
  if (face[i].cell[2] != 0)
  {
   for (y=1;y<=2;y++)
   {
    if (face[i].node[y] != -1)
    {	
     if (n2c[face[i].node[y]].num_of_cells > 0)
     {
      n2c[face[i].node[y]].flag = 0;
      for (j=1;j<=n2c[face[i].node[y]].num_of_cells;j++)
      {
       if (n2c[face[i].node[y]].cell[j] == face[i].cell[2])
       {
        n2c[face[i].node[y]].flag = 1;
        break;
       }
      }
      if (n2c[face[i].node[y]].flag == 0)
      {
       n2c[face[i].node[y]].num_of_cells++;
       n2c[face[i].node[y]].cell[n2c[face[i].node[y]].num_of_cells] = face[i].cell[2];
      }
     }
     else
     {
       n2c[face[i].node[y]].num_of_cells++;
       n2c[face[i].node[y]].cell[n2c[face[i].node[y]].num_of_cells] = face[i].cell[2];  
     } 
    }
   }
  }
 }

 for(i=1;i<=n;i++)
 {
  if(n2c[i].num_of_cells == 0)
  {
   printf("Error in Number of nodes for cell %d\n",i);
   printf("Hence Exiting the program\n");
   return(1);
  }
 }

/* FILE *fp;
 
 fp = fopen("abc.out","w");

 for (i=1;i<=n;i++)
 {
  fprintf(fp,"%d\t",n2c[i].num_of_cells);
  for (j=1;j<=n2c[i].num_of_cells;j++)
   fprintf(fp,"%d\t",n2c[i].cell[j]);
  fprintf(fp,"\n");
 }*/

 printf("*****NODE TO CELL CONNECTIVITY BUILT*****\n\n");
 
 return(0);

}

//---------------------------------------------------------------------------------------------------------------------------------------------------

int writing_vtk()
{
 
 printf("*****WRITING OUTPUT FILE IN .vtk FORMAT*****\n");

 FILE *fp;
 double T_p[n];
 double dist,t_dist;

 fp = fopen("solver_general.vtk","w");

 fprintf(fp,"# vtk DataFile Version 3.0\n");
 fprintf(fp,"VTK datafile generated by Trushant\n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
 fprintf(fp,"POINTS %d double\n",n);

 for (i=1;i<=n;i++)
 {
   fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",node[i].xn,node[i].yn,node[i].zn);
 }

 fprintf(fp,"CELLS %d %d\n",c,count);

 for (i=1;i<=c;i++)
 {
  fprintf(fp,"%d\t",cell[i].num_nodes);
  for (j=1;j<=cell[i].num_nodes;j++)
  {
   fprintf(fp,"%d\t",cell[i].node[j]-1);
  }
  fprintf(fp,"\n");
 }

 fprintf(fp,"CELL_TYPES %d\n",c);

 for (i=1;i<=c;i++)
 {
  if (ncell[i].num_of_faces == 4)
   fprintf(fp,"9\n");
  if (ncell[i].num_of_faces == 3)
   fprintf(fp,"5\n");
 }

 fprintf(fp,"CELL_DATA %d\n",c);
 fprintf(fp,"SCALARS Temperature double\n");
 fprintf(fp,"LOOKUP_TABLE default\n");

 for (i=1;i<=c;i++)
 {
  fprintf(fp,"%24.16E\n",T[i]);
 }

 fprintf(fp,"POINT_DATA %d\n",n);
 fprintf(fp,"SCALARS Temperature double\n");
 fprintf(fp,"LOOKUP_TABLE default\n");

 for (i=1;i<=n;i++)
 {
  T_p[i] = 0.0;
  t_dist = 0.0;
  for (j=1;j<=n2c[i].num_of_cells;j++)
  {
   dist = pow((xc[n2c[i].cell[j]] - node[i].xn),2) + pow((yc[n2c[i].cell[j]] - node[i].yn),2);
   dist = pow(dist,0.5);
   t_dist = t_dist + dist;
   T_p[i] = T_p[i] + dist*T[n2c[i].cell[j]];
  }
  T_p[i] = T_p[i] / t_dist;
  fprintf(fp,"%24.16E\n",T_p[i]);
 }

 fclose(fp);

 printf("*****OUTPUT FILE IN .vtk FORMAT WRITTEN*****\n\n");

 return(0);

}

//---------------------------------------------------------------------------------------------------------------------------------------------------

int initialization_unsteady()
{
 
 printf("Enter the stopping time and dt:");
 scanf("%lf %lf",&t_end,&dt);

 printf("Enter diffusivity:");
 scanf("%lf",&diffusivity);

 printf("Enter relaxation factor:");
 scanf("%lf",&beta);

 printf("1 if you want to write the final file:\n");
 printf("2 if you want to write intermediate files too:\n");
 printf("Enter your choice:");
 scanf("%d",&choice);

 if (choice == 2)
 {
  printf("Number of time steps after which you want to write the output file:");
  scanf("%d",&w_timestep);
 }

 return(0);

}

//--------------------------------------------------------------------------------------------------------------------------------------------------

int writing_multiple_vtk()
{

 printf("*****WRITING OUTPUT FILE IN .vtk FORMAT FOR TIME %lf*****\n",time);

 double T_p[n+1];
 double dist,t_dist;

 char A[30];
 char B[]=".vtk";

 sprintf(A,"%lf",time);
 strcat(A,B);

 FILE *fp;

 fp = fopen(A,"w");

 fprintf(fp,"# vtk DataFile Version 3.0\n");
 fprintf(fp,"VTK datafile generated by Trushant\n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
 fprintf(fp,"POINTS %d double\n",n);

 for (i=1;i<=n;i++)
 {
   fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",node[i].xn,node[i].yn,node[i].zn);
 }

 fprintf(fp,"CELLS %d %d\n",c,count);

 for (i=1;i<=c;i++)
 {
  fprintf(fp,"%d\t",cell[i].num_nodes);
  for (j=1;j<=cell[i].num_nodes;j++)
  {
   fprintf(fp,"%d\t",cell[i].node[j]-1);
  }
  fprintf(fp,"\n");
 }

 fprintf(fp,"CELL_TYPES %d\n",c);

 for (i=1;i<=c;i++)
 {
  if (ncell[i].num_of_faces == 4)
   fprintf(fp,"9\n");
  if (ncell[i].num_of_faces == 3)
   fprintf(fp,"5\n");
 }

 fprintf(fp,"CELL_DATA %d\n",c);
 fprintf(fp,"SCALARS Temperature double\n");
 fprintf(fp,"LOOKUP_TABLE default\n");

 for (i=1;i<=c;i++)
 {
  fprintf(fp,"%24.16E\n",T[i]);
 }
 
 fprintf(fp,"POINT_DATA %d\n",n);
 fprintf(fp,"SCALARS Temperature double\n");
 fprintf(fp,"LOOKUP_TABLE default\n");

 for (i=1;i<=n;i++)
 {
  T_p[i] = 0.0;
  t_dist = 0.0;
  for (j=1;j<=n2c[i].num_of_cells;j++)
  {
   dist = pow((xc[n2c[i].cell[j]] - node[i].xn),2) + pow((yc[n2c[i].cell[j]] - node[i].yn),2);
   dist = pow(dist,0.5);
   t_dist = t_dist + dist;
   T_p[i] = T_p[i] + dist*T[n2c[i].cell[j]];
  }
  T_p[i] = T_p[i] / t_dist;
  fprintf(fp,"%24.16E\n",T_p[i]);
 }

 fclose(fp);

 printf("*****WRITING THE OUTPUT FILE COMPLETED FOR TIME %lf*****\n\n",time);

 return(0);

}

//---------------------------------------------------------------------------------------------------------------------------------------------------

int initialization_lid_driven()
{

 int face_flag[5];

 for (i=1;i<=4;i++)
  face_flag[i] = 0;

 printf("Enter the stopping time and dt:");
 scanf("%lf %lf",&t_end,&dt);

 printf("Enter relaxation factor:");
 scanf("%lf",&beta);

 for (i=1;i<=f;i++)
 {
  if (face[i].bc == 2 && face_flag[face[i].bc] == 0)
  {
   printf("Enter the Reynolds number at which the flow is happening:");
   scanf("%lf",&Re);
   face_flag[face[i].bc] = -1;
  }
 }

 return(0);

}
