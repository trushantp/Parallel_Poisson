#include"global.h"

int initializationSteady()
{

 printf("*****INITIALIZING THE BOUNDARY CONDITIONS FOR THE SOLVER*****\n");

 int face_flag[5];
 int n;

 FILE *fp;
 fp = fopen("input.in","r");
 fscanf(fp,"%d",&n);
 int bc_number[n+1];
 double bc_prop[n+1];

 for (i=1;i<=n;i++)
    fscanf(fp,"%d %lf",&bc_number[i],&bc_prop[i]);

 for (i=1;i<=4;i++)
    face_flag[i] = 0;

 for (i=1;i<=f;i++)
 {
  if ((face[i].bc == 1001 || face[i].bc == 1002 || face[i].bc == 1003 || face[i].bc == 1004 ) && face_flag[face[i].bc-1000] == 0)
  {
    for (j=1;j<=n;j++)
    {
        if (face[i].bc == bc_number[j])
        {
            t[face[i].bc-1000] = bc_prop[j];
            break;
        }
    }
   printf("Dirichlet condition for wall %d is %lf K.\n",face[i].bc-1000,t[face[i].bc-1000]);
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
    for (j=1;j<=n;j++)
    {
        if (face[i].bc == bc_number[j])
        {
            hf[face[i].bc-1000] = bc_prop[j];
            break;
        }
    }
   printf("Non-Homogeneous Neumann Condition for wall %d is %lf.\n",face[i].bc-3000,hf[face[i].bc-3000]);
//    printf("Enter the value of heat flux:");
//    scanf("%lf",&hf[face[i].bc-3000]);
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

 fclose(fp);

 printf("*****INITIALIZED THE BOUNDARY CONDITIONS FOR THE SOLVER*****\n");

 return(0); 

}