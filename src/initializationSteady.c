// This function initializes the boundary condition for each boundary face

#include "global.h"

int initializationSteady()
{

    printf("*****INITIALIZING THE BOUNDARY CONDITIONS FOR THE SOLVER*****\n");

    // Temporary flag for faces to reduce the redundancy of checking each face for different boundary condition label
    int face_flag[5];
    // Total number of boundary conditions
    int n;

    // File pointer
    FILE *fp;
    // Opening the input file and reading it
    fp = fopen("input.in", "r");
    // Reading the total number of boundary conditions and storing it in variable n
    fscanf(fp, "%d", &n);
    // Allocating the memory for storing the bc label
    int bc_number[n + 1];
    // Allocating the memory to store the properties of the boundary condition
    double bc_prop[n + 1], bc_prop2[n + 1];

    // Reading the input.in file for boundary conditions
    for (i = 1; i <= n; i++)
    {
        // Scanning the line for boundary condition label
        fscanf(fp, "%d", &bc_number[i]);

        // Checking which boundary condition label it corresponds to, and then reads the associated boundary condition property
        if ((bc_number[i] == 1001) || (bc_number[i] == 1002) || (bc_number[i] == 1003) || (bc_number[i] == 1004))
            fscanf(fp, "%lf", &bc_prop[i]);
    }

    // Initializing face flag to 0
    for (i = 1; i <= 4; i++)
        face_flag[i] = 0;

    // Looping over all the faces
    for (i = 1; i <= f; i++)
    {
        // Checking which boundary condition label are present, intitializing the boundary condition to that face, and printing out. Flag is used so that it prints only once.
        if ((face[i].bc == 1001 || face[i].bc == 1002 || face[i].bc == 1003 || face[i].bc == 1004) && face_flag[face[i].bc - 1000] == 0)
        {
            for (j = 1; j <= n; j++)
            {
                if (face[i].bc == bc_number[j])
                {
                    t[face[i].bc - 1000] = bc_prop[j];
                    break;
                }
            }
            printf("Dirichlet condition for wall %d is %lf K.\n", face[i].bc - 1000, t[face[i].bc - 1000]);
            face_flag[face[i].bc - 1000] = 1;
        }

        if ((face[i].bc == 2001 || face[i].bc == 2002 || face[i].bc == 2003 || face[i].bc == 2004) && face_flag[face[i].bc - 2000] == 0)
        {
            printf("Homogeneous Neumann Condition for wall %d:\n", face[i].bc - 2000);
            printf("Heat flux = 0\nInsulated boundary condition\n");
            face_flag[face[i].bc - 2000] = 1;
        }

        if ((face[i].bc == 3001 || face[i].bc == 3002 || face[i].bc == 3003 || face[i].bc == 3004) && face_flag[face[i].bc - 3000] == 0)
        {
            for (j = 1; j <= n; j++)
            {
                if (face[i].bc == bc_number[j])
                {
                    hf[face[i].bc - 3000] = bc_prop[j];
                    break;
                }
            }
            printf("Non-Homogeneous Neumann Condition for wall %d is %lf.\n", face[i].bc - 3000, hf[face[i].bc - 3000]);
            face_flag[face[i].bc - 3000] = 1;
        }

        if ((face[i].bc == 4001 || face[i].bc == 4002 || face[i].bc == 4003 || face[i].bc == 4004) && face_flag[face[i].bc - 4000] == 0)
        {
            for (j = 1; j <= n; j++)
            {
                if (face[i].bc == bc_number[j])
                {
                    alpha[face[i].bc - 4000] = bc_prop[j];
                    Tinf[face[i].bc - 4000] = bc_prop2[j];
                    break;
                }
            }
            printf("Robin's or Mixed condition for wall %d:\n", face[i].bc - 4000);
            printf("Alpha(convection coeff./thermal conductivity) is %lf", alpha[face[i].bc - 4000]);
            printf("Ambient temperature is %lf", Tinf[face[i].bc - 4000]);
            face_flag[face[i].bc - 4000] = 1;
        }
    }

    // Closing the input file
    fclose(fp);

    printf("*****INITIALIZED THE BOUNDARY CONDITIONS FOR THE SOLVER*****\n");

    return (0);
}