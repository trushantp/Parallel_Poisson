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
            printf("Dirichlet condition for wall %d is %lf C.\n", face[i].bc - 1000, t[face[i].bc - 1000]);
            face_flag[face[i].bc - 1000] = 1;
        }
    }

    // Closing the input file
    fclose(fp);

    printf("*****INITIALIZED THE BOUNDARY CONDITIONS FOR THE SOLVER*****\n");

    return (0);
}