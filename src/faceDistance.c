// This function calculates the face distance of each face (i.e. edge length for 2-d)
#include "global.h"

int faceDistance()
{

    printf("*****CALCULATING FACE DISTANCE*****\n");

    // Allocates the memory to store face distance
    d_f = (double *)malloc((f + 1) * sizeof(double));

    // Looping over all the face to find the face length from the given nodes connecting the face
    for (i = 1; i <= f; i++)
    {
        d_f[i] = pow(node[face[i].node[1]].xn - node[face[i].node[2]].xn, 2);
        d_f[i] = d_f[i] + pow(node[face[i].node[1]].yn - node[face[i].node[2]].yn, 2);
        d_f[i] = pow(d_f[i], 0.5);
    }

    a_c = (double *)malloc((c + 1) * sizeof(double));
    for (i = 1; i <= c; i++)
    {
        a_c[i] = 0.0;
        for (j = 1; j <= cell[i].num_nodes; j++)
        {
            if (j != cell[i].num_nodes)
            {
                a_c[i] += (node[cell[i].node[j]].xn*node[cell[i].node[j+1]].yn) - (node[cell[i].node[j+1]].xn*node[cell[i].node[j]].yn);
            }
            if (j == cell[i].num_nodes)
            {
                a_c[i] += (node[cell[i].node[j]].xn*node[cell[i].node[1]].yn) - (node[cell[i].node[1]].xn*node[cell[i].node[j]].yn);
            }
        }
        if (a_c[i] < 0)
        {
            a_c[i] = -a_c[i];
        }
        a_c[i] = a_c[i]/2;
    }

    // printf("%lf\n",a_c[1]);
    // printf("%lf\n",a_c[10]);
    // printf("%lf\n",a_c[100]);

    printf("*****FACE DISTANCE CALCULATED*****\n\n");

    return (0);
}