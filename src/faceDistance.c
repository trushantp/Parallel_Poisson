// This function calculates the face distance of each face (i.e. edge length for 2-d)
#include"global.h"

int faceDistance()
{

 printf("*****CALCULATING FACE DISTANCE*****\n");
 
 // Allocates the memory to store face distance
 d_f = (double *) malloc((f+1)*sizeof(double));

 // Looping over all the face to find the face length from the given nodes connecting the face
 for (i=1;i<=f;i++)
 {
  d_f[i] = pow(node[face[i].node[1]].xn - node[face[i].node[2]].xn,2);
  d_f[i] = d_f[i] + pow(node[face[i].node[1]].yn - node[face[i].node[2]].yn,2);
  d_f[i] = pow(d_f[i],0.5);
 }

 printf("*****FACE DISTANCE CALCULATED*****\n\n");

 return(0);

}