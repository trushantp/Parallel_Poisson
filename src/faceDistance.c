#include"global.h"

int faceDistance()
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