#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int i,j,p,q,count;

double len,wid,dx,dy,u=1,Re,rms,vel;

double ap,an,as,aw,ae,ap0,S,beta;

double dt,time=0.0,t_end;
FILE *fp;
int main()
{
 printf("\nEnter the number of columns :");
 scanf("%d",&p);

 printf("Enter the number of rows :");
 scanf("%d",&q);

 printf("\nLength:");
 scanf("%lf",&len);

 printf("Width:");
 scanf("%lf",&wid);

 printf("Enter the velocity at which lid is driven : ");
 scanf("%lf",&vel);
 
 printf("Enter the value of reynolds number of the fluid inside the cavity : ");
 scanf("%lf",&Re);

 printf("Enter the relaxation factor : ");
 scanf("%lf",&beta);
 
 printf("Enter the end time and dt : ");
 scanf("%lf %lf",&t_end,&dt);
 
 dx = len/p;
 dy = wid/q;

 double w[p+1][q+1],w_n[p+1][q+1],xc[p+1][q+1],yc[p+1][q+1];
 
 double psi[p+1][q+1],psi_n[p+1][q+1];
 
 double wbw[q+1],wbn[p+1],wbe[q+1],wbs[p+1];
 
 double vx[p+2][q+2], vy[p+2][q+2];
 
 double dt_crit[p+1][q+1];

 double w_corr[p+1][q+1],R[p+1][q+1],sigma;
 
 for (i=1;i<=p;i++)
 {
  for (j=1;j<=q;j++)
  {
   w[i][j] = 0.0;
   w_n[i][j] = 0.0;
   psi[i][j] = 0.0;
   psi_n[p+1][q+1] = 0.0;
   xc[i][j] = i*dx - dx/2;
   yc[i][j] = j*dy - dy/2;
   S = 0.0;
  }
 }
  
 for (i=0;i<=p+1;i++)
 {
  for(j=0;j<=q+1;j++)
  {
   vx[i][j] = 0.0;
   vy[i][j] = 0.0;
  }
 }

 time = time*vel;
 dt = dt*vel;
 t_end = t_end*vel;

 for (time=0.0;time<=t_end;time=time+dt)
 {
  rms = 10e5;
  count = 0;

  while (rms > 10e-5)
  {

   rms = 0.0;

   for(i=1;i<=p;i++)
   {
    wbw[i] = (2*vy[1][i])/dx;
    wbn[i] = (2*(vx[i][q] - u))/dy;
    wbe[i] = -(2*vy[p][i])/dx;
    wbs[i] = -(2*vx[i][1])/dy;
   }
  	
   for (i=1;i<=p;i++)
   {
    for (j=1;j<=q;j++) 
    {
     ae = ((((vx[i][j] + vx[i+1][j])/4)*dy) - (1/Re))*dt/(dx*dy);
     aw = ((((-vx[i][j] - vx[i-1][j])/4)*dy) - (1/Re))*dt/(dx*dy);
     an = ((((vy[i][j] + vy[i][j+1])/4)*dx) - (1/Re))*dt/(dx*dy); 
     as = ((((-vy[i][j] - vy[i][j-1])/4)*dx) - (1/Re))*dt/(dx*dy);
     ap = ae + an + as + aw + 1 + (8*dt/(Re*dx*dy));
     ap0 = 1;
  
     if (count == 0)
      w[i][j] = w_n[i][j];
			
     if(i==1 && j==1)
     {
      ap=ap - aw -as;
      S = 2*aw*wbw[i] + 2*as*wbs[i];
      sigma = ae*w[i+1][j] + an*w[i][j+1] + S + ap*w[i][j];
     }
     else if(i==1 && j==q)
     {
      ap=ap - aw -an;
      S = 2*aw*wbw[q] + 2*an*wbn[i];
      sigma =ae*w[i+1][j] + as*w[i][j-1] + S + ap*w[i][j];
     }
     else if(i==p && j==1)
     {
      ap=ap - ae -as;
      S = 2*ae*wbe[1] + 2*as*wbs[p];
      sigma = aw*w[i-1][j] + an*w[i][j+1] + S + ap*w[i][j];
     }
     else if(i==p && j==q)
     {
      ap=ap - ae -an;
      S = 2*ae*wbe[q] + 2*an*wbn[p];
      sigma = aw*w[i-1][j] + as*w[i][j-1] + S + ap*w[i][j];
     }
     else if(i==1 && j!=1 && j!=q)
     {
      ap=ap-aw;
      S=2*aw*wbw[j];
      sigma = ae*w[i+1][j] + an*w[i][j+1] + as*w[i][j-1] + S + ap*w[i][j];
     }
     else if(i==p && j!=1 && j!=q)
     {
      ap=ap-ae;
      S=2*ae*wbe[j];
      sigma = aw*w[i-1][j] + an*w[i][j+1] + as*w[i][j-1] + S + ap*w[i][j];
     }
     else if(j==1 && i!=1 && i!=p)
     {
      ap=ap-as;
      S=2*as*wbs[i];
      sigma = aw*w[i-1][j] + ae*w[i+1][j] + an*w[i][j+1] + S + ap*w[i][j];
     }
     else if(j==q && i!=1 && i!=p)
     {
      ap=ap-an;
      S=2*an*wbn[i];
      sigma = aw*w[i-1][j] + ae*w[i+1][j] + as*w[i][j-1] + S + ap*w[i][j];
     }
     else
     {
      S = 0.0;
      sigma = aw*w[i-1][j] + ae*w[i+1][j] + an*w[i][j+1] + as*w[i][j-1] + S + ap*w[i][j];
     }

     R[i][j] = ap0*w_n[i][j] - sigma;

     w_corr[i][j] = R[i][j]*beta/ap;
     w[i][j] = w[i][j] + w_corr[i][j];

     rms=rms + R[i][j]*R[i][j];
    }
   }

   rms=rms/(p*q);
   rms=sqrt(rms);
   count++;
   printf("%d\t%lf\t%lf\n",count,time/vel,rms);
  }
  for (i=1;i<=p;i++)
  {
   for (j=1;j<=q;j++)
    w_n[i][j] = w[i][j];
  }

  count = 0;
  rms = 10e5;
  while (rms >= 10e-5)
  {  
   rms = 0.0;
   for(i=1;i<=p;i++)
   {
    for(j=1;j<=q;j++)
    {
     if(i==1 && j==1)
      psi[i][j] = (psi[i+1][j] + psi[i][j+1] + w_n[i][j]*dx*dy)/6;
   			
     else if(i==1 && j==q)
      psi[i][j] = ( psi[i+1][j] + psi[i][j-1] + w_n[i][j]*dx*dx)/6;
  			
     else if(i==p && j==1)
      psi[i][j] = (psi[i-1][j] + psi[i][j+1] + w_n[i][j]*dx*dx)/6;
  			
     else if(i==p && j==q)
      psi[i][j] = (psi[i-1][j] + psi[i][j-1] + w_n[i][j]*dx*dx)/6;
  			
     else if(i==1 && j!=1 && j!=q)
      psi[i][j] = (psi[i+1][j] + psi[i][j+1] + psi[i][j-1] + w_n[i][j]*dx*dx)/5;
   			
     else if(i==p && j!=1 && j!=q)
      psi[i][j] = (psi[i-1][j] + psi[i][j+1] + psi[i][j-1] + w_n[i][j]*dx*dx)/5;
  			
     else if(i!=1 && i!=p && j==1)
      psi[i][j] = (psi[i+1][j] + psi[i][j+1] + psi[i-1][j] + w_n[i][j]*dx*dx)/5;
  			
     else if(i!=1 && i!=p && j==q)
      psi[i][j] = (psi[i+1][j] + psi[i-1][j] + psi[i][j-1] + w_n[i][j]*dx*dx)/5;
  			
     else
      psi[i][j] = (psi[i-1][j] + psi[i+1][j] + psi[i][j+1] + psi[i][j-1] + w_n[i][j]*dx*dx)/4;
  	
     rms = rms + (psi_n[i][j] - psi[i][j])*(psi_n[i][j] - psi[i][j]);
    }
   }
   rms = rms /(p*q);
   rms = sqrt(rms);
   count++;
   printf("%d %lf\n",count,rms);
   for (i=1;i<=p;i++)
   {
    for (j=1;j<=q;j++)
     psi_n[i][j] = psi[i][j];
   } 
  }
  	
  	for(i=1;i<=p;i++)
  	{
  		for(j=1;j<=q;j++)
  		{
  			if(i==1 && j==1)
  			{
  				vx[i][j] = (psi[i][j+1] + psi[i][j])/(2*dx);
  				vy[i][j] = (-psi[i][j] - psi[i+1][j])/(2*dy);
  				vy[0][1] = -vy[i][j];
  				vy[1][0] = -vy[i][j];
  				vx[0][1] = -vx[i][j];
  				vx[1][0] = -vx[i][j];
  			}
  			
  			else if(i==1 && j==q)
  			{
  				vx[i][j] = (-psi[i][j] - psi[i][j-1])/(2*dx);
  				vy[i][j] = (-psi[i][j] - psi[i+1][j])/(2*dy);
  				vy[0][q] = -vy[i][j];
  				vy[1][q+1] = -vy[i][j];
  				vx[0][q] = -vx[i][j];
  				vx[1][q+1] = -vx[i][j];
  			}
  			
  			else if(i==p && j==1)
  			{
  				vx[i][j] = (psi[i][j+1] + psi[i][j])/(2*dx);
  				vy[i][j] = (psi[i][j] + psi[i-1][j])/(2*dy);
  				vy[p+1][1] = -vy[i][j];
  				vy[p][0] = -vy[i][j];
  				vx[p+1][1] = -vx[i][j];
  				vx[p][0] = -vx[i][j];
  			}
  			
  			else if(i==p && j==q)
  			{
  				vx[i][j] = (-psi[i][j] - psi[i][j-1])/(2*dx);
  				vy[i][j] = (psi[i][j] + psi[i-1][j])/(2*dy);
  				vy[p][q+1] = -vy[i][j];
  				vy[p+1][q] = -vy[i][j];
  				vx[p][q+1] = -vx[i][j];
  				vx[p+1][q] = -vx[i][j];
  			}
  			
  			else if(i==1 && j!=1 && j!=q)
  			{
  				vx[i][j] = (psi[i][j+1] - psi[i][j-1])/(2*dx);
  				vy[i][j] = (-psi[i][j] - psi[i+1][j])/(2*dy);
  				vy[0][j] = -vy[i][j];
  				vx[0][j] = -vx[i][j];
  			}
  			
  			else if(i==p && j!=1 && j!=q)
  			{
  				vx[i][j] = (psi[i][j+1] - psi[i][j-1])/(2*dx);
  				vy[i][j] = (psi[i][j] + psi[i-1][j])/(2*dy);
  				vy[p+1][j] = -vy[i][j];
  				vx[p+1][j] = -vx[i][j];
  			}
  			
  			else if(i!=1 && i!=p && j==1)
  			{
  				vx[i][j] = (psi[i][j+1] + psi[i][j])/(2*dx);
  				vy[i][j] = (psi[i-1][j] - psi[i+1][j])/(2*dy);
  				vy[i][0] = -vy[i][j];
  				vx[i][0] = -vx[i][j];
  			}
  			
  			else if(i!=1 && i!=p && j==q)
  			{
  				vx[i][j] = (-psi[i][j] - psi[i][j-1])/(2*dx);
  				vy[i][j] = (psi[i-1][j] - psi[i+1][j])/(2*dy);
  				vy[i][q+1] = -vy[i][j];
  				vx[i][q+1] = -vx[i][j];
  			}
  			
  			else
  			{
  				vx[i][j] = (psi[i][j+1] - psi[i][j-1])/(2*dx);
  				vy[i][j] = (psi[i-1][j] - psi[i+1][j])/(2*dy);
  			}
  			
  		}
  	}

  	printf("%lf\t%lf\n",dt,time);

}
	fp = fopen("lid_driven_cavity.vtk","w");

 fprintf(fp,"# vtk DataFile Version 3.0\n");
 fprintf(fp,"VTK datafile generated by Trushant\n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
 fprintf(fp,"POINTS %d double\n",(p+1)*(q+1));

 double l=0;

 for (i=0;i<=p;i++)
 {
  for (j=0;j<=q;j++)
   fprintf(fp,"%24.16E\t%24.16E\t%24.16E\n",dx*i,dy*j,l);
 }

 fprintf(fp,"CELLS %d %d\n",p*q,p*q*5);

 int u=0;

 for (i=1;i<=p;i++)
 {
  for (j=1;j<=q;j++)
  {
   fprintf(fp,"4\t%d\t%d\t%d\t%d\n",u+j-1,u+j,u+j+q+1,u+j+q);
  }
  u= q + u + 1;
 }

 fprintf(fp,"CELL_TYPES %d\n",p*q);

 for (i=1;i<=p*q;i++)
  fprintf(fp,"9\n");

/* fprintf(fp,"CELL_DATA %d\n",p*q);
 fprintf(fp,"SCALARS Xvelocity double\n");
 fprintf(fp,"LOOKUP_TABLE default\n");

 for (i=1;i<=p;i++)
 {
  for (j=1;j<=q;j++)
   fprintf(fp,"%24.16E\n",vx[i][j]);
 }
 
fprintf(fp,"POINT_DATA %d\n",(p+1)*(q+1));
fprintf(fp,"SCALARS Xvelocity double\n");
fprintf(fp,"LOOKUP_TABLE default\n");*/

double T_p[p+1][q+1];

/*for (i=0;i<=p;i++)
{
 for (j=0;j<=q;j++)
 {
  if (i==0 && j==0)
   T_p[i][j] = vx[i+1][j+1];
  else if (i == 0 && j !=0 && j!=q)
   T_p[i][j] = (vx[i+1][j]+vx[i+1][j+1])/2;
  else if (i == 0 && j == q)
   T_p[i][j] = vx[i+1][j];
  else if (j ==0 && i!=0 && i!=p)
   T_p[i][j] = (vx[i][j+1]+vx[i+1][j+1])/2;
  else if (j ==0 && i == p)
   T_p[i][j] = vx[i][j+1];
  else if (i ==p && j!=0 && j!=q)
   T_p[i][j] = (vx[i][j]+vx[i][j+1])/2;
  else if (i==p && j==q)
   T_p[i][j] = vx[i][j];
  else if (j==q && i!=0 && i!=p)
   T_p[i][j] = (vx[i][j]+vx[i+1][j])/2;
  else
   T_p[i][j] = (vx[i][j] + vx[i+1][j+1] + vx[i][j+1] + vx[i+1][j])/4;
  fprintf(fp,"%24.16E\n",T_p[i][j]);
 }
}*/

 fprintf(fp,"CELL_DATA %d\n",p*q);
 fprintf(fp,"SCALARS Vorticity double\n");
 fprintf(fp,"LOOKUP_TABLE default\n");

 for (i=1;i<=p;i++)
 {
  for (j=1;j<=q;j++)
   fprintf(fp,"%24.16E\n",w[i][j]);
 }
 
fprintf(fp,"POINT_DATA %d\n",(p+1)*(q+1));
fprintf(fp,"SCALARS Vorticity double\n");
fprintf(fp,"LOOKUP_TABLE default\n");

for (i=0;i<=p;i++)
{
 for (j=0;j<=q;j++)
 {
  if (i==0 && j==0)
   T_p[i][j] = w[i+1][j+1];
  else if (i == 0 && j !=0 && j!=q)
   T_p[i][j] = (w[i+1][j]+w[i+1][j+1])/2;
  else if (i == 0 && j == q)
   T_p[i][j] = w[i+1][j];
  else if (j ==0 && i!=0 && i!=p)
   T_p[i][j] = (w[i][j+1]+w[i+1][j+1])/2;
  else if (j ==0 && i == p)
   T_p[i][j] = w[i][j+1];
  else if (i ==p && j!=0 && j!=q)
   T_p[i][j] = (w[i][j]+w[i][j+1])/2;
  else if (i==p && j==q)
   T_p[i][j] = w[i][j];
  else if (j==q && i!=0 && i!=p)
   T_p[i][j] = (w[i][j]+w[i+1][j])/2;
  else
   T_p[i][j] = (w[i][j] + w[i+1][j+1] + w[i][j+1] + w[i+1][j])/4;
  fprintf(fp,"%24.16E\n",T_p[i][j]);
 }
}

  fclose(fp);

	
	
	return(0);
}
