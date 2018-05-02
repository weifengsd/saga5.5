/*  routine will output array shape with bearing and bow in XY plane    */
/*  for more information see chapter 13 of SAGA manual               */

/*  compile: cc -o estimateHBowf estimateHBowf.c -lsio -lm */


#include <stdio.h>
#include <math.h>

main(argc,argv)
int argc;
char *argv[];
{
 double atof(),atan(),cos(),sin(),pow();
 int nElem,i,np,nc,rl,ih[64],ic[64],debug=0;
 double Y0,depth,bearing,bow,pi; 
 double L,La,Lp;
 double x[64],xRot[64],y[64],yRot[64],dr[64];
 float a[64];

 pi = 4*atan(1.);

 if (argc < 5) {
    printf(" estimateHBow bearing bow depth infile\n");
    printf(" estimateHBow  34.5 15. 213. array.sio\n");
    printf("   where \n");
    printf("         bearing = bearing of array (clockwise from North)\n");
    printf("          bow = bow of array at midpoint (m) \n");
    printf("         array.sio = sio data file containing array coordinates\n");
    printf("                     file is single channel\n");
    printf("           \n");
    printf("   output:       \n");
    printf("          array element positions are output in the \n");
    printf("          format expected by mplfield.  first number of elements\n");
    printf("          then xyz coordinates for each element \n");

    exit(1);
   }


  bearing = atof(argv[1])*pi/180.; /* convert bearing from degrees to radians */
      bow = atof(argv[2]);
    depth = atof(argv[3]);
  
  bearing = -bearing;  /* negate bearing so positive increase goes clockwise */
  ncsetup(&np,&nc,ic,argv[4]);
  nElem = np;
  if (debug) printf(" np: %d  nc: %d\n",np,nc);
  if ( nc != 1 ) {
     printf("Error: array file must be a single channel\n");
     exit(1);
    }
  rdsio(a,np,nc,1,ic,&rl,argv[4],ih);

  if (debug) {
     for (i=0; i < np; i++ ) 
         printf("a[%d] = %f\n",i,a[i]);
        }

 if ( np > 64 ) { 
      printf(" estimateHBowf limited to 64 elements\n");
      exit(1);
     }

 /* compute length of straight array */
 L = a[np-1] - a[0];
 if (debug) printf(" a[%d]: %f - a[0]: %f = L: %f  \n",np,a[np-1],a[0],a[np-1]-a[0]);
 if (debug) printf(" L: %f  \n",(float) L);

 /* load element positions into Y axis */
 for (i=0; i<nElem; i++ ) {
     x[i] = 0.;
     y[i] = a[i];
    dr[i] = a[i+1] - a[i];
    if (debug) printf(" dr[%2d]: %5.2f = %5.2f - %5.2f \n",
                        i,(float) dr[i],(float) a[i+1],(float) a[i]);
   }

 /* compute length along y axis */
 Lp =L - (8/3)*pow(bow,2.)/L; 
 if (debug) printf(" Lp: %f  \n",(float) Lp);
 /* if ( dr < 0 ) Lp = Lp * -1; */

 /* compute new y axis locations */
 for (i=0; i<nElem; i++ ) {
     y[i] = y[i]*(Lp/L); }
 
 /* compute deflection along the x axis for each element (BOW) */
 for (i=0; i<nElem; i++ ) {
     x[i] = 4*bow * (y[i]/Lp - pow((y[i]/Lp),2.)); 
     if (debug) printf(" x[%2d]: %5.2f  y[%2d]: %5.2f\n",i,x[i],i,y[i]);
    }

 /* adjust y axis positions to ensure that elements are */
 /*        exactly 'dr' meters apart along arc          */
 for (i=1; i<nElem; i++ ) {
     if ( dr[i-1] < 0 ) {
        y[i] = y[i-1] - pow( pow(dr[i-1],2.) - pow((x[i]-(x[i-1])),2.),0.5); }
     else {
     if (debug) printf(" y[%2d]:%5.2f = y[%2d]:%5.2f + sqrt(%5.2f + %5.2f))\n",
               i,y[i],i-1,y[i-1],dr[i-1],x[i]-x[i-1]);
        y[i] = y[i-1] + pow( pow(dr[i-1],2.) - pow((x[i]-(x[i-1])),2.),0.5); }
    } 
 /* /* add in Y0 offset */
 /* for (i=0; i<nElem; i++ ) {
 /*     y[i] = y[i] + Y0;
 /*    } */
 

 /* rotate around Z axis (array bearing) */
 for (i=0; i<nElem; i++ ) {
     xRot[i] = x[i]*cos(bearing) - y[i]*sin(bearing);
     yRot[i] = x[i]*sin(bearing) + y[i]*cos(bearing); 
    }
 
 /* for mplFIELD program want x,y,z axis information*/

 printf("%2d\n",nElem);
 for (i=0; i<nElem; i++ )
     printf("           %9.2f %9.2f %9.2f\n",xRot[i],yRot[i],depth);
} 
