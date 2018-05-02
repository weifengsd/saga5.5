/*  routine will output array shape with bearing and bow in XY plane    */
/*  for more information see chapter 13 of SAGA manual               */

/*  compile: cc -o estimateBow estimateBow.c -lm */

#include <stdio.h>
#include <math.h>

main(argc,argv)
int argc;
char *argv[];
{
 double atof(),atan(),cos(),sin(),pow();
 int nElem,i;
 double Y0,depth,dr,bearing,bow,pi; 
 double L,La,Lp;
 double x[64],xRot[64],y[64],yRot[64];

 pi = 4*atan(1.);

 if (argc < 7) {
    printf(" estimateHBow nElem Y0 depth dr bearing bow\n");
    printf(" estimateHBow  27 0. 213. 7.89125 34.5 15.\n");
    printf("   where \n");
    printf("        nElem = number of elements\n");
    printf("           Y0 = Y axis coordinate of first element (m)\n");
    printf("           Z0 = Depth of horizontal array (m)\n");
    printf("           dr = distance between each element (m)\n");
    printf("         bearing = bearing of array (clockwise from North)\n");
    printf("          bow = bow of array at midpoint (m) \n");
    printf("          \n");
    printf("   output:       \n");
    printf("          array element positions are output in the \n");
    printf("          format expected by mplfield.  first number of elements\n");
    printf("          then xyz coordinates for each element \n");

    exit(1);
   }


    nElem = atoi(argv[1]);
       Y0 = atof(argv[2]);
    depth = atof(argv[3]);
       dr = atof(argv[4]);
  bearing = atof(argv[5])*pi/180.; /* convert bearing from degrees to radians */
      bow = atof(argv[6]);

 if ( nElem > 64 ) { 
      printf(" estimateBow limited to 64 elements\n");
      exit(1);
     }
 
 /* compute length of straight array */
 L = (nElem-1) * dr; 

 /* compute element locations on Y axis*/
 for (i=0; i<nElem; i++ )
     y[i] = i * dr;

 /* compute length along y axis */
 Lp =L - (8/3)*pow(bow,2.)/L; 
 /* if ( dr < 0 ) Lp = Lp * -1; */
 /* printf(" L: %f   Lp: %f\n",(float) L,(float) Lp); */

 /* compute new y axis locations */
 for (i=0; i<nElem; i++ ) {
     y[i] = y[i]*(Lp/L); 
   /*  printf(" y[%2d]: %5.2f = y[%2d]: %5.2f * ( %5.2f / %5.2f) \n",  */ 
   /*                 i,y[i],i,y[i]*(L/Lp),Lp,L);                      */
    }
 
 /* compute deflection along the x axis for each element (BOW) */
 for (i=0; i<nElem; i++ ) {
     x[i] = 4*bow * (y[i]/Lp - pow((y[i]/Lp),2.)); }

 /* adjust y axis positions to ensure that elements are */
 /*        exactly 'dr' meters apart along arc          */
 for (i=1; i<nElem; i++ ) {
     if ( dr < 0 ) {
        y[i] = y[i-1] - pow( pow(dr,2.) - pow((x[i]-(x[i-1])),2.),0.5); }
     else {
        y[i] = y[i-1] + pow( pow(dr,2.) - pow((x[i]-(x[i-1])),2.),0.5); }
    }
 /* add in Y0 offset */
 for (i=0; i<nElem; i++ ) {
     y[i] = y[i] + Y0;
    }
 

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
