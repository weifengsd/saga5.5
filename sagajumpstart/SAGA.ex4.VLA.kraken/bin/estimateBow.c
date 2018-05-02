/*  routine will output array shape with bow and tilt in XZ plane    */
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
 double Z0,dr,tilt,bow,pi; 
 double L,La,Lp;
 double x[64],xRot[64],z[64],zRot[64];

 pi = 4*atan(1.);

 if (argc < 6) {
    printf(" estimateBow nElem Z0 dr tilt bow\n");
    printf(" estimateBow  21 100. 5.62475 5.0 10.\n");
    printf("   where \n");
    printf("        nElem = number of elements\n");
    printf("           Z0 = depth of first element (m)\n");
    printf("           dr = distance between each element (m)\n");
    printf("         tilt = tilt of array (degrees)\n");
    printf("          bow = bow of array at midpoint (m) \n");

    exit(1);
   }


 nElem = atoi(argv[1]);
    Z0 = atof(argv[2]);
    dr = atof(argv[3]);
  tilt = atof(argv[4])*pi/180.; /* convert tilt from degrees to radians */
   bow = atof(argv[5]);

 if ( nElem > 64 ) { 
      printf(" estimateBow limited to 64 elements\n");
      exit(1);
     }
 
 /* compute length of straight array */
 L = (nElem-1) * dr; 

 /* compute element depths */
 for (i=0; i<nElem; i++ )
     z[i] = i * dr;

 /* compute length along z axis */
 Lp =L - (8/3)*pow(bow,2.)/L; 
 /* if ( dr < 0 ) Lp = Lp * -1; */

 /* compute new z axis locations */
 for (i=0; i<nElem; i++ ) {
     z[i] = z[i]*(Lp/L); }
 
 /* compute deflection along the x axis for each element (BOW) */
 for (i=0; i<nElem; i++ ) {
     x[i] = 4*bow * (z[i]/Lp - pow((z[i]/Lp),2.)); }

 /* adjust z axis positions to ensure that elements are */
 /*        exactly 'dr' meters apart along arc          */
 for (i=1; i<nElem; i++ ) {
     if ( dr < 0 ) {
        z[i] = z[i-1] - pow( pow(dr,2.) - pow((x[i]-(x[i-1])),2.),0.5); }
     else {
        z[i] = z[i-1] + pow( pow(dr,2.) - pow((x[i]-(x[i-1])),2.),0.5); }
    }

 /* rotate around Y axis (array tilt) */
 for (i=0; i<nElem; i++ ) {
     xRot[i] = x[i]* cos(tilt) + z[i]*sin(tilt);
     zRot[i] = x[i]*-sin(tilt) + z[i]*cos(tilt); 
    }
 
 /* for FIELD program want zaxis information first */
 /* then want deflection in xaxis */

 printf(" %2d",nElem);
 for (i=0; i<nElem; i++ )
     printf(" %7.3f",zRot[i]+Z0);
 printf("\n");

 printf(" %2d",nElem);
 for (i=0; i<nElem; i++ )
     printf(" %07.3f",xRot[i]);
 printf("\n");
} 
