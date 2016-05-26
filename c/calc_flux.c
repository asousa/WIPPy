#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include "consts.h"
#include "get_freqs.c"

/* A cleaned-up version of calcFlux.c.
 * Started 5/25/2016, Austin Sousa -- asousa@stanford.edu
 *  -Added code to get frequencies from directory listing
 *  -Deleted alpha file stuff (no resuming)
 */ 


// -----------------------------------------------
// GLOBAL VARIABLES:
// -----------------------------------------------

double      L_TARG;
int nancounter;


// -----------------------------------------------
// Constants to perform Gauss quadrature integration
// -----------------------------------------------
double t5[] =   {   -0.9061798459, 
            -0.5384693101, 
            0, 
            0.9061798459, 
            0.5384693101    };

double beta5[] = {  0.2369268851,
            0.4786286745, 
            0.56888888889, 
            0.2369268851,
            0.4786286745    };


// -----------------------------------------------
// PROTOTYPES
// -----------------------------------------------

float *getArr(void);
void updateArr(float *arr1, float *arr2);
void compFlux(float *arr, double L, int k, char *dir, char *flux_filename);
void readJ(float J[][100], char *filename);
double getJdiff(float J[][100], double E, double alpha_lc);



/* -----------------------------------------------
 * calc_flux Main
 *
 * - Calculates differential electron flux due to
 *   the deflection matrices (pN*, pS* files), based
 *   on an initial magnetosphere population (flux_filename).
 *   
 * Inputs:
 *  -dir:    Directory to folder of pN, pS files
 *  -L_TARG: L-shell to compute on
 *  -flux_filename: Path to the flux file (EQFLUXMA.dat)
 */
int main(int argc, char *argv[])
{
  FILE *inPtr, *alphaPtr;
  int numFreqs, i, k, m, nin, ei, ti;
  char *dir, sysCmd[512], filename[64], *NS, alphaFile[128];
  char *flux_filename;
  float L, *arr, *arrarr[128]; // <- array of pointers to precip arrays
  float J[100][100];
  double Jext;
  int *freqs;
  int num_freqs;

  nancounter = 0;
  // Get input parameters from command line
  dir   = argv[1];
  L_TARG = atof(argv[2]);
  flux_filename = argv[3];


  
  printf("DIRECTORY: %s\n",dir);
  printf("L_TARG: %f\n",L_TARG);
  printf("Flux file: %s\n",flux_filename);
  

  printf("\n\aWill do %d time steps\n", NUM_STEPS);


  freqs = get_freqs_at(dir, &num_freqs);

  printf("\n\aWill do %d frequencies\n",num_freqs);  

  
  // printf("\n\aWill do %d frequencies\n", numFreqs);
  //numFreqs = 1;
  sprintf( alphaFile, "%s/alpha_%g_%s", dir, L_TARG, "N" );
  printf("alphaFile: %s\n",alphaFile);
  arr = getArr();


  for(i=0; i<num_freqs; i++) {
    // printf("Frequency: %d \n",freqs[i]);
    for(k=0; k<2; k++) {

        NS = (k==0) ? "N" : "S";

          // Initialize the array (if we're the first)
          if(i==0)  arrarr[k]=getArr();

              sprintf(filename,"%s/p%s%d_%g.dat",dir,NS,freqs[i],L_TARG);
              // printf("i: %d, k: %d, filename: %s\n", i, k, filename);
              inPtr = fopen(filename, "r");
              printf("opened %s, pointer is %d\n",filename,inPtr);
              nin = fread(arr, sizeof(float), (NUM_E*NUM_STEPS), inPtr);
              fclose(inPtr);
      
              // Add to rolling sum from previous files
              updateArr( arrarr[k] , arr );

            } // for(k ... )  N/S - hemisphere

          } // freqs
  
  // Compute flux for north and south
  for(k=0; k<2; k++) {
    compFlux(arrarr[k], L_TARG, k, dir, flux_filename);
  } // N/S - hemisphere


  return 0;
}



/*
 * FUNCTION: compFlux
 * ------------------
 * This function will open a file with the appropriate hemisphere and
 * L-shell inserted into the name, calculate the precipitated flux  
 * due to each cell in the dAlpha_RMS matrix and write it into the 
 * file.
 *
 */

void compFlux(float *arr, double L, int k, char *dir, char *flux_filename)
{
  FILE *phiPtr, *QPtr, *NPtr, *alphaPtr;
  int ei, ti, i, nout;
  double mag, P, Q, epsm, alpha_eq, v, crunch;
  double I=0.0, x, g, field, Phi_p, b=1.0, vcm, fcgs;
  double v_tot_arr[NUM_E], E_tot_arr[NUM_E], gamma;
  double Jdiff[NUM_E];
  float Phi_float, alpha, J[100][100];
  char *NS, PhiFile[128], QFile[128], NFile[128], alphaFile[128];

  // float arg_in1;
  // float arg_in2;
  // float arg_in3;
  // float arg1;
  // float arg2;
  // float arg3;

  // Open up Phi file for writing
  if(k==0) {NS = "N";} else {NS="S";}  

  sprintf( PhiFile, "%s/phi_%g_%s", dir, L, NS );
  // sprintf( QFile,   "%s/Q_%g_%s", dir, L, NS );
  // sprintf( NFile,   "%s/N_%g_%s", dir, L, NS );
  // sprintf( alphaFile,   "%s/alpha_%g_%s", dir, L, NS );  

  printf("writing %s\n", PhiFile);
  if( (phiPtr=fopen(PhiFile, "w"))==NULL ) {
    printf("\nProblem opening %s\n", PhiFile);
   exit(0);
   }

  if( (alphaPtr=fopen(alphaFile, "w"))==NULL ) {
   exit(0);
  }


  epsm = (1/L)*(R_E+H_IONO)/R_E;

  crunch  = sqrt(1+3*(1-epsm))/pow(epsm,3) ;

  alpha_eq  = asin(sqrt( 1/crunch ));

  readJ(J, flux_filename);

  // Precalculate energy and velocity values
  for(i=0; i<NUM_E; i++) {
    E_tot_arr[i] = pow(10, (E_EXP_BOT+(DE_EXP/2)+DE_EXP*i) ); //E in eV
    Jdiff[i] = getJdiff( J, E_tot_arr[i], alpha_eq );
    v_tot_arr[i] = C*sqrt(1 - pow( (E_EL/(E_EL+E_tot_arr[i])) ,2) );
  }

  for(ei=0; ei<NUM_E; ei++) {

    if(ALPHA_DISTRIBUTION) {
      v = v_tot_arr[ei];
      vcm = v*100;  // v in cm for distrib fn calculation
      gamma = 1.0/sqrt( 1 - v*v/(C*C) );
      fcgs =  4.9e5/pow( (vcm*gamma) ,4) - 
    8.3e14/pow( (vcm*gamma) ,5) + 
    5.4e23/pow( (vcm*gamma) ,6);
     



      //fcgs = fcgs*50;
      // fcgs = 7.034e26 / pow(vcm,6); //10^8/E^2 distribution

      b = (v*v/M_EL)*pow( sqrt(1 - (v*v)/(C*C)), 3) * 1.6e-8 * fcgs;

      // b = 1e8 / pow(E_tot_arr[i],2);
    } else {
      b = Jdiff[ei]*1000;
    }

    //display values of b and E
    // USEFUL BUT ANNOYING RIGHT NOW - 5/2015
  //  printf("b = %f\tE = %f\n",b,E_tot_arr[ei]);

    for(ti=0; ti<NUM_STEPS; ti++) {
      
      alpha = sqrt( arr[ei*NUM_STEPS+ti] );

      

      nout=fwrite(&alpha, sizeof(float), 1, alphaPtr);      

      mag = 1.4142135623731*alpha;  // sqrt(2)*alpha_RMS = peak
      
      P = mag/2;    //[ alpha_lc - (alpha_lc-mag) ] /2
      Q = alpha_eq - mag/2; //[ alpha_lc + (alpha_lc-mag) ] /2
      
      I = 0.0;
      if(mag != 0) {
       for(i=0; i<5; i++) {
         x = P*t5[i] + Q ;

         // ------------- Austin's edits: Trying to fix NaNs being
         //               thrown by asin() function
   

          if(ALPHA_DISTRIBUTION) {



           g = (P/PI)*sin(2.0*x)*(  asin((x-alpha_eq)/mag)+ (PI/2.0) );

            } else {

           
            g = (P/PI)*sin(2.0*x)*((x - alpha_eq)*( asin((x-alpha_eq)/mag)+ (PI/2.0) ) +
               sqrt(mag*mag-pow((x-alpha_eq),2)));
          }; // Square

       I += ( beta5[i]*g );

       } // for(i ... ) -> Gauss quad integration
      } // if mag != 0
      Phi_p = PI*crunch*b*I;
      Phi_float = ( float ) Phi_p;

      //if isnan(I) { printf("I ISNAN\n"); };


      if isnan(Phi_p) { 
        nancounter=nancounter + 1;
        //printf("Total NaNs: %i\n",nancounter);
      };

      nout=fwrite(&Phi_float, sizeof(float), 1, phiPtr);
      
    } // for(ti ... )
  } // for(ei ... )
  
  fclose(phiPtr);
  fclose(alphaPtr);

  printf("Total NaNs: %i\n",nancounter);

  // Now calculate Q and N
  //
  // Need to integrate over E
  //for(ti=0; ti<NUM_TIMES; ti++) {
  //  for(ei=0; ei<NUM_E; ei++) {
  //    } // ei
  // }  // ti

}



/*
 * FUNCTION: getJdiff
 * ------------------
 * Using the AE8 data stored in J, calculate the differential flux 
 * by taking the (energy) derivative of the integral flux, dividing
 * by the total solid angle and extrapolating down in energy if 
 * need be.
 *
 */
double  getJdiff(float J[][100], double E, double alpha_lc)
{
  int row, i, topCol, botE;
  double J1, J2, I, x1, x2, y1, y2, m, c, x_ext, y_ext, J_ext;

  row = (int)floor((L_TARG+0.11 - J[1][0])/0.1); // to make sure! 
  
  // if(  fabs((double)J[row][0]-L_TARG) > 1e-3   ) 
  //   printf("\nL-shell not matching data\n\a");

  I = PI * cos(alpha_lc) * (PI - 2*alpha_lc);

  // Find column corresponding to highest energy value
  for(i=0; i<100; i++) {
    if(J[0][i+1] < 0.01) { 
      topCol = i; 
      break; 
    }
  }



  // Case 1. E < 100 keV
  // -------------------

  if( E <= 1e5 ) {
 
    // diff flux @ 100 keV and 200 keV
    J1 = 1e-6*fabs(J[row][2] - J[row][1]) / (J[0][2] - J[0][1]); 
    J2 = ((1e-6*fabs(J[row][3] - J[row][2]) / (J[0][3] - J[0][2])) 
      + J1 )/2; // central difference

    // do extrapolation in log-log space for best fit 
    x1 = log10( J[0][1]*1e6 );
    x2 = log10( J[0][2]*1e6 );
    y1 = log10( J1 );
    y2 = log10( J2 );

    m = (y2-y1)/(x2-x1);            // gradient of line
    c = (y1*x2 - y2*x1)/(x2-x1) ;   // offset of line, i.e.
    
    // y = m*x + c
    x_ext = log10( E );
    y_ext = m*x_ext + c;
    J_ext = pow(10, y_ext);

    return (J_ext/I);

  }


  

  // Case 2. E > 7 MeV
  // -----------------

  if( E >= 7e6 ) {
  
    // If flux at 7 Mev = 0, flux above it is zero too
    if( J[row][topCol]==0 )  return 0;

    // Otherwise need to extrapolate as in case 1.
    // diff flux @ 6.5 MeV and 7 MeV
    J2 = 1e-6*fabs( J[row][topCol] - J[row][topCol-1] ) 
      / (J[0][topCol] - J[0][topCol-1]); 

    J1 = ((1e-6*fabs( J[row][topCol-1] - J[row][topCol-2]) / 
       (J[0][topCol-1] - J[0][topCol-2]) ) + J2 )/2; // cdiff

    // do extrapolation in log-log space for best fit 
    x1 = log10( J[0][topCol-1]*1e6 );
    x2 = log10( J[0][topCol]*1e6 );
    y1 = log10( J1 );
    y2 = log10( J2 );

    m = (y2-y1)/(x2-x1);        // gradient of line
    c = (y1*x2 - y2*x1)/(x2-x1) ;   // offset of line, i.e.
                    // y = m*x + c
    x_ext = log10( E );
    y_ext = m*x_ext + c;
    J_ext = pow(10, y_ext);

    if(J_ext < 1e-10 ) J_ext = 0.0;

    return (J_ext/I);
  }


  // Case 3. 100 keV < E < 7 MeV
  if( E<7e6 && E>1e5 ) {


    // Find column corresponding lower energy value
    for(i=1; i<100; i++) {
      if( (J[0][i+1]*1e6) > E ) { 
    botE = i; 
    break; 
      }
    }


    // central diff flux @ lower and higher energies
    J1 = ( (1e-6 * fabs( J[row][botE] - J[row][botE-1] )
        / ( J[0][botE] - J[0][botE-1] ) ) + 
       (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
        / ( J[0][botE+1] - J[0][botE] ) )  ) / 2;

    J2 = ( (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
        / ( J[0][botE+1] - J[0][botE] ) ) + 
       (1e-6 * fabs( J[row][botE+2] - J[row][botE+1] )
        / ( J[0][botE+2] - J[0][botE+1] ) )  ) / 2;

    if(botE == 1)
      J1 =  (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
          / ( J[0][botE+1] - J[0][botE] ) );
    
    if(botE == (topCol-1))
      J2 = (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
        / ( J[0][botE+1] - J[0][botE] ) );
    



    // If J1 = J2 = 0, interpolated value also 0
    if( J1==0 && J2==0 ) return 0;



    // If only J2 = 0, do linear interpolation
    if( J2 == 0 ) {
      J_ext = J1*( ( J[0][botE+1]-(E*1e-6) )/
           ( J[0][botE+1] - J[0][botE] ) );
      return (J_ext/I);
    }



    // Otherwise interpolate as in case 1 (log-log space)

    x1 = log10( J[0][botE]*1e6 );
    x2 = log10( J[0][botE+1]*1e6 );
    y1 = log10( J1 );
    y2 = log10( J2 );

    m = (y2-y1)/(x2-x1);        // gradient of line
    c = (y1*x2 - y2*x1)/(x2-x1) ;   // offset of line, i.e.
                    // y = m*x + c
    x_ext = log10( E );
    y_ext = m*x_ext + c;
    J_ext = pow(10, y_ext);

    return (J_ext/I);
  }

}


// -----------------------------------------------
// UTILITIES
// -----------------------------------------------

/*
 * FUNCTION: readJ
 * ---------------
 * This function simply looks for a file 'filename' and reads 
 * it in.  The columns are energies and the rows are L-shells.
 * The first column is just a list of L-shells and the first row 
 * is just a list of energies (in MeV).
 *
 */
void readJ(float J[][100], char *filename)
{
  // char *filename;
  FILE *filePtr;
  int i,j;

  // filename = "EQFLUXMA.dat";

  if( (filePtr = fopen( filename ,"r")) == NULL ) {
    printf("Error opening the flux file! path: %s\n",filename);
    exit(1);
  }
  

  // INITIALIZE
  for(i=0; i<100; i++) {
    for(j=0; j<100; j++) {
      J[i][j] = 0.0;
    }
  }

  // READ IN VALUES
  for(i=0; i<47; i++) {
    for(j=0; j<31; j++) {
      fscanf(filePtr, "%e", &(J[i][j]));
    }
  }

  fclose(filePtr);
}




/* 
 * FUNCTION: getArr
 * ----------------
 * This function simply allocates dynamically a block of memory the 
 * size of NUM_E * NUM_STEPS of type float, and initializes it.
 * It returns the pointer to the memory.
 *
 */
float *getArr(void)
{
  float *arr;
  int ei, ti;
    printf("calling getArr...\n");

  arr = (float *) malloc( NUM_E * NUM_STEPS * sizeof(float) );
  if(arr == NULL) {
    printf("\nProb assigning mem in calcFlux\n");
    exit(0);
  }

  for(ei=0; ei<NUM_E; ei++) {
    for(ti=0; ti<NUM_STEPS; ti++) {
      arr[ei*NUM_STEPS+ti] = 0.0;
    }
  }
  //printf("Finishing getArr...\n");
  return arr;

}


/* 
 * FUNCTION: updateArr
 * -------------------
 * This function updates the values of arr1 with those of arr2.
 *
 */
void updateArr(float *arr1, float *arr2)
{
  int ei, ti;
  //printf("starting updateArr\n");
  for(ei=0; ei<NUM_E; ei++) {
    for(ti=0; ti<NUM_STEPS; ti++) {
      if(arr2[ei*NUM_STEPS+ti]>0.0)
        arr1[ei*NUM_STEPS+ti] += arr2[ei*NUM_STEPS+ti];
    }
  }
}
