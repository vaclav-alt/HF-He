#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FORNAME(x) x##_
#define MAXIT 100


/* Declarations for the fortran interface */
// FORNAME just adds underscore to the fortran function names
void FORNAME(define_grid)();
void FORNAME(define_spline)();
void FORNAME(rkb_alloc)(int *);
void FORNAME(rkb_cell)();
double FORNAME(over1)(int *i1, int *l1, int *m1, int *i2, int *l2, int *m2);
double FORNAME(dint1)(int *i1, int *l1, int *m1, int *i2, int *l2, int *m2, double *Z);
double FORNAME(dint2)(int *i1, int *l1, int *m1, int *i2, int *l2, int *m2, 
                      int *i3, int *l3, int *m3, int *i4, int *l4, int *m4);
int FORNAME(getnbr)();

// pomocne funkce
// printArray vypise 1D pole
void printArray(double *pole, int size);
// prints out a 2D array, needs to be called like this: printMatrix(array[0], rows, cols)l
void printMatrix(double *pole, int rows, int cols);
// copies the content of one 2D array into another
// square matrices of size "size" are assumed
// contains no safety checks
void copyMatrix(double *from, double *to, int size);

/* End of declarations */

int main()
{
int LMAX2 = 2*0;
// dummy variables
int i1,l1,m1,ib1, i2,l2,m2,ib2, i3,l3,m3,ib3, i4,l4,m4,ib4, NB;
double ZNUC;
int i, j, k, l, run=0;

/* sets up grid points and initializes the values of the spline */
FORNAME(define_grid)();
FORNAME(define_spline)();
int NBR = FORNAME(getnbr)();

// we don't want the first and the last spline
// which is the dimension for our eigenvalue problem
NBR -= 2;
NB = NBR * (LMAX + 1)**2;

// dummy variables neede for LAPACK's DSYGV
int info =0, ITYPE = 1, LWORK;
char JOBZ='V', UPLO='U';
double W[NBR], help[1];
double *WORK;

// variables needed for calculation
double energy, energy2, energy_old, delta,
    prec=1e-12;                // required precision
double h[NBR][NBR];                 // one electron hamiltonian
double S[NBR][NBR];                 // the overlap matrix
double S_back[NBR][NBR];            // it's copy
double DAOI[NBR][NBR][NBR][NBR];    // two electron integrals
double H[NBR][NBR];                 // full hamiltonian
double J[NBR][NBR];                 // J operator matrix
double K[NBR][NBR];                 // K operator matrix
double Hen[NBR][NBR];               // operator for energy calculation
double vecs[NBR][NBR];              // matrix for storing of the eiegenvectors

/* allocate and calculate 2-electron radial integrals
      using Cell integral algorithm */
FORNAME(rkb_alloc)(&LMAX2);
FORNAME(rkb_cell)();

/* test for 1- and 2-electron elements
   bra vector specified by i1,l1,m1 
   ket vector specified by i2,l2,m2 
*/
i1=1; l1=0; m1=0; i2=2; l2=0; m2=0; ZNUC=2.0;
l3=0; m3=0; l4=0; m4=0;

// preparation of th oneElInt and overlap matrices
// printf("Calculating one- and two-electron integrals...");
for (i1=1; i1<=NBR; ++i1) {
	for (i2=1; i2<=NBR; ++i2) {
		// one electron matrix elements
		H[i1-1][i2-1] = FORNAME(dint1)(&i1,&l1,&m1, &i2,&l2,&m2, &ZNUC);
		// its copy
		h[i1-1][i2-1] = H[i1-1][i2-1];
		// overlap matrix
		S[i1-1][i2-1] = FORNAME(over1)(&i1,&l1,&m1, &i2,&l2,&m2);
		// its copy
		S_back[i1-1][i2-1] = S[i1-1][i2-1];

        // two more nested cycles two fill the NBRxNBRxNBRxNBR array DAOI
        // which contains the two electron integrals
		for (i3=1; i3<=NBR; i3++) {
			for(i4=1; i4<=NBR; i4++) {
				// two electron matrix elements
				DAOI[i1-1][i2-1][i3-1][i4-1] = FORNAME(dint2)(&i1,&l1,&m1,&i2,&l2,&m2,&i3,&l3,&m3,&i4,&l4,&m4);
			}
		}
	}
}
// printf("done\n\n");

// do-while cycle - is done at least once. At the end, if the condition in " while(condition)"
// is satisfied, it runs again
do {

// count this particular run
run++;
// dry run of DSYGV - LWORK estimation
// printf("DIAGONALIZATION OF H\n");
// printf("dry run - LWORK estimate...");
LWORK=-1;
DSYGV(&ITYPE, &JOBZ, &UPLO,&NBR, H, &NBR, S, &NBR, W, help, &LWORK, &info);
if (info==0) {
	// in case of success prints out result
	// printf("done\n");
	LWORK = help[0];
	// printf("\tEstimated LWORK: %d\n", LWORK);
}
// printf("\n");

// allocation of array WORK with size LWORK
// printf("\tAllocating WORK array with size LWORK: %d...", LWORK); 
WORK = (double *)malloc(LWORK * sizeof(double));
// printf("done\n");

// matrix diagonalization
// printf("Diagonalizing\n");
DSYGV(&ITYPE, &JOBZ, &UPLO,&NBR, H, &NBR, S, &NBR, W, WORK, &LWORK, &info);
if (info == 0) {
	// in case of success prints out lowest eigenvalue
	// printf("\tDiagonalization finished.\n\tLowest eigenvalue: %f\n", W[0]);
	// prints out all eigenvalues
	// printArray(W,NBR);
} else {
	// printf("\tDSYGV error number %d\n", info);
}

// dump WORK memory space
// printf("Freeing WORK memory space\n\n");
free(WORK);

// restore overlap from backup, save eigenvectors to vecs
// printf("Making copies of matrices...");
copyMatrix(S_back[0],S[0],NBR);
copyMatrix(H[0],vecs[0],NBR);
// printf("done\n\n");

// construction of matrices J, K, H according to formulae presented during the lecture
// printf("Constructing matrices J, K and H...");
for (i=0; i<NBR; i++) {
	for(j=0; j<NBR; j++) {
		J[i][j]=0;
		K[i][j]=0;
		for (k=0; k<NBR; k++) {
			for (l=0; l<NBR; l++) {
                // singlet
				J[i][j] += DAOI[i][j][k][l] * vecs[0][k] * vecs[0][l];
				K[i][j] += DAOI[i][l][k][j] * vecs[0][k] * vecs[0][l];
                // triplet
				// J[i][j] += DAOI[i][j][k][l] * ( vecs[0][k] * vecs[0][l] + vecs[1][k] * vecs[1][l]);
				// K[i][j] += DAOI[i][l][k][j] * ( vecs[0][k] * vecs[0][l] + vecs[1][k] * vecs[1][l]);
			}
		}
        // singlet
		H[i][j] = h[i][j]+2.*J[i][j]-K[i][j];
		Hen[i][j] = H[i][j]+h[i][j];
        // triplet
		// H[i][j] = h[i][j]+J[i][j]-K[i][j];
        // Hen[i][j] = h[i][j]+0.5*(J[i][j]-K[i][j]);
	}
}
// printf("done\n\n");

// energy calculated from c^T.Hel.c - multiplying Hel by eigenvector, which belongs
// to the lowest eigenvalue.
// BEWARE - these eigenvectors are produced by a FORTRAN procedure and are therefore
// save in the ROWS of the vecs array, not the columns
// printf("Calculating energy...");
energy_old = energy;
energy = 0;
energy2 = 0;
for (i=0; i<NBR; i++) {
	for (j=0; j<NBR; j++) {
        // energy2 += vecs[0][i] * h[i][j] * vecs[0][j];
        // triplet
		// energy += Hen[i][j] * (vecs[0][i] * vecs[0][j] + vecs[1][i] * vecs[1][j]);
        // singlet
		energy += Hen[i][j] * vecs[0][i] * vecs[0][j];
	}
}
energy2 += W[0];

//absolute value of energy progress
delta=fabs(energy-energy_old);
// printf("done\n");
//printf("Energy: %.8f\n",energy);
// printf("it. %2d: E1=%.12f\tE2=%.12f\tDIFF=%e\n", run, energy, energy2,fabs(energy-energy2));
printf("it. %2d: E1=%.12f\n", run, energy);
// printf("_____________________\n\n");

// if precision not reached or iteration count not exceeded, another iteration starts
} while (run <= MAXIT && delta >= prec);

 printf("\n======================\n");
 printf("delta: %.4e a prec: %.4e\n", delta, prec);
printf("Calculation stopped after %d iterations.\n",run);
printf("final energy=%.8f (a.u.)=%.8f eV\n", energy, energy*27.211);

// eigenvectors are saved in array vecs

// this verifies, that the indexth eiegen vector is an eigenvector
// index = 0;
// WORK = (double *)malloc(2* NBR * sizeof(double));
// for (i = 0; i < NBR; ++i) {
// 	WORK[i] = 0;
// 	WORK[NBR+i] = 0;
// 	// maticove nasobeni
// 	for (j = 0; j < NBR; ++j) {
// 		// multiplies the eigenvector with matrix H
// 		WORK[i] += H[i][j] * vecs[index][j];
// 		// multiplies the same eigenvector with the overlap matrix and the relevant eigenvalue
// 		WORK[NBR+i] += W[index]* S[i][j] * vecs[index][j];
// 	}
// }
// compares LHS and RHS of equation hC=eSC
// printf("hC: %f a eSC: %f\n", WORK[0], WORK[NBR+0]);
// printf("difference: %f\n", WORK[0]-WORK[NBR+0]);
// printf("eigenvalue: %f\n", W[index]);
// free(WORK);


// KONEC MAIN
}

void printArray(double *pole, int size) {
	int i;
	for (i=0; i<size; i++) {
		printf("%d: %f\n", i+1, *(pole + i));
	}
}


void printMatrix(double *pole, int rows, int cols) {
	int i,j;
	for (i = 0; i < rows; ++i) {
		for (j = 0; j < cols; ++j) {
			printf("%f ", *(pole + i*cols + j));
			// printf("%.3f ", pole[i][j]);
		}
		printf("\n");
	}
}

// void copyMatrix(double *from, double *to, int size) {
// 	memcpy(to,from,size*size*sizeof(double));
// }
void copyMatrix(double *from, double *to, int size) {
	int i, j, off;
	for (i=0; i<size; i++) {
		for (j=0; j<size; j++) {
			off = i*size +j;
			*(to + off)=*(from + off);
		}
	}
}
