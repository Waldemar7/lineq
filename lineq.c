#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXDIM 4095
#define MAXDIM1 4096


int maxrows;
int maxcols;

double mymatrix[MAXDIM][MAXDIM1];


void printmatrix() {
	int row, col;
	for (row=0; row<maxrows; ++row) {
		for (col=0; col<maxcols; ++col) {
			printf("%f\t", mymatrix[row][col]);
		}
		printf("\n");
	}
}


void randomfill() {
	int row, col;
	for (row=0; row<maxrows; ++row) {
		for (col=0; col<maxcols; ++col) {
			mymatrix[row][col] = (double)rand()*(double)rand()/((double)rand()+1);
		}
	}
}


void fix_wiki_example_1() {
	maxrows = 3;
	maxcols = maxrows+1;
	mymatrix[0][0] = 1;
	mymatrix[0][1] = 3;
	mymatrix[0][2] = -2;
	mymatrix[0][3] = 5;
	mymatrix[1][0] = 3;
	mymatrix[1][1] = 5;
	mymatrix[1][2] = 6;
	mymatrix[1][3] = 7;
	mymatrix[2][0] = 2;
	mymatrix[2][1] = 4;
	mymatrix[2][2] = 3;
	mymatrix[2][3] = 8;
}


void fix_wiki_example_2() {
	maxrows = 3;
	maxcols = maxrows+1;
	mymatrix[0][0] = 2;
	mymatrix[0][1] = 1;
	mymatrix[0][2] = -1;
	mymatrix[0][3] = 8;
	mymatrix[1][0] = -3;
	mymatrix[1][1] = -1;
	mymatrix[1][2] = 2;
	mymatrix[1][3] = -11;
	mymatrix[2][0] = -2;
	mymatrix[2][1] = 1;
	mymatrix[2][2] = 2;
	mymatrix[2][3] = -3;
}


void solve(int do_log) {

	int row, col;
	int subrow;
	double factor, divis;

	long int starttime = clock();
	if (do_log) printf("Solve - Start:   %ld\n", clock());

	for (row=0; row<maxrows; ++row) {

		// TODO: check if divis is 0
		divis = mymatrix[row][row];
		for (col=row; col<maxcols; ++col) {
			mymatrix[row][col] /= divis;
		}

		for (subrow=row+1; subrow<maxrows; ++subrow) {
			factor = -(mymatrix[subrow][row]/mymatrix[row][row]);
			for (col=0; col<maxcols; ++col) {
				// printf("-> %d %d %d >%f< %f %f\n", row, subrow, col, factor, mymatrix[subrow][row], mymatrix[row][row]);
				mymatrix[subrow][col] += factor*mymatrix[row][col];
			}
		}
		//printf("Row: %d\n", row);
		//printmatrix();

	}

	if (do_log) printf("Solve - Echelon: %ld\n", clock());

	for (row=maxrows-1; row>0; --row) {

		for (subrow=row-1; subrow>=0; --subrow) {

			factor = -(mymatrix[subrow][row]/mymatrix[row][row]);
			for (col=0; col<maxcols; ++col) {
				// printf("-> %d %d %d >%f< %f %f\n", row, subrow, col, factor, mymatrix[subrow][row], mymatrix[row][row]);
				mymatrix[subrow][col] += factor*mymatrix[row][col];
			}
		}

		// printf("Row: %d\n", row);
		// printmatrix();

	}

	if (do_log) printf("Solve - Done:    %ld\n", clock());
	if (do_log) printf("Solve - Time:    %ld\n", clock()-starttime);

}


void solve_opt(int do_log) {

	int row, col;
	int subrow;
	double factor, divis;

	long int starttime = clock();
	if (do_log) printf("Solve - Start:   %ld\n", clock());

	for (row=0; row<maxrows; ++row) {

		// TODO: check if divis is 0 - and do something different
		divis = mymatrix[row][row];
		mymatrix[row][row] = 1; // optimization: set to 1 instead of would be dividing by itself
		for (col=row+1; col<maxcols; ++col) {
			mymatrix[row][col] /= divis;
		}

		for (subrow=row+1; subrow<maxrows; ++subrow) {
			factor = -(mymatrix[subrow][row]/mymatrix[row][row]);
			mymatrix[subrow][row] = 0; // optimization: normal calculation will set to 0, so we just do it
			for (col=row+1; col<maxcols; ++col) {
				// printf("-> %d %d %d >%f< %f %f\n", row, subrow, col, factor, mymatrix[subrow][row], mymatrix[row][row]);
				mymatrix[subrow][col] += factor*mymatrix[row][col];
			}
		}
		//printf("Row: %d\n", row);
		//printmatrix();

	}

	if (do_log) printf("Solve - Echelon: %ld\n", clock());

	for (row=maxrows-1; row>0; --row) {

		for (subrow=row-1; subrow>=0; --subrow) {
			factor = -(mymatrix[subrow][row]/mymatrix[row][row]);
			// optimization: only last column has to be calculated, , rest we set to zero
			/* for (col=row; col<(maxcols-1); ++col) {
				// printf("-> %d %d %d >%f< %f %f\n", row, subrow, col, factor, mymatrix[subrow][row], mymatrix[row][row]);
				// mymatrix[subrow][col] += factor*mymatrix[row][col];
				mymatrix[subrow][col] = 0;
			} */
			mymatrix[subrow][row] = 0;
			mymatrix[subrow][maxcols-1] += factor*mymatrix[row][maxcols-1];
		}

		// printf("Row: %d\n", row);
		// printmatrix();

	}

	if (do_log) printf("Solve - Done:    %ld\n", clock());
	if (do_log) printf("Solve - Time:    %ld\n", clock()-starttime);

}


int main(int argc, char **argv)
{
	int myseed = 7;

	printf("CLOCKS_PER_SEC is %ld\n", CLOCKS_PER_SEC);
	printf("RAND_MAX is %d\n", RAND_MAX);


	printf("\n--- Wikipedia Exmaple 1 (linear equation system) ---\n");
	fix_wiki_example_1();
	printf("Initial:\n");
	printmatrix();
	solve(0);
	printf("Solved:\n");
	printmatrix();

	fix_wiki_example_1();
	solve(0);
	printf("Solved-Optimized:\n");
	printmatrix();

	printf("\n--- Wikipedia Exmaple 2 (gaussian elimination) ---\n");
	fix_wiki_example_2();
	printf("Initial:\n");
	printmatrix();
	solve(0);
	printf("Solved:\n");
	printmatrix();

	fix_wiki_example_2();
	solve(0);
	printf("Solved-Optimized:\n");
	printmatrix();

	printf("\n-------------- SMALL (check output) ----------------\n");

	maxrows = 8;
	maxcols = maxrows+1;

	printf("Start:\t%ld\n", clock());

	srand(myseed);
	randomfill();
	// printf("Rand:\t%ld\n", clock());

	printf("Initial:\n");
	printmatrix();

	solve(1);
	printf("Solved:\n");
	printmatrix();

	srand(myseed);
	randomfill();

	solve_opt(1);
	printf("Solved-Optimized:\n");
	printmatrix();


	printf("\n-------------- LARGE (no output) ----------------\n");

	maxrows = 1500;
	maxcols = maxrows+1;

	printf("Start:\t%ld\n", clock());

	srand(myseed);
	randomfill();
	printf("Initial:\n");
	//printmatrix();

	solve(1);
	printf("Solved:\n");

	srand(myseed);
	randomfill();
	solve_opt(1);
	printf("Solved-Optimized:\n");




	printf("Done:\t%ld\n", clock());

	//printf("Waiting ...\n");
	//getchar();

	return 0;
}
