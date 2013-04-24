#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tgmath.h>

#define OFT(i,j,mat) (i * mat.m + j)

typedef struct {
	double* context;
	int m;
	int n;
} matrix_t;

void setup(matrix_t *matrix) {
	// create 3x3 matrix
	matrix->m = matrix->n = 3;
	matrix->context = calloc(matrix->m * matrix->n, sizeof(double));
	if (matrix->context == NULL) {
		exit(1);
	}

	// fill that sucker with data
	double data[9] = {
		 4.0,   12.0, -16.0,
		 12.0,  37.0, -43.0,
		-16.0, -43.0,  98.0
	};
	memcpy(
		matrix->context,
		&data,
		sizeof(double) * matrix->m * matrix->n
	);

}

void teardown(matrix_t *matrix) {
	free(matrix->context);
	memset(matrix, 0, sizeof(matrix_t));
}

void dump_matrix(matrix_t *matrix) {
	for (int i=0; i < matrix->m*matrix->n; i++) {
		printf("%6.2f", matrix->context[i]);
		if (0 == (i + 1) % matrix->m) {
			printf("\n");
		} else {
			printf(" ");
		}
	}
	printf("\n");
}

void cholesky(matrix_t *A) {
	matrix_t G = { .m = A->m, .n = A->n };
	double Ljj, sum;
	int i,j,k;

	// create the output
	G.context = calloc(G.m * G.n, sizeof(double));
	memcpy(G.context, A->context, G.m * G.n * sizeof(double));

	// compute cholesky
	for (j=0; j < G.n; j++) {
		sum = 0;
		for(k=0; k < j; k++) {
			sum += G.context[OFT(j,k,G)] * G.context[OFT(j,k,G)];
		}
		G.context[OFT(j,j,G)] = sqrt(A->context[OFT(j,j,(*A))] - sum);
		for (i=0; i < G.m; i++) {
			sum = 0;
			for (k=0; k < j; k++) {
				sum += G.context[OFT(i,k,G)] * G.context[OFT(j,k,G)];
			}
			G.context[OFT(i,j,G)] = (1 / G.context[OFT(j,j,G)]) *
				(A->context[OFT(i,j,(*A))] - sum);
		}
	}

	dump_matrix(&G);
}

int main (int argc, char* argv[]) {
	matrix_t matrix;

	// setup matrix
	setup(&matrix);
	dump_matrix(&matrix);

	// do cholesky in place
	cholesky(&matrix);

	// teardown matrix
	teardown(&matrix);

	return 0;
}
