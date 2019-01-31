/**
 * convert result of convergence to matrices of L=50 in format=1 to matrices of L=30 and format=4
 */

/*
INPUT:
    stored matrices in a text file: "output.matrix"
    counts are there and the initial parameters
    
    Parameters:
	rho	- this is threshold on K (number of matches in a profile
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <errno.h>
#include <error.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "PSSM.h"

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define AVG(a, b)  ((a+b)/2.0)

const char *source_matrices_filename = "output.matrix";
const char *result_matrix_filename = "output_30_format4.matrix";
const int length = 30;
const int rho = 10;


void write_Freq_Matrices(const char *filename, double matrices[][50][26], char *CM, unsigned int counts[], int n) {
	char *amino_acids = "ACDEFGHIKLMNPQRSTVWY";
	FILE *fd = NULL;
	int m=0, i=0, j=0;

	if (n <= 0){
		return;
	}
	
	fd = fopen(filename, "w");
	assert(fd != NULL);
	
	for (m=0; m<n; m++) {
		if (CM != NULL && CM[m] != 'x') {
			continue;
		}
		fprintf(fd, "PROFILE 4\n");
		fprintf(fd, "BEGIN\n");
		fprintf(fd, "MATRIX ID=%d K=%d L=%d\n", m, counts[m], length);
	
		fprintf(fd, "%2d ", length);
		for (j=0; j<strlen(amino_acids); j++) {
			fprintf(fd, "       %c ", amino_acids[j]);
		}
		fprintf(fd, "\n");
		for (i=0; i<length; i++) {
			fprintf(fd, "%2d ", i);
			for (j=0; j<strlen(amino_acids); j++) {
	
				if (matrices[m][i][amino_acids[j] - 'A'] > 1.0) {
					fprintf(stderr, "Matrix %d is defected at position %d aminoacid %d (%c) with frequency = %1.6lf. Correcting to 1.0\n", 
						m, i, j, amino_acids[j], matrices[m][i][amino_acids[j] - 'A']);
					
					// correcting unexpectedly high frequency, resulting probably from unprecise double calculations on flancs.
					// ABC & DEF => (with overlap on B&E, ABC is more central) =>  A (B+E)/2 C => kABC+kDEF introduces an error on A and C 
					matrices[m][i][amino_acids[j] - 'A'] = 1.0;
				}
	    			
	    			fprintf(fd, "%1.6lf ", matrices[m][i][amino_acids[j] - 'A']);
				//assert(matrices[m][i][amino_acids[j] - 'A'] <= 1.0);
			}
			fprintf(fd, "\n");
		}
		fprintf(fd, "END\n");
	}
	fclose(fd);
}

// returns number of loaded matrices
int load_Matrices(const char *filename, double (**pmatrices)[50][26], unsigned int **pcounts, int rho_param) {
	int profile_format = 1;
	char segment[52] = "";
	char *buf = calloc(256, sizeof(char));
	int m=0, n=0, j=0;
	int pos=0, c[20];
	double b[20];
	bool discard_matrix = false;
	FILE *fd = NULL;
	Matrix M;

	if (*pmatrices != NULL) {
		printf("Input error. load_Matrices expects an unallocated pointer to matrices (NULL)\n");
		return 0;
	}

	fd = fopen(filename, "r");
	assert(fd != NULL);

	while (fgets(buf, 256, fd)) {
		if (1 == sscanf(buf, "PROTOTYPE %d", &profile_format)) {
			if (profile_format != 1 && profile_format != 2) {
				printf("Format error. Wrong prototype format version! Expecting version 1.\n");
				return 0;
			}
		}
		if (strstr(buf, "BEGIN")) {
			discard_matrix = true;
		}
		if (strstr(buf, "END")) {
			if (!discard_matrix) {
				m++;
			}
		}
		if (profile_format == 1) {
			if (1 == sscanf(buf, "MATRIX K=%d", &M.K)) {
				if (M.K < rho_param) {
					discard_matrix = true;
				} else {
					discard_matrix = false;
				}
			}
		}
		if (profile_format == 2) {
			if (2 == sscanf(buf, "MATRIX ID=%d K=%d", &M.id, &M.K)) {
				if (M.K < rho_param) {
					discard_matrix = true;
				} else {
					discard_matrix = false;
				}
			}
		}
	}
		
	*pmatrices = (double (*)[50][26])calloc(2*m, sizeof(double)*50*26);
	assert(m == 0 || *pmatrices != NULL);
	*pcounts = (unsigned int *)calloc(2*m, sizeof(unsigned int));
	assert(m == 0 || *pcounts != NULL);

	rewind(fd);
	
	while (fgets(buf, 256, fd)) {
		if (1 == sscanf(buf, "PROTOTYPE %d", &profile_format)) {
			if (profile_format != 1 && profile_format != 2) {
				printf("Format error. Wrong prototype format version! Expecting version 1.\n");
				return 0;
			}
		}
		if (strstr(buf, "BEGIN")) {
			discard_matrix = true;
		}
		if (strstr(buf, "END")) {
			if (!discard_matrix) {
				memcpy((*pmatrices)[n], &M.freq, sizeof(double[50][26]));
				(*pcounts)[n] = M.K;
				n++;
			}
		}
		if (1 == sscanf(buf, "SEGMENT %s", segment)) {
			strncpy(M.initial_segment, segment, 50);
		}

		if (profile_format == 1) {
			if (1 == sscanf(buf, "MATRIX K=%d", &M.K)) {
				if (M.K < rho_param) {
					discard_matrix = true;
				} else {
					discard_matrix = false;
				}
				M.id = n;
			}
			if (21 == sscanf(buf,
			    "%2d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", 
			     &pos,&c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],&c[8],&c[9],&c[10],&c[11],&c[12],&c[13],&c[14],&c[15],&c[16],&c[17],&c[18],&c[19])) {
		    	
				if (pos<0 || pos>50-1) {
					printf("Format error: invalid position %d\n", pos);
					continue;
				}
				for (j=20; j--;) {
					if (c[j]<0) {
						printf("Format error: invalid position %d,%d = %d\n", pos, j, c[j]);
					}

					M.freq[pos][amino_acids[j] - 'A'] = (double)(c[j]) / M.K;
				}
			}
		}
		if (profile_format == 2) {
			if (2 == sscanf(buf, "MATRIX ID=%d K=%d", &M.id, &M.K)) {
				if (M.K < rho_param) {
					discard_matrix = true;
				} else {
					discard_matrix = false;
				}
			}
			if (21 == sscanf(buf,
			    "%2d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
			     &pos,&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8],&b[9],&b[10],&b[11],&b[12],&b[13],&b[14],&b[15],&b[16],&b[17],&b[18],&b[19])) {
		    	
				if (pos<0 || pos>49) {
					printf("Format error: invalid position %d\n", pos);
					continue;
				}
				for (j=20; j--;) {
					if (b[j]<0) {
						printf("Format error: invalid position %d,%d = %lf\n", pos, j, b[j]);
					}

					M.freq[pos][amino_acids[j] - 'A'] = b[j];
				}
			}
		}
	}
	
	fclose(fd);
	free(buf);
	return m;
}


int main (int argc, char *argv[]) {

	double (*matrices)[50][26] = NULL;
	double (**pmatrices)[50][26] = &matrices;
	unsigned int *counts = NULL;
	unsigned int **pcounts = &counts;
	int N = 0;
	
	N = load_Matrices(source_matrices_filename, pmatrices, pcounts, rho);
	#ifndef DNEBUG
	fprintf(stderr, "With rho=%d loaded %d matrices\n", rho, N);
	#endif
		
	if (N > 0) {
	    write_Freq_Matrices(result_matrix_filename, matrices, NULL, counts, N);
	    #ifndef DNEBUG
		fprintf(stderr, "Written profiles in format 4 to file: %s\n", result_matrix_filename);
	    #endif
	}
}