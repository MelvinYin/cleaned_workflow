/**
 *    Search PSSM procedure
 *    Coded in C programming language (C99 standard) and MPI
 *    Copyright (c) 2009 Alexandr Goncearenco
 *    Affiliation: CBU, BCCS, UNIFOB AS, UiB, Berezovsky Group.
 */

/*
INPUT:
    extracted.matrix: matrices in frequency format. the needles
    subject: the haystack. SCOP for example
    search parameters (hardcoded): p-value (or E value recalculated on subject): E = p*N

OUTPUT:
    M id \tab score.float \tab full fasta description of the match

WARNING:
    matrices are numbered according to position in extracted.matrix file

///////////////////////////

    read: matrices into memory from frequency format (2) -> F
    read: subject fasta sequences into memory
    calculate composition from subject sequences
    for each F:
	M = calcPSSM(F, K, composition)
	runP(M, permutation_vector, p-value) return S(p-value)
	runS(M, S(p)):
	    report each match to result file

 DONE
///////////////////////////
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <errno.h>
#include <error.h>

#include "PSSM.h"
 
// parameters
//
// with N = 20 M: p= 5^10-6 => E=10. rho = 10. epsilon=3
//
const int delta = 10;		// cutting step for initial segments, in residues
const double omega = 1.0;	// pseudocount scaling factor in PSSM
//const double omega = 2.0;	// pseudocount scaling factor in PSSM
//const double p = 0.00000072;	// initial p-value, lower threshold
//const double p = 0.000005;	// initial p-value, lower threshold
const double rho = 10;		// threshold ratio shuffled/natural number of matches where score > score(p)
const double mu = 0.05;		// upper threshold for p-value
const double eta = 0.01;	// increment of p-value in clustering
const double epsilon = 3;	// max distance(M1, M2) to consider M2 as converged
const int theta = 30;		// max # iterations

const char *extracted_matrix_filename = "extracted.matrix";
const char *subject_filename = "subject.fasta";
const char *control_filename = "control.fasta";
const char *results_filename = "search_matches.tab";

int load_fasta_sequences(const char *filename, Sequence **psequences, bool removex) {
	FILE *fd = NULL;
	char *seq = NULL;
	char desc[256];
	char buf[256];
	char *line = NULL;
	int N = 0;
	int i = 0;

	fd = fopen(filename, "r");
	if (fd == NULL) {
		error(1, errno, "error while opening file %s", filename);
	}
	assert(fd != NULL);

	for (;;) {
		line = fgets(buf, 256, fd);
		//if (line[0] == 0) continue;
		
		if (line != NULL) {
			if (line[strlen(line)] == 0x0A || line[strlen(line)] == 0x0D)
				line[strlen(line)] = 0x00;
			if (strlen(line) > 1 && (line[strlen(line)-1] == 0x0A || line[strlen(line)-1] == 0x0D))
				line[strlen(line) - 1] = 0x00;
		}
		
		if (line == NULL || line[0] == '>') {
			if (seq != NULL) {
				for (i=0; i<strlen(seq); i++) {
					seq[i] = (char)toupper((int)seq[i]);
					
					if (removex) {
						/* remove all proteins with occurences of non-natural amino acids:
				    		A CDEFGHI KLMN PQRST VW Y
				    		 B       J    O     U  X Z
				    		 66      74   79    85 8890
						*/
						if (90==seq[i] || 74==seq[i] || 79==seq[i] || 85==seq[i] || 88==seq[i] || 66==seq[i]) {
							free(seq);
							seq = NULL;
							break;
						}
					}
				}
				if (seq != NULL) {
					N++;
					*psequences = realloc(*psequences, N * sizeof(Sequence));
					strcpy((*psequences + N - 1)->description, desc);
					(*psequences + N - 1)->sequence = seq;
					seq = NULL;
				}
			}
			if (line == NULL) {
				break;
			}
			strcpy(desc, line + 1);
			continue;
		}
		
		if (strlen(line) == 0) {
			continue;
		}
						
		if (seq != NULL) {
			seq = realloc(seq, strlen(seq) + strlen(line) + 1);
			strcpy(seq + strlen(seq), line);
		} else {
			seq = calloc(strlen(line) + 1, sizeof(char));
			strcpy(seq, line);
		}
	}
	fclose(fd);
	return N;
}

void print_PSSM(double M[50][26], double threshold) {
	int i=0, j=0, a=0;
	printf("PSSM");
	
	for (j=0; j<20; j++) {
		printf(" %c    ", amino_acids[j]);
	}
	printf("\n");
	for (i=0; i<50; i++) {
		printf("%-2d ", i);
		for (j=0; j<20; j++) {
			a = amino_acids[j] - 'A';
			if (M[i][a] > threshold) {
				printf("% 1.2f ", M[i][a]);
			} else {
				printf("  .   ");
			}
		}
		printf("\n");
	}
	printf("\n");
}


// returns number of loaded matrices
int load_Matrices(const char *filename, Matrix **pmatrices) {
	//char segment[52] = "";
	char *buf = calloc(256, sizeof(char));
	int m=0, j=0;
	int pos=0, profile_format;
	int c[20];
	double b[20];
	int matrix_number = 0;
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
				printf("Format error. Wrong profile format version! Expecting version 1 or 2.\n");
				return 0;
			}
		}
		if (strstr(buf, "BEGIN")) {
			// no conditions
		}
		if (strstr(buf, "END")) {
			//printf("Sizeof Matrix = %d\n", m*sizeof(Matrix));
			m++;
			*pmatrices = realloc(*pmatrices, m*sizeof(Matrix));
			memcpy(*pmatrices + m - 1, &M, sizeof(Matrix));
		}
		if (profile_format == 1) {
			if (1 == sscanf(buf, "MATRIX K=%d", &M.K)) {
				M.id = matrix_number++;
			}
			if (21 == sscanf(buf,
			    "%2d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", 
			     &pos,&c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],&c[8],&c[9],&c[10],&c[11],&c[12],&c[13],&c[14],&c[15],&c[16],&c[17],&c[18],&c[19])) {
		    	
				if (pos<0 || pos>49) {
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
				// no conditions
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
	//printf("Finished parsing file with matrices\n");
	fclose(fd);
	free(buf);
	return m;
}


// loads fasta sequences and creates trivial matrices for them
// returns number of loaded matrices, which are actually consensus sequences
int load_MatricesFromSequence(const char *filename, Matrix **pmatrices, double composition[20]) {
	const int profile_format_version = 0;
	char initial_segment[51];
	int i=0, j=0, k=0, l=0, q=0;
	int m = 0;
	int N_sequences = 0;
	char *seq;
	Matrix M;
	Sequence *sequences = NULL, **psequences = &sequences;

	N_sequences = load_fasta_sequences(filename, psequences, false);

	if (*pmatrices != NULL) {
		printf("Input error. load_Matrices expects an unallocated pointer to matrices (NULL)\n");
		return 0;
	}

	for (k = 0; k < N_sequences; k++) {
		seq = (sequences + k)->sequence;

		for (l=0; l + 50 <= strlen(seq); l+= 10) { // + delta = 10

				initial_segment[50] = '\0';
				strncpy(initial_segment, seq + l, 50);
			
				bool flag1 = false;
    				for (i=50;i--;) {
					q = initial_segment[i];
					if (90==q || 74==q || 79==q || 85==q || 88==q || 66==q) {
						flag1 = true;
						break;
					}
				}
				if (flag1) continue;
			
				// new initial segment obtained here!
				for (i=50;i--;) for (j=26;j--;) {
					M.freq[i][j] = 0;
				}
				
				for (i=0; i<50;i++) {
					q = initial_segment[i] - 'A';
				    	M.freq[i][q] = 1; // initial number of amino acids
				}
						
				M.id = m;
				M.K = 1;
				//M.L = 50;

				*pmatrices = realloc(*pmatrices, (m+1)*sizeof(Matrix));
				memcpy(*pmatrices + (m+1) - 1, &M, sizeof(Matrix));
				
				m++;
		}
	}

	free(sequences);
	return m;
}




// returns number of loaded matrices
// matrix containing gaps, and counts
int load_VariableMatricesCount(const char *filename, Matrix **pmatrices, double composition[20]) {
	const int profile_format_version = 3;
	//char segment[52] = "";
	char *buf = calloc(256, sizeof(char));
	int m=0, j=0;
	int pos=0, tmp;
	int L;
	double b[21];
	FILE *fd = NULL;
	Matrix M;
	bool skip_profile = true;

	if (*pmatrices != NULL) {
		printf("Input error. load_Matrices expects an unallocated pointer to matrices (NULL)\n");
		return 0;
	}

	fd = fopen(filename, "r");
	assert(fd != NULL);

	while (fgets(buf, 256, fd)) {
		if (1 == sscanf(buf, "PROFILE %d", &tmp)) {
			if (tmp != profile_format_version) {
				printf("Format error. Wrong profile format version! Expecting version %d.\n", profile_format_version);
				return 0;
			}
		}
		if (strstr(buf, "BEGIN")) {
			for (pos=0; pos<50; pos++) {
				for (j=20;j--;) {
					M.freq[pos][amino_acids[j] - 'A'] = composition[j];
				}
			}
			skip_profile = true;
		}
		if (strstr(buf, "END")) {
			//printf("Sizeof Matrix = %d\n", m*sizeof(Matrix));
			if (skip_profile) {
				printf("Matrix header not recognized, ignoring profile\n");
			} else {
				m++;
				*pmatrices = realloc(*pmatrices, m*sizeof(Matrix));
				memcpy(*pmatrices + m - 1, &M, sizeof(Matrix));
			}
		}
		if (2 == sscanf(buf, "MATRIX K=%d L=%d", &M.K, &L)) {
			// no conditions
			M.id = m;
			//printf("Matrix %d %d %d\n", m, M.K, L);
			skip_profile = false;
		}
		if (22 == sscanf(buf,
		    "%4d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		     &pos,&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8],&b[9],&b[10],&b[11],&b[12],&b[13],&b[14],&b[15],&b[16],&b[17],&b[18],&b[19],&b[20])) {

			if (M.K <= 0) {
			    continue;
			}
		     
			if (pos<0 || pos>49) {
				printf("Format error: invalid position %d (borders 0-49)\n", pos);
				continue;
			}
			for (j=20; j--;) {
				if (b[j]<0) {
					printf("Format error: invalid position %d,%d = %lf\n", pos, j, b[j]);
				}
				M.freq[pos][amino_acids[j] - 'A'] = b[j] / (double)M.K;
			}
		}
	}
	//printf("Finished parsing file with matrices\n");
	fclose(fd);
	free(buf);
	return m;
}
// returns number of loaded matrices
int load_VariableMatricesFreq(const char *filename, Matrix **pmatrices, double composition[20]) {
	const int profile_format_version = 4;
	//char segment[52] = "";
	char *buf = calloc(256, sizeof(char));
	int m=0, j=0;
	int pos=0, tmp;
	int L;
	double b[20];
	FILE *fd = NULL;
	Matrix M;

	if (*pmatrices != NULL) {
		printf("Input error. load_Matrices expects an unallocated pointer to matrices (NULL)\n");
		return 0;
	}

	fd = fopen(filename, "r");
	assert(fd != NULL);

	while (fgets(buf, 256, fd)) {
		if (1 == sscanf(buf, "PROFILE %d", &tmp)) {
			if (tmp != profile_format_version) {
				printf("Format error. Wrong profile format version! Expecting version %d.\n", profile_format_version);
				return 0;
			}
		}
		if (strstr(buf, "BEGIN")) {
			for (pos=0; pos<50; pos++) {
				for (j=20;j--;) {
					M.freq[pos][amino_acids[j] - 'A'] = composition[j];
				}
			}
		}
		if (strstr(buf, "END")) {
			//printf("Sizeof Matrix = %d\n", m*sizeof(Matrix));
			m++;
			*pmatrices = realloc(*pmatrices, m*sizeof(Matrix));
			memcpy(*pmatrices + m - 1, &M, sizeof(Matrix));
		}
		if (3 == sscanf(buf, "MATRIX ID=%d K=%d L=%d", &M.id, &M.K, &L)) {
			// no conditions
		}
		if (21 == sscanf(buf,
		    "%4d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		     &pos,&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8],&b[9],&b[10],&b[11],&b[12],&b[13],&b[14],&b[15],&b[16],&b[17],&b[18],&b[19])) {
		     
			if (pos<0 || pos>49) {
				printf("Format error: invalid position %d (borders 0-49)\n", pos);
				continue;
			}
			for (j=20; j--;) {
				if (b[j]<0 || b[j]>1) {
					printf("Format error: invalid position %d,%d = %lf\n", pos, j, b[j]);
				}
				M.freq[pos][amino_acids[j] - 'A'] = b[j];
			}
		}
	}
	//printf("Finished parsing file with matrices\n");
	fclose(fd);
	free(buf);
	return m;
}


int main (int argc, char *argv[]) {
	//in R: sample(0:49)	
	int random_index[50] = 
	    {46,20,41,29,21,26,37,13,17,36,47,8,25,43,9,48,39,24,1,5,34,45,32,18,14,40,38,28,7,6,11,27,0,44,33,19,16,35,23,4,42,49,2,15,30,12,3,10,22,31};

	int i=0, j=0, k=0, l=0, m=0, q=0, f=0;
	double composition[20]; for (j=20;j--;) composition[j] = 0.0;
	unsigned long aa_count[26]; for (j=26;j--;) aa_count[j] = 0;

	double composition_control[20]; for (j=20;j--;) composition_control[j] = 0.0;
	unsigned long aa_count_control[26]; for (j=26;j--;) aa_count_control[j] = 0;

	double PSSM[50][26];	for (i=50;i--;) for (j=26;j--;) PSSM[i][j] = 0.0;
	double PSSM_noinf[50][26];	for (i=50;i--;) for (j=26;j--;) PSSM_noinf[i][j] = 0.0;
	double Information[50];	for (i=50;i--;) Information[i] = 0.0;
	double S = 0.0;
		
	char ch = '\0';
	char *protein = NULL;
	
	Matrix *matrix = NULL, *matrices = NULL, **pmatrices = &matrices;
	int matrices_count = 0;
	
	Sequence *subject = NULL, **psubject = &subject;
	Sequence *control = NULL, **pcontrol = &control;
	int N_subject = 0;  	// # subject or haystack
	int N_control = 0;  	// # subject or haystack
	int N = 0; 		// # unmasked fragments in all subject sequences
	int N_masked = 0;    	// # masked fragments in all subject sequences
	
	FILE *fd = NULL;
	FILE *fd_aln = NULL;
	FILE *fd_scores = NULL;

	double E = 1.0;
	double p = 0.0;
	int output_distributions = 0;
	int number_of_randomizations = 0;
	int a = 1;
	int file_format = 2;
		
	// process command line parameters
	for (a=1; a<argc-1; a++) {
		if (argv[a] && !strcmp(argv[a], "-f")) {
			if (argv[a+1] && strlen(argv[a+1])) {
			    file_format = atoi(argv[a+1]);
			    printf("OPTION: Format = %d\n", file_format);
			}
		}
		
		//if (argv[a] && !strcmp(argv[a], "-r")) {
		//	if (argv[a+1] && strlen(argv[a+1])) {
		//	    number_of_randomizations = atoi(argv[a+1]);
		//	    printf("OPTION: Randomizations = %d\n", number_of_randomizations);
		//	}
		//}
		if (argv[a] && !strcmp(argv[a], "-E")) {
			if (argv[a+1] && strlen(argv[a+1])) {
			    E = atof(argv[a+1]);
			    printf("OPTION: E-value = %f\n", E);
			}
		}
		if (argv[a] && !strcmp(argv[a], "-d")) {
			output_distributions = 1;
			printf("OPTION: output distributions = TRUE\n");
		}
	}
	
	N_subject = load_fasta_sequences(subject_filename, psubject, false);
	N_control = load_fasta_sequences(control_filename, pcontrol, false);

	unsigned long size_of_control = 0;
	unsigned long size_of_control_aminoacid = 0;
	bool flag1 = false;
	for (k=0;k<N_control;k++) {
		for (i=0;(ch = (control + k)->sequence[i++]);) {
			aa_count_control[ch-'A']++;
			size_of_control++;
		}
		for (i=0; i+50 <= strlen((control + k)->sequence); i++) {
			// is masked?
			for (j=50;j--;) {
				q = (control + k)->sequence[i+j];
				if (90==q || 74==q || 79==q || 85==q || 88==q || 66==q) {
					flag1 = true;
					break;
				}
			}
			if (flag1) {
			        continue;
			}
			N++;
		}
	}

	unsigned long size_of_subject = 0;
	unsigned long size_of_subject_aminoacid = 0;
	bool flag2 = false;
	for (k=0;k<N_subject;k++) {
		for (i=0;(ch = (subject + k)->sequence[i++]);) {
			aa_count[ch-'A']++;
			size_of_subject++;
		}
		for (i=0; i+50 <= strlen((subject + k)->sequence); i++) {
			// is masked?
			for (j=50;j--;) {
				q = (subject + k)->sequence[i+j];
				if (90==q || 74==q || 79==q || 85==q || 88==q || 66==q) {
					flag2 = true;
					break;
				}
			}
			if (flag2) {
				N_masked++;
			        continue;
			}
			N++;
		}
	}
	
	for (j=20; j--;) {
		q = amino_acids[j]-'A';
		size_of_subject_aminoacid += aa_count[q];
	}
	for (j=20; j--;) {
		q = amino_acids[j]-'A';
		composition[j] = (double)aa_count[q] / size_of_subject_aminoacid;
	}
	for (j=20; j--;) {
		#ifndef NDEBUG
		printf("%c: %f\n", amino_acids[j], composition[j]);
		#endif
	}
		
		
	for (j=20; j--;) {
		q = amino_acids[j]-'A';
		size_of_control_aminoacid += aa_count_control[q];
	}
	for (j=20; j--;) {
		q = amino_acids[j]-'A';
		composition_control[j] = (double)aa_count_control[q] / size_of_control_aminoacid;
	}
	for (j=20; j--;) {
		#ifndef NDEBUG
		printf("Control %c: %f\n", amino_acids[j], composition_control[j]);
		#endif
	}

	if (file_format == 2 || file_format == 1) {
		matrices_count = load_Matrices(extracted_matrix_filename, pmatrices);
	}
	if (file_format == 3) {
		matrices_count = load_VariableMatricesCount(extracted_matrix_filename, pmatrices, composition);
	}	
	if (file_format == 4) {
		matrices_count = load_VariableMatricesFreq(extracted_matrix_filename, pmatrices, composition);
	}		
	if (file_format == 0) {
		matrices_count = load_MatricesFromSequence(extracted_matrix_filename, pmatrices, composition);
	}		
	
	p = E / N;


	fd = fopen(results_filename, "w");
	assert(fd != NULL);

	#ifndef NDEBUG
	printf("searchPSSM: Matrices: %d, Subject sequences: %d, Masked segments: %d, unmasked: %d. p=%E E=%f\n"
	       "---------------------------------------------------------------\n",
	       matrices_count, N_subject, N_masked, N, p, E);
	#endif

	//assert(E > 0.0);	
	if (E < 1.0) E=1.0;

	char scores_filename[50];
	int last_informative_position = 0;
	double information_sum = 0.0;
				
	for (m=0;m<matrices_count;m++) {
		matrix = matrices + m;
		
			for (i=50;i--;) for (j=26;j--;) {
				PSSM[i][j] = 0.0;
				PSSM_noinf[i][j] = 0.0;
			}
			for (i=50;i--;) {
				Information[i] = 0.0;
			}
			
			information_sum = 0.0;
			
			double sum_aa = 0;
			double fragment_freq = 0.0;
			double pc = 0.0;
			double H_sum = 0.0; // attention: it is in nats, not in bits
			
			print_PSSM(matrix->freq, 0.0);
			
			last_informative_position = 0;
			for (i=0;i<50;i++) {
				H_sum = 0.0;
				for (j=20;j--;) {
					q = amino_acids[j] - 'A';
					if (matrix->freq[i][q] > 0.0) {
						H_sum += matrix->freq[i][q] * (log(matrix->freq[i][q])/log(2));
					}
				}
				Information[i] = log(20)/log(2) + H_sum;
				if (Information[i] < 1.0) {
				    Information[i] = 0.0;
				} else {
				    last_informative_position = i;
				}
				
				information_sum += Information[i];
				
				printf("%d %f %d\n", i, Information[i], last_informative_position);
			}
			// hack!
			last_informative_position = 49;
			
			for (j=20;j--;) {
				q = amino_acids[j] - 'A';
			    	
			    	sum_aa = 0;
			    	for (i=50;i--;) {
				    sum_aa += matrix->freq[i][q] * matrix->K;
				}
				fragment_freq = sum_aa/(50.0 * matrix->K);
			    	//compensating by observed frequencies in aligned fragment
				pc = 0.001;
				if (fragment_freq < composition[j]) {
					pc = composition[j] - fragment_freq;// /composition[j];
				}
				printf("PC[%d][%c] = %f\n", i, q+'A', pc);
				for (i=50;i--;) {
					PSSM[i][q] =  Information[i] * 
							    log( ((matrix->freq[i][q]*matrix->K + omega*pc) / (matrix->K + 20.0*omega*pc)) 
							    / composition[j] );

					PSSM_noinf[i][q] =   log( ((matrix->freq[i][q]*matrix->K + omega*pc) / (matrix->K + 20.0*omega*pc)) 
							    / composition[j] );
			    	}
			}

			List *list = NULL, *item = NULL, *curr = NULL, *prev = NULL;
			int list_length = 0;
			double minS = 0.0;


			print_PSSM(PSSM_noinf, 0);

			////////////////////////// SHUFFLED PROTEOME. NO INFORMATION ///////////////////////////////////////////////
			
			list = NULL;
			item = NULL;
			curr = NULL;
			prev = NULL;
			
			list_length = 0;
			minS = 0.0;

			if (output_distributions) {
				snprintf(scores_filename, 50, "scores_%d_shuffled_proteome_noinf.tab", matrix->id); 
				fd_scores = fopen(scores_filename, "w");
	    		}

			for (f=0; f<N_control; f++) {
				protein = (control + f)->sequence;
				
				int protein_length = strlen(protein);
				for (l=0; l+50 <= protein_length; l++) {
					S = 0.0;
					for (i=0;i<50;i++) {
						q = protein[l+i];
						if (90==q || 74==q || 79==q || 85==q || 88==q || 66==q) {
							S = -10000000; // -inf
							break;
						} else {
							S +=  PSSM_noinf[i][q-'A'];
						}
					}
					// S /= 50;
					
					if (output_distributions && S > -10000000) {
						fprintf(fd_scores, "%f\n", S);
					}
					
					if (list_length < E || S > minS) {
						//printf("S=%f, minS=%f, list=%d, E=%d\n", S, minS, list_length, (int)FP);
						curr = list;
						prev = NULL;
						while (curr != NULL && curr->score < S) {
							prev = curr;
							curr = curr->next;
						}
						item = malloc(sizeof(List));
						item->score = S;
						item->next = curr;
						if (list == item->next) {
							list = item;
						}
						if (prev != NULL) {
							prev->next = item;
						}						
						list_length++;
												
						if (list_length > E) {
							//truncate one from list head
							prev = list;
							list = list->next;
							free(prev);
							list_length--;
						}
					}
					minS = list->score;
				}			
			}
			while (list != NULL) {
				prev = list;
				list = list->next;
				free(prev);
				list_length--;
			}
			assert(list_length == 0);
			
			if (output_distributions) {
				fclose(fd_scores);
			}
			//minS += 5;
						
			
			// min score
			#ifndef NDEBUG
			printf("Matrix %d: Shuffled proteome. No information |E(s>%f)| <= %f \n", matrix->id, minS, E);
			#endif


			///////////////////////////////// SHUFFLED PROFILE. NO INFORMATION /////////////////////////////////////

			list = NULL;
			item = NULL;
			curr = NULL;
			prev = NULL;
			
			list_length = 0;
			minS = 0.0;

			if (output_distributions) {
				snprintf(scores_filename, 50, "scores_%d_shuffled_profile_noinf.tab", matrix->id); 
				fd_scores = fopen(scores_filename, "w");
	    		}

			for (f=0; f<N_subject; f++) {
				protein = (subject + f)->sequence;
				
				int protein_length = strlen(protein);
				for (l=0; l+50 <= protein_length; l++) {
					S = 0.0;
					for (i=0;i<50;i++) {
						// compare 50 residues in one fragment from proteome with reshuffled matrix position
						q = protein[l+i];
						if (90==q || 74==q || 79==q || 85==q || 88==q || 66==q) {
							S = -10000000; // -inf
							break;
						} else {
							S +=  PSSM_noinf[random_index[i]][q-'A'];
						}
					}
					//S /= 50;

					
					if (output_distributions && S > -10000000) {
						fprintf(fd_scores, "%f\n", S);
					}
					
					if (list_length < E || S > minS) {
						//printf("S=%f, minS=%f, list=%d, E=%d\n", S, minS, list_length, (int)FP);
						curr = list;
						prev = NULL;
						while (curr != NULL && curr->score < S) {
							prev = curr;
							curr = curr->next;
						}
						item = malloc(sizeof(List));
						item->score = S;
						item->next = curr;
						if (list == item->next) {
							list = item;
						}
						if (prev != NULL) {
							prev->next = item;
						}						
						list_length++;
												
						if (list_length > E) {
							//truncate one from list head
							prev = list;
							list = list->next;
							free(prev);
							list_length--;
						}
					}
					minS = list->score;
				}			
			}
			while (list != NULL) {
				prev = list;
				list = list->next;
				free(prev);
				list_length--;
			}
			assert(list_length == 0);
			
			if (output_distributions) {
				fclose(fd_scores);
			}
			
			//minS += 5;
						
			
			// min score
			#ifndef NDEBUG
			printf("Matrix %d: Shuffled profile. no information: |E(s>%f)| <= %f \n", matrix->id, minS, E);
			#endif




			print_PSSM(PSSM, 0);


			////////////////////////// SHUFFLED PROTEOME ///////////////////////////////////////////////
			
			list = NULL;
			item = NULL;
			curr = NULL;
			prev = NULL;
			
			list_length = 0;
			minS = 0.0;

			if (output_distributions) {
				snprintf(scores_filename, 50, "scores_%d_shuffled_proteome.tab", matrix->id); 
				fd_scores = fopen(scores_filename, "w");
	    		}

			for (f=0; f<N_control; f++) {
				protein = (control + f)->sequence;
				
				int protein_length = strlen(protein);
				for (l=0; l+50 <= protein_length; l++) {
					S = 0.0;
					for (i=0;i<50;i++) {
						q = protein[l+i];
						if (90==q || 74==q || 79==q || 85==q || 88==q || 66==q) {
							S = -10000000; // -inf
							break;
						} else {
							S +=  PSSM[i][q-'A'];
						}
					}
					
					S /= information_sum;
					//S /= 50;
					
					if (output_distributions && S > -10000000) {
						fprintf(fd_scores, "%f\n", S);
					}
					
					if (list_length < E || S > minS) {
						//printf("S=%f, minS=%f, list=%d, E=%d\n", S, minS, list_length, (int)FP);
						curr = list;
						prev = NULL;
						while (curr != NULL && curr->score < S) {
							prev = curr;
							curr = curr->next;
						}
						item = malloc(sizeof(List));
						item->score = S;
						item->next = curr;
						if (list == item->next) {
							list = item;
						}
						if (prev != NULL) {
							prev->next = item;
						}						
						list_length++;
												
						if (list_length > E) {
							//truncate one from list head
							prev = list;
							list = list->next;
							free(prev);
							list_length--;
						}
					}
					minS = list->score;
				}			
			}
			while (list != NULL) {
				prev = list;
				list = list->next;
				free(prev);
				list_length--;
			}
			assert(list_length == 0);
			
			if (output_distributions) {
				fclose(fd_scores);
			}
			//minS += 5;
						
			
			// min score
			#ifndef NDEBUG
			printf("Matrix %d: CONTROL BY SHUFFLING PROTEOME GIVES |E(s>%f)| <= %f \n", matrix->id, minS, E);
			#endif


			///////////////////////////////// SHUFFLED PROFILE /////////////////////////////////////

			list = NULL;
			item = NULL;
			curr = NULL;
			prev = NULL;
			
			list_length = 0;
			minS = 0.0;

			if (output_distributions) {
				snprintf(scores_filename, 50, "scores_%d_shuffled_profile.tab", matrix->id); 
				fd_scores = fopen(scores_filename, "w");
	    		}

			for (f=0; f<N_subject; f++) {
				protein = (subject + f)->sequence;
				
				int protein_length = strlen(protein);
				for (l=0; l+50 <= protein_length; l++) {
					S = 0.0;
					for (i=0;i<50;i++) {
						// compare 50 residues in one fragment from proteome with reshuffled matrix position
						q = protein[l+i];
						if (90==q || 74==q || 79==q || 85==q || 88==q || 66==q) {
							S = -10000000; // -inf
							break;
						} else {
							S +=  PSSM[random_index[i]][q-'A'];
						}
					}
					S /= information_sum;
					//S /= 50;
					
					if (output_distributions && S > -10000000) {
						fprintf(fd_scores, "%f\n", S);
					}
					
					if (list_length < E || S > minS) {
						//printf("S=%f, minS=%f, list=%d, E=%d\n", S, minS, list_length, (int)FP);
						curr = list;
						prev = NULL;
						while (curr != NULL && curr->score < S) {
							prev = curr;
							curr = curr->next;
						}
						item = malloc(sizeof(List));
						item->score = S;
						item->next = curr;
						if (list == item->next) {
							list = item;
						}
						if (prev != NULL) {
							prev->next = item;
						}						
						list_length++;
												
						if (list_length > E) {
							//truncate one from list head
							prev = list;
							list = list->next;
							free(prev);
							list_length--;
						}
					}
					minS = list->score;
				}			
			}
			while (list != NULL) {
				prev = list;
				list = list->next;
				free(prev);
				list_length--;
			}
			assert(list_length == 0);
			
			if (output_distributions) {
				fclose(fd_scores);
			}
			
			//minS += 5;
						
			
			// min score
			#ifndef NDEBUG
			printf("Matrix %d: Shuffled profile: |E(s>%f)| <= %f \n", matrix->id, minS, E);
			#endif



			//////////////////////////// REAL SCORES /////////////////////////////////////////////////


			
			char match_header[256];
			char match_sequence[51];
			char alignment_filename[50];
			snprintf(alignment_filename, 50, "search_aln_%d.fasta", matrix->id); 
			
			fd_aln = fopen(alignment_filename, "w");
			assert(fd_aln != NULL);
			
			if (output_distributions) {
				snprintf(scores_filename, 50, "scores_%d.tab", matrix->id); 
				fd_scores = fopen(scores_filename, "w");
			}
			
			// Use S(p) as threshold on natural matrix matches:
			for (f=0; f<N_subject; f++) {
				protein = (subject + f)->sequence;
				int protein_length = strlen(protein);
				for (l=0; l+50 <= protein_length; l++) {
					S = 0.0;
					for (i=0;i<50;i++) {
						q = protein[l+i];
						if (90==q || 74==q || 79==q || 85==q || 88==q || 66==q) {
							S = -10000000; // -inf
							break;
						} else {
							S += PSSM[i][q-'A'];
						}
					}
					S /= information_sum;
					//S /= 50;
					
					if (output_distributions && S > -10000000) {
						fprintf(fd_scores, "%f\n", S);
					}

					if (S > minS) {
						#ifndef NDEBUG
						printf("match: [%d] [%2.2f] [%-90.90s]\n", matrix->id, S, (subject + f)->description);
						#endif
						
						fprintf(fd, "%d\t%2.2f\t%s\n", matrix->id, S, (subject + f)->description);
						snprintf(match_header, 120, (subject + f)->description);
						fprintf(fd_aln, ">%s\n", match_header);
						snprintf(match_sequence, last_informative_position+2, protein + l);
						fprintf(fd_aln, "%s\n", match_sequence);
					}		
				}
			}
			
			fclose(fd_aln);
			if (output_distributions) {
				fclose(fd_scores);
			}
	}

	fclose(fd);
	return(0);
}

