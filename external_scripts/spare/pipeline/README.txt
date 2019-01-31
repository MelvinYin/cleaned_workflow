
Profile derivation protocol
---------------------------

1. Input data - proteome sequence in fasta format

Non-redundant sequences required, 50% redundancy or less, for instance UniRef50 from Uniprot.
CD-HIT tool can be used to reduce redundancy. For example:
cdhit -c 0.5 -n 3 -i all.fasta -o all_NR50.fasta

It is recommended to mask segments with low complexity with SEG, for example:
./seg uniref50.fasta -x -n >uniref50.masked.fasta

2. Build software

mpich2 is the recommended MPI library. It works with OpenMPI, but I had some problems with clustering in certain configurations of OpenMPI.
Tweak Makefile if needed. Then run:
make clean
make all
This version compiles on Linux. With minor modifications it can also compile on Mac.

There are three executables:
converge - converges profiles
clustering - hierarchical clustering and merging of converged profiles
search - search profiles against sequence database


3. Converge profiles

converge has the following parameters:
 -E 1 Evalue = 1. Could be 1 or more (integer)
 -r 1 The number of randomizations
 -f 4 Matrix format: 1,2,3,4
 -p <proteome filename>
 -i <initial segments or initial matrix filename>
 -L 30 Length or profiles. Should be 50 or 30.
 -B use BLOSUM62 matrix
 -K 1 Do not converge profiles with less than K hits
 -o <output file name>

The error checking in parameter processing is minimal, in case of doubt check the source code.
mpirun -np 2 specifies 2 processors. Fore more parameters man mpirun.

Example:
mpirun -np 2 ./converge -B -E 1 -r 1 -f 4 -p proteome.fasta -i initial.fasta

As a result, there will be *.matrix files with converged profiles

4. Clustering 

    It is only needed to run clustering if there are too many profiles and many of them are potentially very similar.
    Otherwise, this step can be skipped.
    
    clustering has no parameters. It looks for output.matrix 
    and writes files to clustering.output directory, which should exist beforehand.
    Clustering accepts matrix files in format 1.

5. Profile logo
    profile_logo.py is a python script that requires Weblogo 3 to be installed
    
    typically it is called like this:
    python profile_logo.py combined.matrix all 2.0
    
    which means that all profiles contained in combined.matrix file will be processed and Shannon information 2.0 bits will be used as a threshold in signatures.
    files will be generated in profile.logo directory
    A file with profile signatures will also be generated using 2.0 bits as a threshold. Positions with information less than 2.0 will not be shown as letters.
    
6. Searching the sequence database
    Search can also be run as a parallel program with MPI:
    
    search's two important parameters are 
    -E 1  E-value
    -f 4  matrix file format (4 recommended)
    
    This code will run search on 8 CPU using E-value <= 1 and format 4.
    mpirun -np 8 ./search -E 1 -f 4
    
    Search expects input files with fixed filenames:
    extracted.matrix  - is the file containing all matrices that will be used as queries in search
    composition.csv - Amino acid composition of the searched proteome in csv format (example file provided)
    subject.fasta - is the sequence database in fasta format that will be searched as subject
    
    The program will produce results in tab-delimited format saved as search_matches.tab
    The fields are as follows: profile, score, evalue, position, fasta header of the sequence hit




APPENDIX. PROFILE file formats
------------------------------

PROTOTYPE 1
BEGIN
SEGMENT HVHPKDLISGEITPIERRGYPAPIVNHNLRQKQFKALYNQLKAAIAEPEA
MATRIX K=637 N=19305087 P=0.00000005 S=0.708697 W=1.000000
50     A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y
 0    57    11     2     6     9     4     4    91    13    20    11     1     3     5     4    54   166   161     1    14
 1    37     9    15    17     1     9    13    92    36    38     6    28    23     7    36    52    35   149     7    27
 2    43     4    65   129     2    20    10    13    22    12     3    25   111    11    18    78    46    16     0     9
.
.
.
END

--------------------
PROTOTYPE 2
BEGIN
MATRIX ID=0 K=637
50        A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
 0 0.089482 0.017268 0.003140 0.009419 0.014129 0.006279 0.006279 0.142857 0.020408 0.031397 0.017268 0.001570 0.004710 0.007849 0.006279 0.084772 0.260597 0.252747 0.001570 0.021978
 1 0.058085 0.014129 0.023548 0.026688 0.001570 0.014129 0.020408 0.144427 0.056515 0.059655 0.009419 0.043956 0.036107 0.010989 0.056515 0.081633 0.054945 0.233909 0.010989 0.042386
.
.
.
END

--------------------
PROFILE 3
BEGIN
MATRIX K=45 L=2
   2    A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y    .
   0    1    0    0    2    0    1    1    1    1    0    1    2    0    0    0    0    0    1    0    1   33
   1    1    0    0    0    2    0    0    0    1    1    0    1    2    0    1    0    2    1    0    0   33
.
.
.
END

--------------------
PROFILE 4
BEGIN
MATRIX ID=0 K=30 L=30
30        A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
 0 0.089482 0.017268 0.003140 0.009419 0.014129 0.006279 0.006279 0.142857 0.020408 0.031397 0.017268 0.001570 0.004710 0.007849 0.006279 0.084772 0.260597 0.252747 0.001570 0.021978
 1 0.058085 0.014129 0.023548 0.026688 0.001570 0.014129 0.020408 0.144427 0.056515 0.059655 0.009419 0.043956 0.036107 0.010989 0.056515 0.081633 0.054945 0.233909 0.010989 0.042386
.
.
.
.
END

