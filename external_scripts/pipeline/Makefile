CC=mpicc
#CFLAGS=-Wall -O2 -pipe -fforce-addr -fomit-frame-pointer -funroll-loops -std=c99 -pedantic
CFLAGS=-Wall -O2 -pipe -funroll-loops -std=c99 -pedantic
CXXFLAGS=${CFLAGS}

all: mpiconverge clustering search

converge: convergePSSM.c
	${CC} -D NDEBUG ${CFLAGS} convergePSSM.c -lm -o converge

mpiconverge: convergePSSM.c
	${CC} -D MPI -D NDEBUG ${CFLAGS} convergePSSM.c -lm -o converge

debug: clean convergePSSM.c
	${CC} -g ${CFLAGS} convergePSSM.c -lm -o converge

mpidebug: clean convergePSSM.c
	${CC} -g ${CFLAGS} -D MPI  convergePSSM.c -lm -o converge

clustering: clean clustering.c
	${CC} -D MPI -g ${CFLAGS} clustering.c -lm -o clustering

clean:
	rm -f converge
	rm -f clustering
	rm -f search

search: searchPSSM.c
	rm -f search
	${CC} -D MPI -D NDEBUG ${CFLAGS} searchPSSM.c -lm -o search

searchdebug: searchPSSM.c
	rm -f search
	${CC} -D MPI -g ${CFLAGS} searchPSSM.c -lm -o search

searchcontrol: searchPSSM_control.c
	rm -f search_control
	${CC} -g ${CFLAGS} searchPSSM_control.c -lm -o search_control

convert_matrix_50_30: convert_matrix_50_30.c
	rm -f convert_matrix_50_30
	cc -D NDEBUG ${CFLAGS} convert_matrix_50_30.c -lm -o convert_matrix_50_30
    