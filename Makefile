
ELSI_PATH = /home/vama/soft/elsi/elsi-interface
ELSI_LIB = ${ELSI_PATH}/lib
ELSI_INC = ${ELSI_PATH}/include

LIBS = -L${ELSI_LIB}  -lelsi -lelpa -L${ELSI_LIB} -lcheck_singularity -L${ELSI_LIB} -lOMM -lMatrixSwitch -lpspblas -ltomato  

INC = -I${ELSI_INC} -I${ELSI_PATH}/src/ELSI


FFLAGS = -g  -O0 -fpp
LDFLAGS = -g -mkl=cluster -O0
MPIF90 = mpiifort

all: test_standard_ev_real.x test_generalized_ev_real.x

test_standard_ev_real.x: test_standard_ev_real.o
	${MPIF90} ${LDFLAGS}  -o test_standard_ev_real.x test_standard_ev_real.o  ${LIBS}

test_standard_ev_real.o: test_standard_ev_real.f90
	${MPIF90} ${FFLAGS}  -c ${INC}  test_standard_ev_real.f90

test_generalized_ev_real.x: test_generalized_ev_real.o
	${MPIF90} ${LDFLAGS}  -o test_generalized_ev_real.x test_generalized_ev_real.o  ${LIBS}

test_generalized_ev_real.o: test_generalized_ev_real.f90
	${MPIF90} ${FFLAGS}  -c ${INC}  test_generalized_ev_real.f90


clean:
	rm test*.o test*.x 
