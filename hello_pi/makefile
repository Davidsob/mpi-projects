EXECS=mpi_pi
MPICXX?=mpicxx
TARGET=Pi
all: ${TARGET}

${TARGET}: ${EXECS}.cpp
	${MPICXX} -o ${TARGET} ${EXECS}.cpp


clean:
	rm $(TARGET)
