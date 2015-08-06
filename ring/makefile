EXECS=mpi_ring
MPICXX?=mpicxx
TARGET=HotPotato
all: ${TARGET}

${TARGET}: ${EXECS}.cpp
	${MPICXX} -o ${TARGET} ${EXECS}.cpp


clean:
	rm $(TARGET)
