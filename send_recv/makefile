EXECS=mpi_send_recv
MPICXX?=mpicxx
TARGET=Send_And_Receive
all: ${TARGET}

${TARGET}: ${EXECS}.cpp
	${MPICXX} -o ${TARGET} ${EXECS}.cpp


clean:
	rm $(TARGET)
