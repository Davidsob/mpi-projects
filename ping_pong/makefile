EXECS=mpi_ping_pong
MPICXX?=mpicxx
TARGET=PingPong
all: ${TARGET}

${TARGET}: ${EXECS}.cpp
	${MPICXX} -o ${TARGET} ${EXECS}.cpp


clean:
	rm $(TARGET)
