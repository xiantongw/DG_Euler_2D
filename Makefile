CC = icpc

STDFLAG = -std=c++11
OPTFLAG = -O0 -g

CPPFLAG = ${STDFLAG} ${OPTFLAG}

SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

OBJECTS = ${BUILD_DIR}/TriMesh.o ${BUILD_DIR}/utils.o ${BUILD_DIR}/main.o 

solver : ${OBJECTS}
	${CC} ${OBJECTS} -o solver.exe

${BUILD_DIR}/TriMesh.o: ${SRC_DIR}/TriMesh.cpp ${INCLUDE_DIR}/TriMesh.h
	${CC} ${CPPFLAG} -c ${SRC_DIR}/TriMesh.cpp -o ${BUILD_DIR}/TriMesh.o

${BUILD_DIR}/utils.o: ${SRC_DIR}/utils.cpp ${INCLUDE_DIR}/utils.h
	${CC} ${CPPFLAG} -c ${SRC_DIR}/utils.cpp -o ${BUILD_DIR}/utils.o

${BUILD_DIR}/main.o: ${SRC_DIR}/main.cpp
	${CC} ${CPPFLAG} -c ${SRC_DIR}/main.cpp -o ${BUILD_DIR}/main.o 

clean:
	rm -rf solver.exe
	rm -rf ${OBJECTS}
