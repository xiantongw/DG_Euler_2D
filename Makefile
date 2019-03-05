CC = icpc

STDFLAG = -std=c++11
OPTFLAG = -O3
INCLUDEPATH = -I/opt/local/include

CPPFLAG = ${STDFLAG} ${OPTFLAG} ${INCLUDEPATH} -Wno-unknown-pragmas

SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

OBJECTS = ${BUILD_DIR}/TriMesh.o ${BUILD_DIR}/utils.o ${BUILD_DIR}/geometry.o\
			 ${BUILD_DIR}/lagrange.o ${BUILD_DIR}/main.o ${BUILD_DIR}/InvertMatrix.o \
			 ${BUILD_DIR}/ConstructCurveMesh.o

solver : ${OBJECTS}
	${CC} ${OBJECTS} -o solver.exe

${BUILD_DIR}/TriMesh.o: ${SRC_DIR}/TriMesh.cpp ${INCLUDE_DIR}/TriMesh.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/TriMesh.cpp -o ${BUILD_DIR}/TriMesh.o

${BUILD_DIR}/utils.o: ${SRC_DIR}/utils.cpp ${INCLUDE_DIR}/utils.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/utils.cpp -o ${BUILD_DIR}/utils.o

${BUILD_DIR}/lagrange.o: ${SRC_DIR}/lagrange.cpp ${INCLUDE_DIR}/lagrange.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/lagrange.cpp -o ${BUILD_DIR}/lagrange.o

${BUILD_DIR}/geometry.o: ${SRC_DIR}/geometry.cpp ${INCLUDE_DIR}/geometry.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/geometry.cpp -o ${BUILD_DIR}/geometry.o

${BUILD_DIR}/InvertMatrix.o: ${SRC_DIR}/InvertMatrix.cpp ${INCLUDE_DIR}/InvertMatrix.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/InvertMatrix.cpp -o ${BUILD_DIR}/InvertMatrix.o

${BUILD_DIR}/ConstructCurveMesh.o: ${SRC_DIR}/ConstructCurveMesh.cpp ${INCLUDE_DIR}/ConstructCurveMesh.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/ConstructCurveMesh.cpp -o ${BUILD_DIR}/ConstructCurveMesh.o

${BUILD_DIR}/main.o: ${SRC_DIR}/main.cpp | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/main.cpp -o ${BUILD_DIR}/main.o

${BUILD_DIR}:
	mkdir -p build

clean:
	rm -rf solver.exe
	rm -rf ${BUILD_DIR}
