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
		${BUILD_DIR}/ConstructCurveMesh.o ${BUILD_DIR}/GetQuadraturePointsWeight2D.o \
		${BUILD_DIR}/GetQuadraturePointsWeight1D.o ${BUILD_DIR}/solver.o \
		${BUILD_DIR}/euler.o ${BUILD_DIR}/Collective.o

OBJECTS_POSTPROC = ${BUILD_DIR}/TriMesh.o ${BUILD_DIR}/utils.o ${BUILD_DIR}/geometry.o\
		${BUILD_DIR}/lagrange.o ${BUILD_DIR}/InvertMatrix.o \
		${BUILD_DIR}/ConstructCurveMesh.o ${BUILD_DIR}/GetQuadraturePointsWeight2D.o \
		${BUILD_DIR}/GetQuadraturePointsWeight1D.o ${BUILD_DIR}/solver.o \
		${BUILD_DIR}/euler.o ${BUILD_DIR}/Collective.o

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

${BUILD_DIR}/euler.o: ${SRC_DIR}/euler.cpp ${INCLUDE_DIR}/euler.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/euler.cpp -o ${BUILD_DIR}/euler.o

${BUILD_DIR}/InvertMatrix.o: ${SRC_DIR}/InvertMatrix.cpp ${INCLUDE_DIR}/InvertMatrix.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/InvertMatrix.cpp -o ${BUILD_DIR}/InvertMatrix.o

${BUILD_DIR}/ConstructCurveMesh.o: ${SRC_DIR}/ConstructCurveMesh.cpp ${INCLUDE_DIR}/ConstructCurveMesh.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/ConstructCurveMesh.cpp -o ${BUILD_DIR}/ConstructCurveMesh.o

${BUILD_DIR}/GetQuadraturePointsWeight2D.o: ${SRC_DIR}/GetQuadraturePointsWeight2D.cpp ${INCLUDE_DIR}/GetQuadraturePointsWeight2D.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/GetQuadraturePointsWeight2D.cpp -o ${BUILD_DIR}/GetQuadraturePointsWeight2D.o

${BUILD_DIR}/GetQuadraturePointsWeight1D.o: ${SRC_DIR}/GetQuadraturePointsWeight1D.cpp ${INCLUDE_DIR}/GetQuadraturePointsWeight1D.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/GetQuadraturePointsWeight1D.cpp -o ${BUILD_DIR}/GetQuadraturePointsWeight1D.o

${BUILD_DIR}/solver.o: ${SRC_DIR}/solver.cpp ${INCLUDE_DIR}/solver.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/solver.cpp -o ${BUILD_DIR}/solver.o

${BUILD_DIR}/Collective.o: ${SRC_DIR}/Collective.cpp ${INCLUDE_DIR}/Collective.h | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/Collective.cpp -o ${BUILD_DIR}/Collective.o

${BUILD_DIR}/main.o: ${SRC_DIR}/main.cpp | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/main.cpp -o ${BUILD_DIR}/main.o

${BUILD_DIR}:
	mkdir -p build

${BUILD_DIR}/PostProc.o: ${SRC_DIR}/PostProc.cpp | ${BUILD_DIR}
	${CC} ${CPPFLAG} -c ${SRC_DIR}/PostProc.cpp -o ${BUILD_DIR}/PostProc.o

postproc: ${OBJECTS_POSTPROC} ${BUILD_DIR}/PostProc.o
	${CC} ${OBJECTS_POSTPROC} ${BUILD_DIR}/PostProc.o -o postproc.exe

rundir:
	mkdir -p run
	cd run; ln -s ../solver.exe .; ln -s ../postproc.exe .; cp ../PARAM* .; cp -r ../mesh .

clean:
	rm -rf *.exe
	rm -rf ${BUILD_DIR}*

