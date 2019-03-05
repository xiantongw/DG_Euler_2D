#ifndef CONSTRUCTCURVEMESH_H
#define CONSTRUCTCURVEMESH_H

#include <iostream>
#include <algorithm>
#include <boost/numeric/ublas/io.hpp>

#include "../include/TriMesh.h"
#include "../include/utils.h"
#include "../include/lagrange.h"
#include "../include/geometry.h"


TriMesh ConstructCurveMesh(TriMesh& mesh, TriMesh& curved_mesh, double (*pBumpFunction)(double), string boundary_name, int p);


#endif