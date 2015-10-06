#pragma once

#include <cmath>
#include <kore/typedefs.h>
#include <kore/array/array.h>
#include <kore/array/array_linalg.h>
#include <Eigen/Core>

namespace pome
{

typedef double Real;
typedef std::complex<double> Complex;
typedef kore::int64_t Integer;

typedef Eigen::VectorXcd ComplexVector;
typedef Eigen::VectorXd RealVector;
typedef Eigen::MatrixXcd ComplexMatrix;
typedef Eigen::MatrixXd RealMatrix;


}
