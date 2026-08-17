// CARDfree.hpp
#ifndef CARDFREE_HPP
#define CARDFREE_HPP

#include <pybind11/pybind11.h>
#include <armadillo>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>

namespace py = pybind11;

py::dict CARDfree(const Eigen::Ref<const Eigen::MatrixXd>& XinputIn,
                    const Eigen::Ref<const Eigen::MatrixXd>& UIn,
                    const Eigen::Ref<const Eigen::MatrixXd>& WIn,
                    double phiIn,
                    int max_iterIn,
                    double epsilonIn,
                    const Eigen::Ref<const Eigen::MatrixXd>& initV,
                    const Eigen::Ref<const Eigen::VectorXd>& initb,
                    double initSigma_e2,
                    const Eigen::Ref<const Eigen::VectorXd>& initLambda);

#endif // CARDFREE_HPP
