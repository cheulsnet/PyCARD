// CARDref.cpp
#include "CARDref.hpp"
#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>

namespace py = pybind11;

//********************************************************************************************
//  Change log
//  1. use armadillo library for arma::
//  2. use pybind11 to work with python
//  3. need to conversion function from Eigen::MatrixXd to arma::mat without copying
//     it's required to return the result in a dictionary format
//********************************************************************************************

template<typename T>
std::vector<std::vector<T>> armaToVector(const arma::Mat<T>& matrix) {
    std::vector<std::vector<T>> result(matrix.n_rows, std::vector<T>(matrix.n_cols));
    for (arma::uword i = 0; i < matrix.n_rows; ++i) {
        for (arma::uword j = 0; j < matrix.n_cols; ++j) {
            result[i][j] = matrix(i, j);
        }
    }
    return result;
}

// Template function for converting Eigen::MatrixXd to arma::mat without copying
template <typename Derived>
arma::Mat<typename Derived::Scalar> eigenToArma(const Eigen::MatrixBase<Derived>& matrix) {
    using Scalar = typename Derived::Scalar;
    return arma::Mat<Scalar>(
        const_cast<Scalar*>(matrix.derived().data()), matrix.rows(), matrix.cols(), false, true
    );
}

//*******************************************************************//
//              spatially informed deconvolution:CARD                        //
//*******************************************************************//
//' SpatialDeconv function based on Conditional Autoregressive model
//' @param XinputIn The input of normalized spatial data
//' @param UIn The input of cell type specific basis matrix B
//' @param WIn The constructed W weight matrix from Gaussian kernel
//' @param phiIn The phi value
//' @param max_iterIn Maximum iterations
//' @param epsilonIn epsilon for convergence 
//' @param initV Initial matrix of cell type compositions V
//' @param initb Initial vector of cell type specific intercept
//' @param initSigma_e2 Initial value of residual variance
//' @param initLambda Initial vector of cell type sepcific scalar. 
//'
//' @return A list
//'

py::dict CARDref(const Eigen::Ref<const Eigen::MatrixXd>& XinputIn,
                 const Eigen::Ref<const Eigen::MatrixXd>& UIn,
                 const Eigen::Ref<const Eigen::MatrixXd>& WIn,
                 double phiIn,
                 int max_iterIn,
                 double epsilonIn,
                 const Eigen::Ref<const Eigen::MatrixXd>& initV,
                 const Eigen::Ref<const Eigen::VectorXd>& initb,
                 double initSigma_e2,
                 const Eigen::Ref<const Eigen::VectorXd>& initLambda) {

    try {
        // read in the data
        arma::mat Xinput = eigenToArma(XinputIn);
        arma::mat U = eigenToArma(UIn);
        arma::mat W = eigenToArma(WIn);
        double phi = phiIn;
        int max_iter = max_iterIn;
        double epsilon = epsilonIn;
        arma::mat V = eigenToArma(initV);
        arma::vec b = eigenToArma(initb);
        double sigma_e2 = initSigma_e2;
        arma::vec lambda = eigenToArma(initLambda);
        
        // initialize some useful items
        int nSample = (int)Xinput.n_cols; // number of spatial sample points
        int mGene = (int)Xinput.n_rows; // number of genes in spatial deconvolution
        int k = (int)U.n_cols; // number of cell type

        arma::mat L = arma::zeros<arma::mat>(nSample, nSample);
        arma::mat D = arma::zeros<arma::mat>(nSample, nSample);
        arma::mat V_old = arma::zeros<arma::mat>(nSample, k);
        arma::mat UtU = arma::zeros<arma::mat>(k, k);
        arma::mat VtV = arma::zeros<arma::mat>(k, k);
        arma::vec colsum_W = arma::zeros<arma::vec>(nSample);
        arma::mat UtX = arma::zeros<arma::mat>(k,nSample);
        arma::mat XtU = arma::zeros<arma::mat>(nSample,k);
        arma::mat UtXV = arma::zeros<arma::mat>(k,k);
        arma::mat temp = arma::zeros<arma::mat>(k,k);
        arma::mat part1 = arma::zeros<arma::mat>(nSample,k);
        arma::mat part2 = arma::zeros<arma::mat>(nSample,k);
        arma::vec updateV_k = arma::zeros<arma::vec>(k);
        arma::vec updateV_den_k = arma::zeros<arma::vec>(k);
        arma::vec vecOne = arma::ones<arma::vec>(nSample);
        arma::vec diag_UtU = arma::zeros<arma::vec>(k);
        bool logicalLogL = false;
        double obj = 0;
        double obj_old = 0;
        double normNMF = 0;
        double logX = 0;
        double logV = 0;
        double alpha = 1.0;
        double beta = nSample / 2.0;
        double logSigmaL2 = 0.0;
        double accu_L = 0.0;
        double trac_xxt = accu(Xinput % Xinput);

        // initialize values
        // constant matrix caculations for increasing speed
        UtX = U.t() * Xinput;
        XtU = UtX.t();
        colsum_W = sum(W,1);
        D =  diagmat(colsum_W); // diagnol matrix whose entries are column
        L = D -  phi*W;         // graph laplacian
        accu_L = accu(L);
        UtXV = UtX * V;
        VtV = V.t() * V;
        UtU = U.t() * U;
        diag_UtU = UtU.diag();

        // calculate initial objective function 
        normNMF = trac_xxt - 2.0 * trace(UtXV) + trace(UtU * VtV);
        logX = -(double)(mGene * nSample) * 0.5 * log(sigma_e2) - 0.5 * (double)(normNMF / sigma_e2);

        temp = (V.t() - b * vecOne.t()) * L * (V - vecOne * b.t());
        logV = - (double)(nSample) * 0.5 * sum(log(lambda )) - 0.5 * (sum(temp.diag() / lambda ));
        logSigmaL2 = -(alpha + 1.0) * sum(log(lambda)) - sum(beta / lambda);
        obj_old = logX + logV + logSigmaL2;
        V_old = V;

        // iteration starts
        for(int i = 1; i <= max_iter; ++i) {
            logV = 0.0;

            b = sum(V.t() * L, 1) / accu_L;
            lambda = (temp.diag() / 2.0 + beta ) / (double(nSample) / 2.0 + alpha + 1.0);  
            part1 = sigma_e2 * (D * V + phi * colsum_W * b.t());
            part2 = sigma_e2 * (phi * W * V + colsum_W * b.t());

            for(int nCT = 0; nCT < k; ++nCT){
                updateV_den_k = lambda(nCT) * (V.col(nCT) * diag_UtU(nCT) + (V * UtU.col(nCT) - V.col(nCT) * diag_UtU(nCT))) +  part1.col(nCT);
                updateV_k = (lambda(nCT) * XtU.col(nCT) + part2.col(nCT)) / updateV_den_k;
                V.col(nCT) %= updateV_k;
            }

            UtXV = UtX * V;
            VtV = V.t() * V;
            normNMF = trac_xxt - 2.0 * trace(UtXV) + trace(UtU * VtV);
            sigma_e2 = normNMF / (double)(mGene * nSample);
            temp = (V.t() - b * vecOne.t()) * L * (V - vecOne * b.t());
            logX = -(double)(nSample * mGene) * 0.5 * log(sigma_e2) - 0.5 * (double)(normNMF / sigma_e2);
            logV = - (double)(nSample) * 0.5 * sum(log(lambda))- 0.5 * (sum(temp.diag() / lambda )); 
            logSigmaL2 = -(alpha + 1.0) * sum(log(lambda)) - sum(beta / lambda);
            obj = logX + logV + logSigmaL2;
            logicalLogL = (obj > obj_old) && (abs(obj - obj_old) * 2.0 / abs(obj + obj_old) < epsilon);

            if(isnan(obj) || (sqrt(accu((V - V_old) % (V - V_old)) / double(nSample * k))  < epsilon) || logicalLogL){
                if(i > 5){ // run at least 5 iterations 
                    break;
                }
            } else {
                obj_old = obj;
                V_old = V;
            }
        }

        // Convert the arma::Mat<double> V to a nested std::vector<std::vector<double>>
        std::vector<std::vector<double>> V_vec = armaToVector(V);

        // Convert the arma::Mat<double> lambda to a nested std::vector<std::vector<double>>
        std::vector<std::vector<double>> lambda_vec = armaToVector(lambda);

        // Convert the arma::Mat<double> b to a nested std::vector<std::vector<double>>
        std::vector<std::vector<double>> b_vec = armaToVector(b);

        //Convert the result to a py::dict
        py::dict py_result;
        py_result["V"] = V_vec;
        py_result["sigma_e2"] = sigma_e2;
        py_result["lambda"] = lambda_vec;
        py_result["b"] = b_vec;
        py_result["Obj"] = obj;
        
        return py_result;

    } catch (std::exception &ex) {
        // Handle exceptions
        std::cerr << "C++ exception: " << ex.what() << std::endl;

    } catch (...) {
        std::cerr << "Unknown C++ exception occurred." << std::endl;
        throw; // Rethrow the caught exception
    }

    // If an exception occurs, return an empty dictionary
    return py::dict();  // Return an empty dict in case of an error
}
