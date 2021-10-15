/*
The C++ class of the classical Jacobi's method implementation for finding eigenvalues and
eigenvectors of a symetric matrix which represented as 1-dimension array.
*/

#include <vector>

#ifndef JACOBIEIGENCLASS_H
#define JACOBIEIGENCLASS_H

class JacobiEigensClass {

private:
    const double epsilon;
    const int N;
    int maxIterations = 100;
    int totalIterations;
    int currentIteration;
    int currentIndexOfMaxElement;
    double currentMaxElement;

    double* varVector;
    //double *rotVector;

public:
    double* originVector = NULL;    //input, original matrix in array representation
    double* eigenVecs = NULL;       //output array with eigenvalues
    double* eigenVals = NULL;       //output array with eigenvector
    int num_of_valu = 0;

    JacobiEigensClass(const std::vector<double>& inputMatrix, const double epsilon);
    ~JacobiEigensClass();
    void calculateEigensByJacobiMethod();

private:
    void findMaxElementInVarVector();
    void rotateVarVector();

};

#endif