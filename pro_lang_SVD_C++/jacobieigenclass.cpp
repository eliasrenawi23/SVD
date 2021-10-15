#include <cmath>
#include <vector>

#include "jacobieigenclass.h"

//using namespace std;   
JacobiEigensClass::JacobiEigensClass(const std::vector<double>& inputMatrix, const double epsilon) : epsilon(epsilon), N(static_cast<int>(sqrt(inputMatrix.size()))) {
    maxIterations = N * N * N;
    originVector = new double[N * N];
    eigenVecs = new double[N * N];
    eigenVals = new double[N];
    for (int n = 0; n < N * N; n++) originVector[n] = inputMatrix[n];

}

JacobiEigensClass::~JacobiEigensClass() {
    delete[] originVector;
    delete[] eigenVecs;
    delete[] eigenVals;
}

void JacobiEigensClass::calculateEigensByJacobiMethod() {

    varVector = new double[N * N];
    currentIteration = 0;

    for (int n = 0; n < N * N; n++) varVector[n] = originVector[n];

    // Setting up the eigenvector matrix
    for (int i = 0; i < N * N; i++) {
        eigenVecs[i] = 0;
    }
    for (int i = 0; i < N; i++) {
        eigenVecs[i * N + i] = 1.0;
    }

    findMaxElementInVarVector();
    while (fabs(varVector[currentIndexOfMaxElement]) > epsilon && currentIteration < maxIterations) {
        rotateVarVector();
        findMaxElementInVarVector();
        currentIteration++;
    }
    for (int n = 0; n < N; n++) eigenVals[n] = varVector[n * N + n];
    delete[] varVector;
}

// Function to find the maximum matrix element
void JacobiEigensClass::findMaxElementInVarVector() {
    currentMaxElement = 0;
    int indexOfMaxElement = 0;
    int vectorLength = N * N;
    int n = 1;
    for (int i = 1; i < vectorLength; i++) {
        if (fabs(currentMaxElement) < fabs(varVector[i])) {
            currentMaxElement = varVector[i];
            currentIndexOfMaxElement = i;
        }
        if (i == n * N - 1) {
            n++;
            i = i + n;
        }
    }
}

//function to do operation rotation
void JacobiEigensClass::rotateVarVector() {
    double c, s;
    int i_max = static_cast<int>(currentIndexOfMaxElement / N);
    int j_max = currentIndexOfMaxElement - static_cast<int>(currentIndexOfMaxElement / N) * N;

    //Angle calculations.
    // Function to find the values of cos and sin
    if (varVector[j_max * N + i_max] != 0.0) {
        double tau, t;
        tau = (varVector[j_max * N + j_max] - varVector[i_max * N + i_max]) / (2.0 * varVector[j_max * N + i_max]);
        if (tau > 0) {
            t = 1.0 / (tau + sqrt(1.0 + tau * tau));
        }
        else {
            t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
        }
        c = 1 / sqrt(1 + t * t);
        s = c * t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }

    //The rotor (eigenVecs) accumulation
    double* tmpRowIFromEigenVecs = new double[N];
    double* tmpRowJFromEigenVecs = new double[N];

    for (int j = 0; j < N; j++) {
        tmpRowIFromEigenVecs[j] = eigenVecs[i_max * N + j];
        tmpRowJFromEigenVecs[j] = eigenVecs[j_max * N + j];

        eigenVecs[i_max * N + j] = tmpRowIFromEigenVecs[j] * c - tmpRowJFromEigenVecs[j] * s;
        eigenVecs[j_max * N + j] = tmpRowIFromEigenVecs[j] * s + tmpRowJFromEigenVecs[j] * c;
    }

    delete[] tmpRowIFromEigenVecs;
    delete[] tmpRowJFromEigenVecs;
    

    // varVector calculations
    double a_ii = varVector[i_max * N + i_max];
    double a_jj = varVector[j_max * N + j_max];
    double a_il, a_jk;

    // changing the matrix elements with indices i_max and j_max
    varVector[i_max * N + i_max] = c * c * a_ii - 2.0 * s * c * varVector[i_max * N + j_max] + s * s * a_jj;
    varVector[j_max * N + j_max] = s * s * a_ii + 2.0 * s * c * varVector[i_max * N + j_max] + c * c * a_jj;
    varVector[i_max * N + j_max] = 0.0;
    varVector[j_max * N + i_max] = 0.0;

    // change the remaining elements
    for (int l = 0; l < N; l++) {
        if (l != i_max && l != j_max) {
            a_il = varVector[i_max * N + l];
            a_jk = varVector[j_max * N + l];

            varVector[i_max * N + l] = c * a_il - s * a_jk;
            varVector[l * N + i_max] = varVector[i_max * N + l];

            varVector[j_max * N + l] = s * a_il + c * a_jk;
            varVector[l * N + j_max] = varVector[j_max * N + l];
        }
    }
}