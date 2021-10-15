//test file

#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
using std::cerr;
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
#include "jacobieigenclass.h"
#include <string>
#include <cstdlib>
#include <sstream>
#include <chrono>   //for time 

#define eps 0.0001  //for jacobi method



using namespace std;

int main() {

    //get the time that had be moved at current time
    auto t1 = std::chrono::high_resolution_clock::now();
    
    ifstream inputdata;
    string line;
    vector<double> input;
    int row = 0;   //number of rows in input matrix
    int col = 0;   //number of colmuns in input matrix
    inputdata.open("data_8.txt",ios::in); // opens the file
    if (!inputdata)
    { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    while (getline(inputdata, line))
    {
        col = 0;
       istringstream buf(line);
       while (!buf.eof())
       {
           double float_num;
           buf >> float_num;
           if (buf.bad())
           {
               std::cerr << "Number formatting error!\n";
               exit(1);
           }
           input.push_back(float_num);
           char comma = 0;
           buf >> comma;
           if (!buf.eof())
           {
               if (buf.fail())
               {
                   std::cerr << "Could not read comma.\n";
                   exit(1);
               }
               if (comma != ',')
               {
                   std::cerr << "Found no comma but '" << comma << "' instead !\n";
                   return 1;
               }
           }
           col++;
       }
       row++;
    }

    int NX = input.size();


    double** A = new double* [row];
    for (int i = 0; i < row; ++i)
        A[i] = new double[col];

    int k = 0;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            A[i][j] = input[k];
            k++;
        }
    }
    
    


    
    
    //********* for test a simple matrix *****************

    /*
    double A[2][4] = {1,0,1,0,0,1,0,1};
    int row = (sizeof(A) / sizeof(A[0]));
    int col = (sizeof(A) / sizeof(A[0][0])) / row;
    */
    
    double** AT;
    double** ATA;
    double** AAT;
    double** segma;


    //define sema matrix,
    segma = new double* [row];
    for (int i = 0; i < row; ++i)
        segma[i] = new double[col];

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            segma[i][j] = 0;
        }
    }

    
    AT = new double* [col];
    for (int i = 0; i < col; ++i)
        AT[i] = new double[row];
    //define a transpze
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j)
        {
            AT[j][i] = A[i][j];
        }

    //define AAT MATRIX
    AAT = new double* [row];
    for (int i = 0; i < row; ++i)
        AAT[i] = new double[row];

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < row; ++j) {
            AAT[i][j] = 0;
        }
    }

    //get AAT matrix by A * AT
    for (int i = 0; i < row; i++) {
        for (int j= 0; j < row; j++) {
            AAT[i][j] = 0;
            for (int k = 0; k < col; k++)
                AAT[i][j] += A[i][k] * AT[k][j];
 
        }
    }


    vector<double> vecAAT;
    //Push AAT matrix in vector for eigen operation
    for(int i = 0; i < row ;i++)
        for (int j = 0; j < row; j++)
        {
            vecAAT.push_back(AAT[i][j]);
        }

    //create a new class of JacobiEigensClass ,,
    JacobiEigensClass* obj = new JacobiEigensClass(vecAAT, eps);

    //get the eigen valu and eigen vectors for AAT matrix 
    obj->calculateEigensByJacobiMethod();

    int N = sqrt(vecAAT.size());

    cout << endl << "Eigenvectors for MATRIX U: " << endl;

    //print eigen vector for Matrix U
    for (int i = 0; i < N; i++)
    {
        cout << "[";
        for (int j = 0; j < N; j++)
        {
            cout << (double)obj->eigenVecs[i * N + j] << ", ";
        }
        cout << "]" << endl;
    }
    
    //set eigen valu in segma vector
    int k = 0;
    vector<double> segmavec1;
    //int s = sizeof(obj->eigenVals)/sizeof(obj->eigenVals[0]);
    for (int i = 0; i < N; i++)
    {
        if (obj->eigenVals[k] != 0)
            segmavec1.push_back(obj->eigenVals[k]);
        k++;
    }
    //set the singular value in the segma matrix
    k = 0;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if ((i == j) && (k < segmavec1.size()))
            {
                segma[i][j] = sqrt(segmavec1[k]);
                k++;
            }
        }
    }


    cout << "----------------------------------------------------------------------" << endl;
    
    //define ATA matrix
    ATA = new double* [col];
    for (int i = 0; i < col; ++i)
        ATA[i] = new double[col];

    //set ATA matrix = {0}
    for (int i = 0; i < col; ++i) {
        for (int j = 0; j < col; ++j) {
            ATA[i][j] = 0;
        }
    }

    //get ATA by A * AT
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
            ATA[i][j] = 0;
            for (int k = 0; k < row; k++)
                ATA[i][j] += AT[j][k] * A[k][i]  ;

        }
    }


    vector<double> vecATA;
    //set tha ATA matrix in vector for eigen operation 
    for (int i = 0; i < col; i++)
        for (int j = 0; j < col; j++)
        {
            vecATA.push_back(ATA[i][j]);
        }
    //create a new class of jacobe class method
    JacobiEigensClass* obj1 = new JacobiEigensClass(vecATA, 0.0001);

    //get the eigen value and eigen vector for the matrix ATA
    obj1->calculateEigensByJacobiMethod();

    int N1 = sqrt(vecATA.size());

    cout << endl << "Eigenvectors for MATRIX V: " << endl;

    //print the vectors of the matrix V
    for (int i = 0; i < N1; i++)
    {
        cout << "[";
        for (int j = 0; j < N1; j++)
        {
            cout << (double)obj1->eigenVecs[i * N1 + j] << ", ";
        }
        cout << "]" << endl;
    }

    //set eigen valu in segma vector
    k = 0;
    vector<double> segmavec2;
    //int s1 = sizeof(obj1->eigenVals) / sizeof(obj1->eigenVals[0]);
    for (int i = 0; i < N1; i++)
    {
        if (obj1->eigenVals[k] != 0)
            segmavec2.push_back(obj1->eigenVals[k]);
        k++;
    }

    //check if there a new sigular value of the segma matrix and set them in the segma matrix
    k = 0;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if ((i == j) && segmavec2[k] != 0 && k < segmavec2.size())
            {
                if (segma[i][j] != sqrt(segmavec2[k]))
                {
                    segma[i][j] = sqrt(segmavec2[k]);
                }
                k++;
            }
            
        }
    }

    cout << "---------------------------------Print Segma Matrix---------------------------------------"<<endl;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            cout << (double)segma[i][j];
            cout << ',';
        }
        cout <<""<< endl;
    }

    // get the time that had be moved from last time 
    auto t2 = std::chrono::high_resolution_clock::now();
    //get the time of that the project have to find U,SEGMA AND V
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    std::cout << duration;

    cout << endl;

    return 0;
}