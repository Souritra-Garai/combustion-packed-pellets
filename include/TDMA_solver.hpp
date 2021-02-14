#ifndef __TDMA_SOLVER__
#define __TDMA_SOLVER__

#include <vector>

using namespace std;

namespace TDMA_solver
{
    class solver
    {
        private:

            // size of x vector for the linear algebra problem A.x = b
            size_t n;

            // vectors / arrays carrying values of diagonals 
            // of decomposed upper and lower matrices
            
            // Upper matrix
            // Main Diagonal elements - total size n
            vector<float> Upper_Matrix_Diagonal;
            // Diagonal immediately to right of main diagonal
            // total size n; Upper_Matrix_Diagonal_plus_1[n-1] is ignored
            vector<float> Upper_Matrix_Diagonal_plus_1;
            
            // Lower matrix
            // Diagonal innediately to left of main diagonal
            // total size n; Lower_Matrix_Diagonal_less_1[0] is ignored
            vector<float> Lower_Matrix_Diagonal_less_1;

            // flag to specify whether everything has been
            // initialised for solving
            bool ready_to_solve;

            // helper funcitons to solve the matrix equation in two steps
            vector<float> solve_Ldb(vector<float> b);
            vector<float> solve_Uxd(vector<float> d);

            // helper function to decompose the A matrix in Upper and Lower matrices
            void LU_Decomposition(vector<float>, vector<float>, vector<float>);

        public:

            // Constructor
            // m is the size of x vector in A.x = b
            solver(int m);

            // function to initialise the solver
            // F is the main diagonal of A
            // E is the diagonal immediately left of the main diagonal
            // G is the diagonal immediately right of the main diagonal
            // F, E and G all have n elements
            // E[0] is ignored; G[n-1] is ignored
            void setup_banded_matrix(   vector<float> E,
                                        vector<float> F,
                                        vector<float> G   );
            
            // overloaded function solve to solve the equation A.x = b
            // b may be 1D vector or 2D vector with n rows
            vector<float> solve(vector<float> b);
            vector<vector<float>> solve(vector<vector<float>> B);
    };
}

#endif