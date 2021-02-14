#include <iostream>
#include <stdio.h>
#include <vector>
#include <omp.h>

#include "TDMA_solver.h"

#define NUM_THREADS 8

int main()
{
    vector<float> E{1, 1};
    vector<float> F{1,-1};
    vector<float> G{1, 1};

    vector<vector<float>> b{ {10, 5}, {2, 0}, {8, 4}, {7, 9} };
    vector<vector<float>> x;

    omp_set_num_threads(NUM_THREADS);

    try
    {
        TDMA_solver::solver mySolver(2);

        mySolver.setup_banded_matrix(E, F, G);

        x = mySolver.solve(b);

        for (int i=0; i<2; i++) 
        {
            for (int j=0; j<b.size(); j++) printf("%2.1f\t", x[j][i]);
            printf("\n");
        }
    }
    catch(const char * error_msg)
    {
        cerr << error_msg << endl;
    }

    return 0;
}