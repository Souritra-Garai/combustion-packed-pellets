#include "TDMA_solver.hpp"

using namespace TDMA_solver;

solver::solver(int m)
{
    n = m;

    Upper_Matrix_Diagonal.resize(n, 0.0);
    Upper_Matrix_Diagonal_plus_1.resize(n, 0.0);
    Lower_Matrix_Diagonal_less_1.resize(n, 0.0);

    ready_to_solve = false;
}

void solver::LU_Decomposition( vector<float> E, vector<float> F, vector<float> G)
{
    Upper_Matrix_Diagonal[0] = F[0];
    Upper_Matrix_Diagonal_plus_1[0] = G[0];
    
    for (int i=1; i<n; i++)
    {
        Lower_Matrix_Diagonal_less_1[i] = E[i] / Upper_Matrix_Diagonal[i-1];

        Upper_Matrix_Diagonal[i] = F[i] - G[i-1]*Lower_Matrix_Diagonal_less_1[i];

        Upper_Matrix_Diagonal_plus_1[i] = G[i];
    }
}

void solver::setup_banded_matrix(   vector<float> E,
                                    vector<float> F,
                                    vector<float> G   )
{
    if (E.size() != n) throw "Size of E must be equal to n";
    if (F.size() != n) throw "Size of F must be equal to n";
    if (G.size() != n) throw "Size of G must be equal to n";

    LU_Decomposition(E, F, G);

    for (int i=0; i<n; i++)

        if (Upper_Matrix_Diagonal[i] == 0)

            throw "Upper Matrix main diagonal has a zero!!";

    ready_to_solve = true;    
}

vector<float> solver::solve_Ldb(vector<float> b)
{
    vector<float> d(n, 0.0);

    d[0] = b[0];

    for (int i=1; i<n; i++)

        d[i] = b[i] - d[i-1] * Lower_Matrix_Diagonal_less_1[i];

    return d;
}

vector<float> solver::solve_Uxd(vector<float> d)
{
    vector<float> x(n, 0.0);

    x[n-1] = d[n-1] / Upper_Matrix_Diagonal[n-1];

    for (int i=n-2; i>=0; i--)

        x[i] = ( d[i] - Upper_Matrix_Diagonal_plus_1[i]*x[i+1] ) / Upper_Matrix_Diagonal[i];

    return x;    
}

vector<float> solver::solve(vector<float> b)
{
    if (b.size() != n) throw "Size of b must be equal to n";

    if (!ready_to_solve) throw "Solver has not been set up yet!!";

    return solve_Uxd(solve_Ldb(b));
}

vector<vector<float>> solver::solve(vector<vector<float>> B)
{
    vector<vector<float>> X(B.size(), vector<float>(n, 0.0));

    #pragma omp parallel for schedule(static) default(none) shared(B, X)

        for (int i=0; i<B.size(); i++)

            X[i] = solve(B[i]);

    return X;
}