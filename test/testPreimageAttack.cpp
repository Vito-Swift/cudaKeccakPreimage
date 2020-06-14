//
// Created by vitowu on 6/5/20.
//

#include "keccak.h"
#include "solver.h"

int main(int argc, char **argv) {
    KeccakSolver solver;
    keccakSolverInit(&solver, argc, argv);
    keccakSolverLoop(&solver);
    keccakSolverExit(&solver);
}
