//
// Created by vitowu on 3/13/20.
//

#include "solver.h"

__host__ void
keccakSolverInit(KeccakSolver *keccakSolver, int argc, char **argv) {
    options_init(&keccakSolver->options);
    options_parse(&keccakSolver->options, argc, argv);
}

__host__ void
keccakSolverLoop(KeccakSolver *keccakSolver) {

}

__host__ void
keccakSolverExit(KeccakSolver *keccakSolver) {
    options_free(&keccakSolver->options);
}