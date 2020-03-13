/**
 * @desc: Keccak 4-round pre-image attack implementation
 *  ╦╔═┌─┐┌─┐┌─┐┌─┐┬┌─  ╔═╗┌─┐┬ ┬  ┬┌─┐┬─┐
 *  ╠╩╗├┤ │  │  ├─┤├┴┐  ╚═╗│ ││ └┐┌┘├┤ ├┬┘
 *  ╩ ╩└─┘└─┘└─┘┴ ┴┴ ┴  ╚═╝└─┘┴─┘└┘ └─┘┴└─
 * @language: C++/CUDA
 */

#include <cuda_runtime.h>

#include <iostream>
#include <stdlib.h>
#include <inttypes.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>

#include "keccak.h"
#include "solver.h"

int main(int argc, char **argv) {
    KeccakSolver keccakSolver;

    keccakSolverInit(&keccakSolver, argc, argv);
    keccakSolverLoop(&keccakSolver);
    keccakSolverExit(&keccakSolver);

    return 0;
}
