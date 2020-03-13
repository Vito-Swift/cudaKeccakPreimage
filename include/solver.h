//
// Created by vitowu on 3/13/20.
//

#ifndef KECCAKSOLVER_SOLVER_H
#define KECCAKSOLVER_SOLVER_H

#include "options.h"
#include "cuda_utils.h"

typedef struct KeccakSolver {

  bool verbose;

  Options options;

  uint8_t* device_linsys_buffer;
  uint8_t *device_mq_buffer;
  uint8_t *device_c_constr_buffer;
  uint32_t *device_output_buffer;

} KeccakSolver;

void keccakSolverInit(KeccakSolver *keccakSolver, int argc, char** argv);
void keccakSolverLoop(KeccakSolver *keccakSolver);
void keccakSolverExit(KeccakSolver *keccakSolver);

#endif //KECCAKSOLVER_SOLVER_H
