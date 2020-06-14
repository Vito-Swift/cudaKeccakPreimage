//
// Created by vitowu on 3/13/20.
//

#ifndef KECCAKSOLVER_SOLVER_H
#define KECCAKSOLVER_SOLVER_H

#include "options.h"
#include "cuda_utils.h"
#include "threadpool.h"
#include "kmath.h"

typedef struct KeccakSolver {
    bool verbose;
    Options options;
    MathSystem mathSystem;
    uint8_t *device_mq_buffer;
    uint8_t *device_output_buffer;
    threadpool_t *threadpool;
    uint32_t **test_message;
} KeccakSolver;

typedef struct tmqarg_t {
    uint64_t guessingBits;
    uint8_t *mqbuffer;
    uint32_t *mq2lin;
    uint32_t *lin2mq;
    uint8_t *lindep;
    MathSystem *mathSystem;
} tmqarg_t;

typedef struct tcheckarg_t {
    uint8_t *result_buffer;
    uint8_t *lindep;
    uint32_t *mq2lin;
    uint8_t *mqbuffer;
    uint32_t *lin2mq;
    uint32_t *minDiff;
    bool *preimage_found;
} tcheckarg_t;

void keccakSolverInit(KeccakSolver *keccakSolver, int argc, char **argv);

void keccakSolverLoop(KeccakSolver *keccakSolver);

void keccakSolverExit(KeccakSolver *keccakSolver);

#endif //KECCAKSOLVER_SOLVER_H
