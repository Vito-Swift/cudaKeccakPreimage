//
// Created by vitowu on 3/13/20.
//

#ifndef KECCAKSOLVER_KECCAK_H
#define KECCAKSOLVER_KECCAK_H

#include <cuda_runtime.h>
#include <inttypes.h>

#include "cuda_utils.h"

__device__
bool
gpu_VerifyKeccakResult(uint32_t A[25]);

/* function: ROR32
 * usage: perform right-rotate permutation
 */
__device__ __forceinline__
uint32_t
ROR32(uint32_t x, uint32_t a) {
    return (x >> a) | (x << (32 - a));
}

inline void
printHash(uint32_t hash[3]) {
    for (uint32_t i = 0; i < 3; i++) {
        PRINTF_STAMP("Hash[%d]: 0x%08x\n", i, hash[i]);
    }
}

inline void
printStatus(uint32_t A[5][5]) {
    uint32_t i, j;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 5; j++) {
            PRINTF_STAMP("A[%d][%d]: 0x%08x\n", i, j, A[i][j]);
        }
    }
}

bool
cpu_VerifyKeccakResult(const uint32_t A[5][5], uint32_t *minDiff);

void
cpu_DumpKeccakResult(const uint32_t A[5][5], uint32_t hash[3]);

void
cpu_VerifyRound2(const uint32_t A[5][5]);

#endif //KECCAKSOLVER_KECCAK_H
