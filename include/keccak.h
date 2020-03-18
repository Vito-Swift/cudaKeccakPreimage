//
// Created by vitowu on 3/13/20.
//

#ifndef KECCAKSOLVER_KECCAK_H
#define KECCAKSOLVER_KECCAK_H

#include <cuda_runtime.h>
#include <inttypes.h>

#include "cuda_utils.h"

__device__ static __forceinline__
bool gpu_VerifyKeccakResult(uint32_t A[25]);

/* function: ROR32
 * usage: perform right-rotate permutation
 */
__device__ __forceinline__
uint32_t
ROR32(uint32_t x, uint32_t a) {
    return (x >> a) | (x << (32 - a));
}

bool cpu_VerifyKeccakResult(const uint32_t A[5][5]);

#endif //KECCAKSOLVER_KECCAK_H
