//
// Created by vitowu on 3/14/20.
//

#include <iostream>
#include <stdint.h>

#include "keccak.h"

__device__ static __forceinline__
uint32_t
ctz(uint32_t c) {
    return __clz(__brev(c));
}

__global__ void
testMqTime() {
    uint32_t count = 0;
    uint32_t bound = (1U << 31) - 1;
    uint32_t fp_idx, pre_fp_idx;
    uint32_t a[31];
    uint32_t func_eval = 0x0U;

    while (count < bound) {
        if (count == 1 << 28)
            printf("hi\n");
        count++;
        fp_idx = ctz(count);
        if (count & (count - 1)) {
            pre_fp_idx = ctz(count ^ (1U << fp_idx));
            a[fp_idx] ^= pre_fp_idx;
        }
        func_eval ^= 1;
    }
}

int main() {
//    dim3 tpb(64);
//    cudaDeviceSynchronize();
//    testMqTime<<<64, tpb>>>();
//    cudaDeviceSynchronize();
    uint32_t A[5][5];
    cpu_VerifyKeccakResult(A);
    return 0;
}