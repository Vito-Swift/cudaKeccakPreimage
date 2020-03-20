//
// Created by vitowu on 3/21/20.
//

#ifndef KECCAKSOLVER_MQ_H
#define KECCAKSOLVER_MQ_H

#include <cuda_runtime.h>

#include "cuda_utils.h"
#include "params.h"

__device__ static void __forceinline__
fast_exhaustive(uint8_t* mqsystem, uint8_t* solution);

#endif //KECCAKSOLVER_MQ_H
