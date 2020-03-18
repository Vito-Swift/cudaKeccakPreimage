//
// Created by vitowu on 3/18/20.
//

#ifndef KECCAKSOLVER_MATH_H
#define KECCAKSOLVER_MATH_H

#define kidx(i, j, k) \
    ((((5 * i) + j) * 32) + k)

#define cbinom2(n) \
    ( ((n) * ((n)+1)) / 2)

#define deg2midx1(vnum, var1_idx) \
    (cbinom2(vnum) + (var1_idx))

// pre-requirements: var2_idx > var1_idx
#define deg2midx2(var1_idx, var2_idx) \
    ((var2_idx) > (var1_idx) ? (cbinom2(var2_idx) + (var1_idx)) : (cbinom2(var1_idx) + (var2_idx)))

#include <stdint.h>
#include <cstring>
#include <iostream>

typedef struct MathSystem {
  uint32_t *round3_lin_dep;
  uint32_t *round3_mq2lin;
  uint32_t *round3_lin2mq;
  uint8_t *round3_append_system;
  uint8_t *round3_mq_system;

  uint32_t *round4_lin_dep;
} MathSystem;

void
initMathSystem(MathSystem *system);

void
freeMathSystem(MathSystem *system);

#endif //KECCAKSOLVER_MATH_H
