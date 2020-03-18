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

#include "params.h"

typedef struct MathSystem {
  uint8_t *round3_lin_dep[800];
  uint32_t *round3_mq2lin;
  uint32_t *round3_lin2mq;
  uint8_t *round3_append_system[AMQ_LIN_EQNUM];
  uint8_t *round3_mq_system[MQ_EQ_NUM];
  uint8_t *round3_iter_system[LIN_ITER_EQNUM];
  uint8_t *round4_lin_dep;
} MathSystem;

void
initMathSystem(MathSystem *system);

void
extractRound3LinearDependency(MathSystem *system, uint8_t lin_system[LIN_CONST_EQNUM][801]);

void
reduceRound3AppendSystem(MathSystem *system, uint8_t **append_system);

void
reduceRound3MQSystem(MathSystem *system, uint8_t **mqsystem);

void
reduceIterativeConstraints(MathSystem *system, uint8_t iterative_constr[LIN_ITER_EQNUM][801]);

void
guessingBitsToMqSystem(MathSystem *system,
                       uint64_t guessingBits,
                       uint8_t *mqbuffer,
                       uint32_t *mq2lin,
                       uint8_t *lin_dep);

void
freeMathSystem(MathSystem *system);

#endif //KECCAKSOLVER_MATH_H
