//
// Created by vitowu on 3/18/20.
//

#ifndef KECCAKSOLVER_MATH_H
#define KECCAKSOLVER_MATH_H

#include <stdint.h>
#include <cstring>

typedef struct MathSystem {
  uint32_t* round3_lin_dep;
  uint8_t round3_append_system;
} MathSystem;

#endif //KECCAKSOLVER_MATH_H
