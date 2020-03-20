//
// Created by vitowu on 3/21/20.
//

#include "mq.h"

__device__ static uint32_t __forceinline__
ctz(uint32_t c) {
    return __clz(__brev(c));
}

__device__ static void __forceinline__
diff_eq(uint8_t *func, uint32_t idx, uint8_t *result) {
    memset(result, 0x0, (MQ_VAR_NUM + 1) * sizeof(uint8_t));

    result[MQ_VAR_NUM] = func[MQ_XVAR_NUM - 1 - (MQ_VAR_NUM - idx)];

    uint32_t cursor = MQ_XVAR_NUM - 1 - MQ_VAR_NUM - 1;
    uint32_t i, j;
    uint32_t bound = (0 == idx) ? 1 : idx;
    for (i = MQ_VAR_NUM - 1; i >= bound; --i) {
        if (i == idx) {
            for (j = 1; j <= i; ++j) {
                result[i - j] ^= func[cursor - j];
            }
        } else {
            result[i] ^= func[cursor - (i - idx)];
        }
        cursor -= i + 1;
    }
}

__device__ static void __forceinline__
find_partial_derivs(uint8_t *mqsystem, uint8_t derivs[MQ_EQ_NUM][MQ_VAR_NUM][MQ_VAR_NUM + 1]) {
    uint32_t eq_idx, var_idx;
    for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
        for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
            diff_eq(mqsystem + eq_idx * MQ_XVAR_NUM, var_idx, derivs[eq_idx][var_idx]);
        }
    }
}

__device__ static void __forceinline__
reduce_sys(uint8_t *mqsystem) {
    uint32_t eq_idx, var_idx, i, sqr_term_idx;
    for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
        for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
            for (i = 0, sqr_term_idx = 0; i < var_idx; i++) {
                sqr_term_idx += i + 2;
            }
            mqsystem[eq_idx * MQ_XVAR_NUM + (MQ_XVAR_NUM - 2 - (MQ_VAR_NUM - 1 - var_idx))] ^=
                mqsystem[eq_idx * MQ_XVAR_NUM + sqr_term_idx];
            mqsystem[eq_idx * MQ_XVAR_NUM + sqr_term_idx] = 0;
        }
    }
}

__device__ static void __forceinline__
fast_exhaustive(uint8_t *mqsystem, uint8_t *solution) {
    uint64_t pdiff_eval[MQ_VAR_NUM];
    uint64_t func_eval = 0x0UL;
    uint64_t pre_fp_idx;
    uint32_t count = 0;
    uint32_t fp_idx;
    const uint32_t bound = (0x1U << MQ_VAR_NUM) - 1;
    uint64_t pdiff2[MQ_VAR_NUM][MQ_VAR_NUM];

    reduce_sys(mqsystem);
    uint8_t derivs[MQ_EQ_NUM][MQ_VAR_NUM][MQ_VAR_NUM + 1];
    find_partial_derivs(mqsystem, derivs);

    uint32_t eq_idx, var_idx, i, term;
    for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
        memset(pdiff2[var_idx], 0x0, sizeof(uint64_t) * MQ_VAR_NUM);
        for (i = 0; i < MQ_VAR_NUM; i++) {
            for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
                term = derivs[eq_idx][var_idx][i];
                pdiff2[i][var_idx] |= term << eq_idx;
            }
        }
    }

    memset(pdiff_eval, 0x0, sizeof(uint64_t) * MQ_VAR_NUM);
    for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
        for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
            if (var_idx == 0) {
                term = derivs[eq_idx][0][MQ_VAR_NUM];
            } else {
                term = derivs[eq_idx][var_idx][MQ_VAR_NUM] ^ derivs[eq_idx][var_idx][var_idx - 1];
            }
            pdiff_eval[var_idx] |= term << eq_idx;
        }
    }

    // brute forcing
    for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
        term = mqsystem[eq_idx * MQ_XVAR_NUM + MQ_XVAR_NUM - 1];
        func_eval |= term << eq_idx;
    }
    while (func_eval && count < bound) {
        count++;
        fp_idx = ctz(count);

        if (count & (count - 1)) {
            pre_fp_idx = ctz(count ^ (0x1U << fp_idx));
            pdiff_eval[fp_idx] ^= pdiff2[fp_idx][pre_fp_idx];
        }

        func_eval ^= pdiff_eval[fp_idx];
    }

    if (!func_eval) {
        for (var_idx = 0; var_idx < MQ_VAR_NUM; var_idx++) {
            solution[var_idx] = (uint8_t) (((count ^ count >> 1) >> var_idx) & 0x1U);
        }
    } else {
        memset(solution, 0x0, MQ_VAR_NUM * sizeof(uint8_t));
    }
}