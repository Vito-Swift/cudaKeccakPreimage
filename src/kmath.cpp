//
// Created by vitowu on 3/18/20.
//

#include "kmath.h"

#include "params.h"
#include "utils.h"

#define eqvar(i, j, size) \
        ((i) * (size) + (j))

#define _VERIFICATION

void
initMathSystem(MathSystem *system) {
    uint32_t i;

    // linear dependency: x_k = x_n_1 + x_n_2 + ... + x_n_j
    for (i = 0; i < 800; i++) {
        system->round3_lin_dep[i] = (uint8_t *) malloc((IMQ_VAR_NUM + 1) * sizeof(uint8_t));
        memset(system->round3_lin_dep[i], 0, (IMQ_VAR_NUM + 1) * sizeof(uint8_t));
    }

    // append system in order 2
    for (i = 0; i < AMQ_LIN_EQNUM; i++) {
        system->round3_append_system[i] = (uint8_t *) malloc(IMQ_XVAR_NUM * (sizeof(uint8_t)));
        memset(system->round3_append_system[i], 0, IMQ_XVAR_NUM * sizeof(uint8_t));
    }

    // mq system in order 2
    for (i = 0; i < MQ_EQ_NUM; i++) {
        system->round3_mq_system[i] = (uint8_t *) malloc(IMQ_XVAR_NUM * sizeof(uint8_t));
        memset(system->round3_append_system[i], 0, IMQ_XVAR_NUM * sizeof(uint8_t));
    }

    // iterative system in order 1
    for (i = 0; i < LIN_ITER_EQNUM; i++) {
        system->round3_iter_system[i] = (uint8_t *) malloc((IMQ_VAR_NUM + 1) * sizeof(uint8_t));
        memset(system->round3_iter_system[i], 0, (IMQ_VAR_NUM + 1) * sizeof(uint8_t));
    }

    // transfer mq var index (94) to linear var index (94)
    system->round3_mq2lin = (uint32_t *) malloc(IMQ_VAR_NUM * sizeof(uint32_t));
    memset(system->round3_mq2lin, 0, IMQ_VAR_NUM * sizeof(uint32_t));

    // reversely transfer linear index to mq var index
    system->round3_lin2mq = (uint32_t *) malloc(800 * sizeof(uint32_t));
    memset(system->round3_lin2mq, 0, 800 * sizeof(uint32_t));

    // linear dependency
    system->round4_lin_dep = (uint8_t *) malloc(IMQ_VAR_NUM * (AMQ_VAR_NUM + 1) * sizeof(uint8_t));
    memset(system->round4_lin_dep, 0, IMQ_VAR_NUM * (AMQ_VAR_NUM + 1) * sizeof(uint8_t));
}

void
extractRound3LinearDependency(MathSystem *system, uint8_t lin_system[LIN_CONST_EQNUM][801]) {
    uint8_t tmprow[801];
    uint8_t lin_pivot[800] = {0};
    for (uint32_t i = 0; i < LIN_CONST_EQNUM; i++) {
        uint32_t j = 0;
        uint32_t k = 0;

        // index pivot
        for (j = 0; j < 800; j++) {
            for (k = i; k < LIN_CONST_EQNUM; k++) {
                if (lin_system[k][j])
                    break;
            }
            if (k != LIN_CONST_EQNUM)
                break;
        }
        if (j == 800)
            continue;

        // row swapping
        memcpy(tmprow, lin_system[k], 801);
        memcpy(lin_system[k], lin_system[i], 801);
        memcpy(lin_system[i], tmprow, 801);

        lin_pivot[j] = 1;

        for (k = 0; k < LIN_CONST_EQNUM; k++) {
            if (lin_system[k][j] && k != i)
                for (uint32_t v = j; v < 801; v++)
                    lin_system[k][v] ^= lin_system[i][v];
        }
    }
    uint32_t fvar_offset = 0;
    for (uint32_t i = 0; i < 800; i++) {
        if (lin_pivot[i] == 0) {
            system->round3_mq2lin[fvar_offset] = i;
            system->round3_lin2mq[i] = fvar_offset;
            fvar_offset++;
        } else {
            system->round3_lin2mq[i] = DEP_PLACEMENT;
        }
    }

    // note: must do index transferring first before calc linear dependency
    for (uint32_t i = 0; i < LIN_CONST_EQNUM; i++) {
        for (uint32_t j = 0; j < 800; j++) {
            if (lin_system[i][j] == 1) {
                for (uint32_t k = j + 1; k < 800; k++) {
                    if (lin_system[i][k] == 1) {
                        system->round3_lin_dep[j][system->round3_lin2mq[k]] = 1;
                    }
                }
                system->round3_lin_dep[j][IMQ_VAR_NUM] = lin_system[i][800];
                break;
            }
        }
    }
}

void
reduceRound3AppendSystem(void *r3aparg) {
    r3aparg_t *args = (r3aparg_t *) r3aparg;
    MathSystem *system = args->mathSystem;
    uint32_t eq_idx = args->eq_idx;
    uint8_t **append_system = args->append_system;

    uint32_t var_idx, i, j;
    uint32_t multix_1;
    uint32_t multix_2;

    uint8_t **reduced_system = system->round3_append_system;
    uint8_t **lin_dep = system->round3_lin_dep;
    uint32_t *lin2mq = system->round3_lin2mq;

    PRINTF_STAMP("\t\treducing append system:\tequation %d...\n", eq_idx);
    memset(reduced_system[eq_idx], 0, IMQ_XVAR_NUM);
    for (multix_1 = 0; multix_1 < 800; multix_1++) {
        bool is_multix1_dep = (bool) (lin2mq[multix_1] == DEP_PLACEMENT);

        for (multix_2 = multix_1; multix_2 < 800; multix_2++) {
            var_idx = deg2midx2(multix_1, multix_2);
            if (append_system[eq_idx][var_idx]) {
                bool is_multix2_dep = (bool) (lin2mq[multix_2] == DEP_PLACEMENT);

                if (is_multix1_dep && is_multix2_dep) {
                    for (i = 0; i < IMQ_VAR_NUM; i++) {
                        for (j = 0; j < IMQ_VAR_NUM; j++) {
                            if (lin_dep[multix_1][i] && lin_dep[multix_2][j])
                                reduced_system[eq_idx][deg2midx2(i, j)] ^= 1;
                        }
                    }

                    if (lin_dep[multix_1][IMQ_VAR_NUM]) {
                        for (i = 0; i < IMQ_VAR_NUM; i++) {
                            if (lin_dep[multix_2][i])
                                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, i)] ^= 1;
                        }
                    }

                    if (lin_dep[multix_2][IMQ_VAR_NUM]) {
                        for (i = 0; i < IMQ_VAR_NUM; i++) {
                            if (lin_dep[multix_1][i])
                                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, i)] ^= 1;
                        }
                    }

                    if (lin_dep[multix_1][IMQ_VAR_NUM] && lin_dep[multix_2][IMQ_VAR_NUM])
                        reduced_system[eq_idx][IMQ_XVAR_NUM - 1] ^= 1;

                } else if (is_multix1_dep) {
                    uint32_t tmp_idx = system->round3_lin2mq[multix_2];
                    for (i = 0; i < IMQ_VAR_NUM; i++) {
                        if (lin_dep[multix_1][i]) {
                            reduced_system[eq_idx][deg2midx2(tmp_idx, i)] ^= 1;
                        }
                    }

                    if (lin_dep[multix_1][IMQ_VAR_NUM])
                        reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, tmp_idx)] ^= 1;

                } else if (is_multix2_dep) {
                    uint32_t tmp_idx = system->round3_lin2mq[multix_1];
                    for (i = 0; i < IMQ_VAR_NUM; i++) {
                        if (lin_dep[multix_2][i]) {
                            reduced_system[eq_idx][deg2midx2(tmp_idx, i)] ^= 1;
                        }
                    }

                    if (lin_dep[multix_2][IMQ_VAR_NUM])
                        reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, tmp_idx)] ^= 1;

                } else {
                    uint32_t tmp_idx1 = system->round3_lin2mq[multix_1];
                    uint32_t tmp_idx2 = system->round3_lin2mq[multix_2];
                    reduced_system[eq_idx][deg2midx2(tmp_idx1, tmp_idx2)] ^= 1;
                }
            }
        }
        if (append_system[eq_idx][deg2midx1(800, multix_1)]) {
            if (is_multix1_dep) {
                for (i = 0; i < IMQ_VAR_NUM + 1; i++) {
                    if (lin_dep[multix_1][i])
                        reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, i)] ^= 1;
                }
            } else {
                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, lin2mq[multix_1])] ^= 1;
            }
        }
    }
    reduced_system[eq_idx][IMQ_XVAR_NUM - 1] ^= append_system[eq_idx][BMQ_XVAR_NUM - 1];
}

void
reduceRound3MQSystem(void *r3mqarg) {
    r3mqarg_t *args = (r3mqarg_t *) r3mqarg;
    MathSystem *system = args->mathSystem;
    uint32_t eq_idx = args->eq_idx;
    uint8_t **mqsystem = args->mqsystem;

    uint32_t var_idx, i, j;

    uint8_t **reduced_system = system->round3_mq_system;
    uint8_t **lin_dep = system->round3_lin_dep;
    uint32_t *lin2mq = system->round3_lin2mq;

    PRINTF_STAMP("\t\treducing mq system:\tequation %d...\n", eq_idx);
    memset(reduced_system[eq_idx], 0, IMQ_XVAR_NUM);
    uint32_t multix_1;
    uint32_t multix_2;
    for (multix_1 = 0; multix_1 < 800; multix_1++) {
        bool is_multix1_dep = (bool) (lin2mq[multix_1] == DEP_PLACEMENT);

        for (multix_2 = multix_1; multix_2 < 800; multix_2++) {
            var_idx = deg2midx2(multix_1, multix_2);
            if (mqsystem[eq_idx][var_idx]) {
                bool is_multix2_dep = (bool) (lin2mq[multix_2] == DEP_PLACEMENT);

                if (is_multix1_dep && is_multix2_dep) {
                    for (i = 0; i < IMQ_VAR_NUM; i++) {
                        for (j = 0; j < IMQ_VAR_NUM; j++) {
                            if (lin_dep[multix_1][i] && lin_dep[multix_2][j])
                                reduced_system[eq_idx][deg2midx2(i, j)] ^= 1;
                        }
                    }

                    if (lin_dep[multix_1][IMQ_VAR_NUM]) {
                        for (i = 0; i < IMQ_VAR_NUM; i++) {
                            if (lin_dep[multix_2][i])
                                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, i)] ^= 1;
                        }
                    }

                    if (lin_dep[multix_2][IMQ_VAR_NUM]) {
                        for (i = 0; i < IMQ_VAR_NUM; i++) {
                            if (lin_dep[multix_1][i])
                                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, i)] ^= 1;
                        }
                    }

                    if (lin_dep[multix_1][IMQ_VAR_NUM] && lin_dep[multix_2][IMQ_VAR_NUM])
                        reduced_system[eq_idx][IMQ_XVAR_NUM - 1] ^= 1;

                } else if (is_multix1_dep) {
                    uint32_t tmp_idx = lin2mq[multix_2];
                    for (i = 0; i < IMQ_VAR_NUM; i++) {
                        if (lin_dep[multix_1][i]) {
                            reduced_system[eq_idx][deg2midx2(tmp_idx, i)] ^= 1;
                        }
                    }

                    if (lin_dep[multix_1][IMQ_VAR_NUM])
                        reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, tmp_idx)] ^= 1;

                } else if (is_multix2_dep) {
                    uint32_t tmp_idx = lin2mq[multix_1];
                    for (i = 0; i < IMQ_VAR_NUM; i++) {
                        if (lin_dep[multix_2][i]) {
                            reduced_system[eq_idx][deg2midx2(tmp_idx, i)] ^= 1;
                        }
                    }

                    if (lin_dep[multix_2][IMQ_VAR_NUM])
                        reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, tmp_idx)] ^= 1;

                } else {
                    uint32_t tmp_idx1 = lin2mq[multix_1];
                    uint32_t tmp_idx2 = lin2mq[multix_2];
                    reduced_system[eq_idx][deg2midx2(tmp_idx1, tmp_idx2)] ^= 1;
                }
            }
        }
        if (mqsystem[eq_idx][deg2midx1(800, multix_1)]) {
            if (is_multix1_dep) {
                for (i = 0; i < IMQ_VAR_NUM + 1; i++) {
                    if (lin_dep[multix_1][i])
                        reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, i)] ^= 1;
                }
            } else {
                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, lin2mq[multix_1])] ^= 1;
            }
        }
    }
    reduced_system[eq_idx][IMQ_XVAR_NUM - 1] ^= mqsystem[eq_idx][BMQ_XVAR_NUM - 1];
}

void
reduceIterativeConstraints(MathSystem *system, uint8_t iterative_constr[LIN_ITER_EQNUM][801]) {
    uint32_t i, j, k;
    PRINTF_STAMP("\t\treducing round 3 iterative constraints...\n");
    for (i = 0; i < LIN_ITER_EQNUM; i++) {
        for (j = 0; j < 800; j++) {
            if (iterative_constr[i][j]) {
                if (system->round3_lin2mq[j] == DEP_PLACEMENT) {
                    for (k = 0; k < IMQ_VAR_NUM + 1; k++)
                        system->round3_iter_system[i][k] ^= system->round3_lin_dep[j][k];
                } else {
                    system->round3_iter_system[i][system->round3_lin2mq[j]] ^= 1;
                }
            }
        }
        system->round3_iter_system[i][IMQ_VAR_NUM] ^= iterative_constr[i][800];
    }
}

void
guessingBitsToMqSystem(const MathSystem *system,
                       const uint64_t guessingBits,
                       uint8_t *mqbuffer,
                       uint32_t *mq2lin,
                       uint32_t *lin2mq,
                       uint8_t *lin_dep) {
    uint32_t i = 0;
    uint32_t j = 0;
    uint32_t k = 0;

    /* create a temporary copy of reduced iterative system */
    uint8_t it_tmp[LIN_ITER_EQNUM][IMQ_VAR_NUM + 1];
    for (i = 0; i < LIN_ITER_EQNUM; i++)
        memcpy(it_tmp[i], system->round3_iter_system[i], IMQ_VAR_NUM + 1);

    /* assign the rhs coef according to guessing bits */
    for (i = 0; i < LIN_ITER_EQNUM; i++) {
        uint8_t c = (uint8_t) ((guessingBits >> i) & 0x1U);
        it_tmp[i][IMQ_VAR_NUM] ^= c;
    }

    /* perform gaussian elimination on iterative system */
    uint8_t lin_pivot[IMQ_VAR_NUM] = {0};
    uint8_t tmprow[IMQ_VAR_NUM + 1] = {0};
    for (i = 0; i < LIN_ITER_EQNUM; i++) {
        for (j = 0; j < IMQ_VAR_NUM; j++) {
            for (k = i; k < LIN_ITER_EQNUM; k++) {
                if (it_tmp[k][j])
                    break;
            }
            if (k != LIN_ITER_EQNUM)
                break;
        }
        if (j == IMQ_VAR_NUM)
            continue;

        memcpy(tmprow, it_tmp[k], IMQ_VAR_NUM + 1);
        memcpy(it_tmp[k], it_tmp[i], IMQ_VAR_NUM + 1);
        memcpy(it_tmp[i], tmprow, IMQ_VAR_NUM + 1);

        lin_pivot[j] = 1;

        for (k = 0; k < LIN_ITER_EQNUM; k++) {
            if (it_tmp[k][j] && k != i)
                for (uint32_t v = j; v < IMQ_VAR_NUM + 1; v++)
                    it_tmp[k][v] ^= it_tmp[i][v];
        }
    }

    /* record free variables after adding iterative constraints */
    uint32_t fvar_offset = 0;
    // mq2lin
    uint32_t fvar_index[AMQ_VAR_NUM] = {0};
    // lin2mq
    uint32_t rfvar_index[IMQ_VAR_NUM] = {0};
    for (i = 0; i < IMQ_VAR_NUM; i++) {
        if (lin_pivot[i] == 0) {
            fvar_index[fvar_offset] = i;
            rfvar_index[i] = fvar_offset;
            fvar_offset++;
        } else {
            rfvar_index[i] = DEP_PLACEMENT;
        }
    }

    /* calculating linear dependency */
    uint8_t tmp_lindep[IMQ_VAR_NUM][AMQ_VAR_NUM + 1];
    for (i = 0; i < IMQ_VAR_NUM; i++)
        memset(tmp_lindep[i], 0, (AMQ_VAR_NUM + 1) * sizeof(uint8_t));

    for (i = 0; i < LIN_ITER_EQNUM; i++) {
        for (j = 0; j < IMQ_VAR_NUM; j++) {
            if (it_tmp[i][j] == 1) {
                for (k = j + 1; k < IMQ_VAR_NUM; k++) {
                    if (it_tmp[i][k] == 1)
                        tmp_lindep[j][rfvar_index[k]] = 1;
                }
                tmp_lindep[j][AMQ_VAR_NUM] = it_tmp[i][IMQ_VAR_NUM];
                break;
            }
        }
    }

    /* use new variables to eliminate append system */
    uint8_t append_lin_system[AMQ_LIN_EQNUM][AMQ_XVAR_NUM];
    uint32_t eq_idx, var_idx;
    for (eq_idx = 0; eq_idx < AMQ_LIN_EQNUM; eq_idx++) {
        memset(append_lin_system[eq_idx], 0, AMQ_XVAR_NUM * sizeof(uint8_t));
        uint32_t multix_1;
        uint32_t multix_2;
        for (multix_1 = 0; multix_1 < IMQ_VAR_NUM; multix_1++) {
            bool is_multix1_dep = (bool) (rfvar_index[multix_1] == DEP_PLACEMENT);

            for (multix_2 = multix_1; multix_2 < IMQ_VAR_NUM; multix_2++) {
                var_idx = deg2midx2(multix_1, multix_2);
                if (system->round3_append_system[eq_idx][var_idx]) {
                    bool is_multix2_dep = (bool) (rfvar_index[multix_2] == DEP_PLACEMENT);

                    if (is_multix1_dep && is_multix2_dep) {
                        for (i = 0; i < AMQ_VAR_NUM; i++) {
                            for (j = 0; j < AMQ_VAR_NUM; j++) {
                                append_lin_system[eq_idx][deg2midx2(i, j)] ^=
                                    (tmp_lindep[multix_1][i] && tmp_lindep[multix_2][j]);
                            }
                        }

                        if (tmp_lindep[multix_1][AMQ_VAR_NUM]) {
                            for (i = 0; i < AMQ_VAR_NUM; i++) {
                                append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, i)] ^= tmp_lindep[multix_2][i];
                            }
                        }

                        if (tmp_lindep[multix_2][AMQ_VAR_NUM]) {
                            for (i = 0; i < AMQ_VAR_NUM; i++) {
                                append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, i)] ^= tmp_lindep[multix_1][i];
                            }
                        }

                        append_lin_system[eq_idx][AMQ_XVAR_NUM - 1] ^=
                            (tmp_lindep[multix_1][AMQ_VAR_NUM] && tmp_lindep[multix_2][AMQ_VAR_NUM]);

                    } else if (is_multix1_dep) {
                        uint32_t tmp_idx = rfvar_index[multix_2];
                        for (i = 0; i < AMQ_VAR_NUM; i++) {
                            append_lin_system[eq_idx][deg2midx2(tmp_idx, i)] ^= tmp_lindep[multix_1][i];
                        }

                        append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, tmp_idx)] ^= tmp_lindep[multix_1][AMQ_VAR_NUM];

                    } else if (is_multix2_dep) {
                        uint32_t tmp_idx = rfvar_index[multix_1];
                        for (i = 0; i < AMQ_VAR_NUM; i++) {
                            append_lin_system[eq_idx][deg2midx2(tmp_idx, i)] ^= tmp_lindep[multix_2][i];
                        }

                        append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, tmp_idx)] ^= tmp_lindep[multix_2][AMQ_VAR_NUM];
                    } else {
                        uint32_t tmp_idx1 = rfvar_index[multix_1];
                        uint32_t tmp_idx2 = rfvar_index[multix_2];
                        append_lin_system[eq_idx][deg2midx2(tmp_idx1, tmp_idx2)] ^= 1;
                    }
                }
            }
            if (system->round3_append_system[eq_idx][deg2midx1(IMQ_VAR_NUM, multix_1)]) {
                if (is_multix1_dep) {
                    for (i = 0; i < AMQ_VAR_NUM + 1; i++) {
                        append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, i)] ^= tmp_lindep[multix_1][i];
                    }
                } else {
                    append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, rfvar_index[multix_1])] ^= 1;
                }
            }
        }
        append_lin_system[eq_idx][AMQ_XVAR_NUM - 1] ^= system->round3_append_system[eq_idx][IMQ_XVAR_NUM - 1];
        for (i = 0; i < AMQ_VAR_NUM; i++) {
            for (j = i; j < AMQ_VAR_NUM; j++)
                assert(append_lin_system[eq_idx][deg2midx2(i, j)] == 0);
        }
    }


    /* use linearized append system to reduce 10 variables in mq system */
    uint8_t append_lin_part[AMQ_LIN_EQNUM][AMQ_VAR_NUM + 1];
    for (eq_idx = 0; eq_idx < AMQ_LIN_EQNUM; eq_idx++)
        memcpy(append_lin_part[eq_idx], &append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, 0)], AMQ_VAR_NUM + 1);

    uint8_t append_lin_pivot[AMQ_VAR_NUM] = {0};
    uint8_t append_tmprow[AMQ_VAR_NUM + 1] = {0};
    for (i = 0; i < AMQ_LIN_EQNUM; i++) {
        for (j = 0; j < AMQ_VAR_NUM; j++) {
            for (k = i; k < AMQ_LIN_EQNUM; k++) {
                if (append_lin_part[k][j])
                    break;
            }
            if (k != AMQ_LIN_EQNUM)
                break;
        }
        if (j == AMQ_VAR_NUM)
            continue;

        memcpy(append_tmprow, append_lin_part[k], AMQ_VAR_NUM + 1);
        memcpy(append_lin_part[k], append_lin_part[i], AMQ_VAR_NUM + 1);
        memcpy(append_lin_part[i], append_tmprow, AMQ_VAR_NUM + 1);

        append_lin_pivot[j] = 1;

        for (k = 0; k < AMQ_LIN_EQNUM; k++) {
            if (append_lin_part[k][j] && k != i)
                for (uint32_t v = j; v < AMQ_VAR_NUM + 1; v++)
                    append_lin_part[k][v] ^= append_lin_part[i][v];
        }
    }

    /* back substitution and update linear dependency */
    uint32_t mqfvar_offset = 0;
    // mq to amq
    uint32_t mqfvar_idx[MQ_VAR_NUM] = {0};
    // amq to mq
    uint32_t mqrfvar_idx[AMQ_VAR_NUM] = {0};
    for (i = 0; i < AMQ_VAR_NUM; i++) {
        if (append_lin_pivot[i] == 0) {
            mqfvar_idx[mqfvar_offset] = i;
            mqrfvar_idx[i] = mqfvar_offset;
            mqfvar_offset++;
        } else {
            mqrfvar_idx[i] = DEP_PLACEMENT;
        }
    }

    uint8_t append_lindep[AMQ_VAR_NUM][MQ_VAR_NUM + 1];
    for (i = 0; i < AMQ_VAR_NUM; i++)
        memset(append_lindep[i], 0, (MQ_VAR_NUM + 1) * sizeof(uint8_t));

    for (i = 0; i < AMQ_LIN_EQNUM; i++) {
        for (j = 0; j < AMQ_VAR_NUM; j++) {
            if (append_lin_part[i][j] == 1) {
                for (k = j + 1; k < AMQ_VAR_NUM; k++) {
                    if (append_lin_part[i][k] == 1) {
                        append_lindep[j][mqrfvar_idx[k]] = 1;
                    }
                }
                append_lindep[j][MQ_VAR_NUM] = append_lin_part[i][AMQ_VAR_NUM];
                break;
            }
        }
    }

    uint8_t total_lindep[800][MQ_VAR_NUM + 1];
    for (i = 0; i < 800; i++)
        memset(total_lindep[i], 0, (MQ_VAR_NUM + 1) * sizeof(uint8_t));

    // mq to lin
    uint32_t total_fvar_idx[MQ_VAR_NUM] = {0};
    for (i = 0; i < MQ_VAR_NUM; i++) {
        total_fvar_idx[i] = system->round3_mq2lin[fvar_index[mqfvar_idx[i]]];
    }

    // lin to mq
    uint32_t total_rfvar_idx[800] = {0};
    // reduction path: 800 (org) -> 94 (it) -> 41 (ap) -> 31 (mq)
    // 1. release 10 variables in (ap) to (org), s.t. 10 variables in
    // (org) can be represented by 31 variables
    for (i = 0; i < AMQ_VAR_NUM; i++) {
        if (mqrfvar_idx[i] == DEP_PLACEMENT) {
            uint32_t org_idx = system->round3_mq2lin[fvar_index[i]];
            total_rfvar_idx[org_idx] = DEP_PLACEMENT;
            for (j = 0; j < MQ_VAR_NUM + 1; j++) {
                total_lindep[org_idx][j] = append_lindep[i][j];
            }
        }
    }

    // 2. release 10 variables in (ap) to (it), s.t. all 94 variables
    // in it can be represented by 31 variables
    // two cases:
    //      i. variable in MQ system
    //      ii. variable in append system
    for (i = 0; i < IMQ_VAR_NUM; i++) {
        if (rfvar_index[i] == DEP_PLACEMENT) {
            // variable not in MQ system and not in append system
            uint32_t org_idx = system->round3_mq2lin[i];
            total_rfvar_idx[org_idx] = DEP_PLACEMENT;
            for (j = 0; j < AMQ_VAR_NUM; j++) {
                if (tmp_lindep[i][j]) {
                    if (mqrfvar_idx[j] == DEP_PLACEMENT) {
                        // variable in append system
                        uint32_t append_idx = system->round3_mq2lin[fvar_index[j]];
                        for (k = 0; k < MQ_VAR_NUM + 1; k++)
                            total_lindep[org_idx][k] ^= total_lindep[append_idx][k];
                    } else {
                        // variable in MQ system
                        total_lindep[org_idx][mqrfvar_idx[j]] ^= 1;
                    }
                }
            }
            total_lindep[org_idx][MQ_VAR_NUM] ^= tmp_lindep[i][AMQ_VAR_NUM];
        }
    }

    // 3. release 707 variables in (org), s.t. all 707 variables can be
    // represented by 31 variables
    // three cases:
    //      i. variable in MQ system (resolved)
    //      ii. variable in (ap) system (resolved)
    //      iii. variable in (it) system (resolved)
    for (i = 0; i < 800; i++) {
        if (system->round3_lin2mq[i] == DEP_PLACEMENT) {
            // variable in 707 dependent variables
            total_rfvar_idx[i] = DEP_PLACEMENT;
            for (j = 0; j < IMQ_VAR_NUM; j++) {
                if (system->round3_lin_dep[i][j]) {
                    if (rfvar_index[j] == DEP_PLACEMENT) {
                        // variable in (it) system or in (ap) system
                        uint32_t append_idx = system->round3_mq2lin[j];
                        for (k = 0; k < MQ_VAR_NUM + 1; k++)
                            total_lindep[i][k] ^= total_lindep[append_idx][k];
                    } else if (mqrfvar_idx[rfvar_index[j]] == DEP_PLACEMENT) {
                        uint32_t append_idx = system->round3_mq2lin[j];
                        for (k = 0; k < MQ_VAR_NUM + 1; k++)
                            total_lindep[i][k] ^= total_lindep[append_idx][k];
                    } else {
                        // variable in MQ system
                        total_lindep[i][mqrfvar_idx[rfvar_index[j]]] ^= 1;
                    }
                }
            }
            total_lindep[i][MQ_VAR_NUM] ^= system->round3_lin_dep[i][IMQ_VAR_NUM];
        }
    }

    /* reduce round3 mq system */
    uint8_t mqsystem[MQ_EQ_NUM][MQ_XVAR_NUM];
    for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
        memset(mqsystem[eq_idx], 0, MQ_XVAR_NUM);
        uint32_t multix_1;
        uint32_t multix_2;
        for (multix_1 = 0; multix_1 < IMQ_VAR_NUM; multix_1++) {
            bool is_multix1_dep;
            if (rfvar_index[multix_1] == DEP_PLACEMENT)
                is_multix1_dep = true;
            else if (mqrfvar_idx[rfvar_index[multix_1]] == DEP_PLACEMENT)
                is_multix1_dep = true;
            else
                is_multix1_dep = false;

            for (multix_2 = multix_1; multix_2 < IMQ_VAR_NUM; multix_2++) {
                var_idx = deg2midx2(multix_1, multix_2);
                if (system->round3_mq_system[eq_idx][var_idx]) {
                    bool is_multix2_dep;
                    if (rfvar_index[multix_2] == DEP_PLACEMENT)
                        is_multix2_dep = true;
                    else if (mqrfvar_idx[rfvar_index[multix_2]] == DEP_PLACEMENT)
                        is_multix2_dep = true;
                    else
                        is_multix2_dep = false;

                    if (is_multix1_dep && is_multix2_dep) {
                        uint32_t tmpidx_1 = system->round3_mq2lin[multix_1];
                        uint32_t tmpidx_2 = system->round3_mq2lin[multix_2];
                        for (i = 0; i < MQ_VAR_NUM; i++) {
                            for (j = 0; j < MQ_VAR_NUM; j++) {
                                mqsystem[eq_idx][deg2midx2(i, j)] ^=
                                    (total_lindep[tmpidx_1][i] & total_lindep[tmpidx_2][j]);
                            }
                        }

                        if (total_lindep[tmpidx_1][MQ_VAR_NUM]) {
                            for (i = 0; i < MQ_VAR_NUM; i++) {
                                mqsystem[eq_idx][deg2midx1(MQ_VAR_NUM, i)] ^= total_lindep[tmpidx_2][i];
                            }
                        }

                        if (total_lindep[tmpidx_2][MQ_VAR_NUM]) {
                            for (i = 0; i < MQ_VAR_NUM; i++) {
                                mqsystem[eq_idx][deg2midx1(MQ_VAR_NUM, i)] ^= total_lindep[tmpidx_1][i];
                            }
                        }

                        if (total_lindep[tmpidx_1][MQ_VAR_NUM] && total_lindep[tmpidx_2][MQ_VAR_NUM])
                            mqsystem[eq_idx][MQ_XVAR_NUM - 1] ^= 1;

                    } else if (is_multix1_dep) {
                        uint32_t tmpidx_1 = system->round3_mq2lin[multix_1];
                        uint32_t tmpidx_2 = mqrfvar_idx[rfvar_index[multix_2]];
                        for (i = 0; i < MQ_VAR_NUM; i++) {
                            mqsystem[eq_idx][deg2midx2(tmpidx_2, i)] ^= total_lindep[tmpidx_1][i];
                        }

                        mqsystem[eq_idx][deg2midx1(MQ_VAR_NUM, tmpidx_2)] ^= total_lindep[tmpidx_1][MQ_VAR_NUM];

                    } else if (is_multix2_dep) {
                        uint32_t tmpidx_1 = mqrfvar_idx[rfvar_index[multix_1]];
                        uint32_t tmpidx_2 = system->round3_mq2lin[multix_2];
                        for (i = 0; i < MQ_VAR_NUM; i++) {
                            mqsystem[eq_idx][deg2midx2(tmpidx_1, i)] ^= total_lindep[tmpidx_2][i];
                        }

                        mqsystem[eq_idx][deg2midx1(MQ_VAR_NUM, tmpidx_1)] ^= total_lindep[tmpidx_2][MQ_VAR_NUM];

                    } else {
                        uint32_t tmpidx_1 = mqrfvar_idx[rfvar_index[multix_1]];
                        uint32_t tmpidx_2 = mqrfvar_idx[rfvar_index[multix_2]];
                        mqsystem[eq_idx][deg2midx2(tmpidx_1, tmpidx_2)] ^= 1;
                    }
                }
            }
            if (system->round3_mq_system[eq_idx][deg2midx1(IMQ_VAR_NUM, multix_1)]) {
                if (is_multix1_dep) {
                    for (i = 0; i < MQ_VAR_NUM + 1; i++) {
                        mqsystem[eq_idx][deg2midx1(MQ_VAR_NUM, i)] ^= total_lindep[system->round3_mq2lin[multix_1]][i];
                    }
                } else {
                    mqsystem[eq_idx][deg2midx1(MQ_VAR_NUM, mqrfvar_idx[rfvar_index[multix_1]])] ^= 1;
                }
            }
        }
        mqsystem[eq_idx][MQ_XVAR_NUM - 1] ^= system->round3_mq_system[eq_idx][IMQ_XVAR_NUM - 1];
    }

#ifdef _VERIFICATION
    // validate append system reduction result
    // assume all variable in AMQ_System are set to 1
    uint8_t append_eval[IMQ_VAR_NUM] = {0};
    for (i = 0; i < IMQ_VAR_NUM; i++) {
        if (rfvar_index[i] == DEP_PLACEMENT) {
            for (j = 0; j < AMQ_VAR_NUM + 1; j++)
                append_eval[i] ^= tmp_lindep[i][j];
        } else {
            append_eval[i] = 1;
        }
    }
    for (eq_idx = 0; eq_idx < AMQ_LIN_EQNUM; eq_idx++) {
        uint32_t verifyResult1 = 0;
        uint32_t verifyResult2 = 0;
        for (i = 0; i < AMQ_VAR_NUM + 1; i++) {
            verifyResult1 ^= append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, i)];
        }
        for (i = 0; i < IMQ_VAR_NUM; i++) {
            for (j = i; j < IMQ_VAR_NUM; j++) {
                if (system->round3_append_system[eq_idx][deg2midx2(i, j)])
                    verifyResult2 ^= (append_eval[i] & append_eval[j]);
            }
            if (system->round3_append_system[eq_idx][deg2midx1(IMQ_VAR_NUM, i)])
                verifyResult2 ^= append_eval[i];
        }
        if (system->round3_append_system[eq_idx][IMQ_XVAR_NUM - 1])
            verifyResult2 ^= 1;

        if (verifyResult1 != verifyResult2) {
            EXIT_WITH_MSG("reduce round4 append system is not equal with the original system, "
                          "please check the implementation.\n");
        }
    }

    uint8_t mq_eval[IMQ_VAR_NUM] = {0};
    for (i = 0; i < IMQ_VAR_NUM; i++) {
        if (total_rfvar_idx[system->round3_mq2lin[i]] == DEP_PLACEMENT) {
            for (j = 0; j < MQ_VAR_NUM + 1; j++)
                mq_eval[i] ^= total_lindep[system->round3_mq2lin[i]][j];
        } else {
            mq_eval[i] = 1;
        }
    }
    for (eq_idx = 0; eq_idx < MQ_EQ_NUM; eq_idx++) {
        uint32_t verifyResult1 = 0;
        uint32_t verifyResult2 = 0;
        for (i = 0; i < MQ_VAR_NUM; i++) {
            for (j = i; j < MQ_VAR_NUM; j++) {
                if (mqsystem[eq_idx][deg2midx2(i, j)])
                    verifyResult1 ^= 1;
            }
            if (mqsystem[eq_idx][deg2midx1(MQ_VAR_NUM, i)])
                verifyResult1 ^= 1;
        }
        if (mqsystem[eq_idx][MQ_XVAR_NUM - 1])
            verifyResult1 ^= 1;

        for (i = 0; i < IMQ_VAR_NUM; i++) {
            for (j = i; j < IMQ_VAR_NUM; j++) {
                if (system->round3_mq_system[eq_idx][deg2midx2(i, j)])
                    verifyResult2 ^= (mq_eval[i] & mq_eval[j]);
            }
            if (system->round3_mq_system[eq_idx][deg2midx1(IMQ_VAR_NUM, i)])
                verifyResult2 ^= mq_eval[i];
        }
        if (system->round3_mq_system[eq_idx][IMQ_XVAR_NUM - 1])
            verifyResult2 ^= 1;

        if (verifyResult1 != verifyResult2) {
            EXIT_WITH_MSG("reduced round4 mq system is not equal with the original system, "
                          "please check the implementation.\n");
        }
    }
#endif

    /* copy mq system to memory */
    for (i = 0; i < MQ_VAR_NUM; i++) {
        total_rfvar_idx[total_fvar_idx[i]] = i;
    }
    for (i = 0; i < 800; i++)
        memcpy(lin_dep + i * (MQ_VAR_NUM + 1), total_lindep[i], (MQ_VAR_NUM + 1) * sizeof(uint8_t));
    for (i = 0; i < MQ_EQ_NUM; i++)
        memcpy(mqbuffer + i * (MQ_XVAR_NUM), mqsystem[i], (MQ_XVAR_NUM) * sizeof(uint8_t));
    memcpy(lin2mq, total_rfvar_idx, 800 * sizeof(uint32_t));
    memcpy(mq2lin, total_fvar_idx, MQ_VAR_NUM * sizeof(uint32_t));
}

void
freeMathSystem(MathSystem *system) {
    uint32_t i;
    for (i = 0; i < 800; i++)
        SFREE(system->round3_lin_dep[i]);
    for (i = 0; i < AMQ_LIN_EQNUM; i++)
        SFREE(system->round3_append_system[i]);
    for (i = 0; i < MQ_EQ_NUM; i++)
        SFREE(system->round3_mq_system[i]);
    for (i = 0; i < LIN_ITER_EQNUM; i++)
        SFREE(system->round3_iter_system[i]);

    SFREE(system->round3_mq2lin);
    SFREE(system->round3_lin2mq);
    SFREE(system->round4_lin_dep);
}
