//
// Created by vitowu on 3/18/20.
//

#include "math.h"

#include "params.h"
#include "utils.h"

#define eqvar(i, j, size) \
        ((i) * (size) + (j))

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

    PRINTF_STAMP("reducing append system:\tequation %d...\n", eq_idx);
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
                        if (lin_dep[multix_1][i]) {
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

    uint8_t **reduced_system = system->round3_append_system;
    uint8_t **lin_dep = system->round3_lin_dep;
    uint32_t *lin2mq = system->round3_lin2mq;

    PRINTF_STAMP("reducing mq system:\tequation %d...\n", eq_idx);
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
                        if (lin_dep[multix_1][i]) {
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
    PRINTF_STAMP("reducing round 3 iterative constraints...\n");
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
        uint8_t c = (uint8_t) ((guessingBits >> i) & 1);
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
    uint32_t fvar_index[AMQ_VAR_NUM];
    // lin2mq
    uint32_t rfvar_index[IMQ_VAR_NUM];
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
        memset(append_lin_system[eq_idx], 0, AMQ_XVAR_NUM);
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
                                if (tmp_lindep[multix_1][i] && tmp_lindep[multix_2][j])
                                    append_lin_system[eq_idx][deg2midx2(i, j)] ^= 1;
                            }
                        }

                        if (tmp_lindep[multix_1][AMQ_VAR_NUM]) {
                            for (i = 0; i < AMQ_VAR_NUM; i++) {
                                if (tmp_lindep[multix_2][i])
                                    append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, i)] ^= 1;
                            }
                        }

                        if (tmp_lindep[multix_2][AMQ_VAR_NUM]) {
                            for (i = 0; i < AMQ_VAR_NUM; i++) {
                                if (tmp_lindep[multix_1][i])
                                    append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, i)] ^= 1;
                            }
                        }

                        if (tmp_lindep[multix_1][AMQ_VAR_NUM] && tmp_lindep[multix_2][AMQ_VAR_NUM])
                            append_lin_system[eq_idx][AMQ_XVAR_NUM - 1] ^= 1;

                    } else if (is_multix1_dep) {
                        uint32_t tmp_idx = rfvar_index[multix_2];
                        for (i = 0; i < AMQ_VAR_NUM; i++) {
                            if (tmp_lindep[multix_1][i]) {
                                append_lin_system[eq_idx][deg2midx2(tmp_idx, i)] ^= 1;
                            }
                        }

                        if (tmp_lindep[multix_1][AMQ_VAR_NUM])
                            append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, tmp_idx)] ^= 1;

                    } else if (is_multix2_dep) {
                        uint32_t tmp_idx = rfvar_index[multix_1];
                        for (i = 0; i < AMQ_VAR_NUM; i++) {
                            if (tmp_lindep[multix_2][i]) {
                                append_lin_system[eq_idx][deg2midx2(tmp_idx, i)] ^= 1;
                            }
                        }

                        if (tmp_lindep[multix_2][AMQ_VAR_NUM])
                            append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, tmp_idx)] ^= 1;
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
                        if (tmp_lindep[multix_1][i])
                            append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, i)] ^= 1;
                    }
                } else {
                    append_lin_system[eq_idx][deg2midx1(AMQ_VAR_NUM, rfvar_index[multix_1])] ^= 1;
                }
            }
        }
        append_lin_system[eq_idx][AMQ_XVAR_NUM - 1] ^= system->round3_append_system[eq_idx][IMQ_XVAR_NUM - 1];
        for (i = 0; i < AMQ_VAR_NUM; i++)
            for (j = i; j < AMQ_VAR_NUM; j++)
                if (append_lin_system[eq_idx][deg2midx2(i, j)])
                    printf("%d %d\n", i, j);
        exit(0);
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

    /* back substitution */
    uint32_t mq_fvar_offset = 0;

    /* update linear dependency */


    /* reduce round3 mq system */


    /* copy mq system to memory */

    exit(0);
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
