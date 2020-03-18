//
// Created by vitowu on 3/18/20.
//

#include "math.h"

#include "params.h"

#define DEP_PLACEMENT 0xffffffff

#define eqvar(i, j, size) \
        ((i) * (size) + (j))

void
initMathSystem(MathSystem *system) {
    uint32_t i;

    // linear dependency: x_k = x_n_1 + x_n_2 + ... + x_n_j
    for (i = 0; i < 800; i++)
        system->round3_lin_dep[i] = (uint32_t *) malloc((IMQ_VAR_NUM + 1) * sizeof(uint32_t));

    // append system in order 2
    for (i = 0; i < AMQ_LIN_EQNUM; i++)
        system->round3_append_system[i] = (uint8_t *) malloc(IMQ_XVAR_NUM * (sizeof(uint8_t)));

    // mq system in order 2
    for (i = 0; i < MQ_EQ_NUM; i++)
        system->round3_mq_system[i] = (uint8_t *) malloc(IMQ_XVAR_NUM * sizeof(uint8_t));

    // transfer mq var index (94) to linear var index (94)
    system->round3_mq2lin = (uint32_t *) malloc(IMQ_VAR_NUM * sizeof(uint32_t));

    // reversely transfer linear index to mq var index
    system->round3_lin2mq = (uint32_t *) malloc(800 * sizeof(uint32_t));

    // linear dependency
    system->round4_lin_dep = (uint32_t *) malloc(IMQ_VAR_NUM * (AMQ_VAR_NUM + 1) * sizeof(uint32_t));
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

    for (uint32_t i = 0; i < LIN_CONST_EQNUM; i++) {
        for (uint32_t j = 0; j < 800; j++) {
            if (lin_system[i][j]) {
                // memcpy(lin_dep[j], it_tmp[i], 801);
                // lin_dep[j][j] = 0;
                memcpy(system->round3_lin_dep[j], lin_system[i], 801);
                system->round3_lin_dep[j][j] = 0;
                break;
            }
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
}

void
reduceRound3AppendSystem(MathSystem *system, uint8_t **append_system) {
    memset(system->round3_append_system, 0, IMQ_XVAR_NUM * AMQ_LIN_EQNUM);
    uint32_t eq_idx, var_idx, i, j;

    uint8_t **reduced_system = system->round3_append_system;
    uint32_t **lin_dep = system->round3_lin_dep;
    uint32_t *lin2mq = system->round3_lin2mq;

    for (eq_idx = 0; eq_idx < AMQ_LIN_EQNUM; eq_idx++) {
        uint32_t multix_1;
        uint32_t multix_2;
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
                        for (i = 0; i < 800; i++) {
                            if (lin_dep[multix_1][i]) {
                                reduced_system[eq_idx][deg2midx2(tmp_idx, i)] ^= 1;
                            }

                            if (lin_dep[multix_1][IMQ_VAR_NUM])
                                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, tmp_idx)] ^= 1;
                        }

                    } else if (is_multix2_dep) {
                        uint32_t tmp_idx = system->round3_lin2mq[multix_1];
                        for (i = 0; i < 800; i++) {
                            if (lin_dep[multix_1][i]) {
                                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, tmp_idx)] ^= 1;
                            }

                            if (lin_dep[multix_2][IMQ_VAR_NUM])
                                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, tmp_idx)] ^= 1;
                        }

                    } else {
                        uint32_t tmp_idx1 = system->round3_lin2mq[multix_1];
                        uint32_t tmp_idx2 = system->round3_lin2mq[multix_2];
                        reduced_system[eq_idx][deg2midx2(tmp_idx1, tmp_idx2)] ^= 1;
                    }
                }
            }
        }
    }
}

void
reduceRound3MQSystem(MathSystem *system, uint8_t **mqsystem) {
    memset(system->round3_mq_system, 0, IMQ_XVAR_NUM * AMQ_LIN_EQNUM);
    uint32_t eq_idx, var_idx, i, j;

    uint8_t **reduced_system = system->round3_mq_system;
    uint32_t **lin_dep = system->round3_lin_dep;
    uint32_t *lin2mq = system->round3_lin2mq;

    for (eq_idx = 0; eq_idx < AMQ_LIN_EQNUM; eq_idx++) {
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
                        for (i = 0; i < 800; i++) {
                            if (lin_dep[multix_1][i]) {
                                reduced_system[eq_idx][deg2midx2(tmp_idx, i)] ^= 1;
                            }

                            if (lin_dep[multix_1][IMQ_VAR_NUM])
                                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, tmp_idx)] ^= 1;
                        }

                    } else if (is_multix2_dep) {
                        uint32_t tmp_idx = system->round3_lin2mq[multix_1];
                        for (i = 0; i < 800; i++) {
                            if (lin_dep[multix_1][i]) {
                                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, tmp_idx)] ^= 1;
                            }

                            if (lin_dep[multix_2][IMQ_VAR_NUM])
                                reduced_system[eq_idx][deg2midx1(IMQ_VAR_NUM, tmp_idx)] ^= 1;
                        }

                    } else {
                        uint32_t tmp_idx1 = system->round3_lin2mq[multix_1];
                        uint32_t tmp_idx2 = system->round3_lin2mq[multix_2];
                        reduced_system[eq_idx][deg2midx2(tmp_idx1, tmp_idx2)] ^= 1;
                    }
                }
            }
        }
    }
}

void
freeMathSystem(MathSystem *system) {
    uint32_t i;
    for (i = 0; i < 800; i++)
        free(system->round3_lin_dep[i]);
    for (i = 0; i < AMQ_LIN_EQNUM; i++)
        free(system->round3_append_system[i]);
    for (i = 0; i < MQ_EQ_NUM; i++)
        free(system->round3_mq_system);

    free(system->round3_mq2lin);
    free(system->round3_lin2mq);
    free(system->round4_lin_dep);
}
