//
// Created by vitowu on 3/11/20.
//

#include <iostream>
#include <cstring>

#define kidx(i, j, k) \
    ((((5 * i) + j) * 32) + k)

#define cbinom2(n) \
    ( ((n) * ((n)+1)) / 2)

#define deg2midx1(vnum, var1_idx) \
    (cbinom2(vnum) + (var1_idx))

// pre-requirements: var2_idx > var1_idx
#define deg2midx2(var1_idx, var2_idx) \
    (cbinom2(var2_idx) + (var1_idx))

#define BMQ_VAR_NUM 800
#define BMQ_EQ_NUM 38
#define BMQ_XVAR_NUM ((((BMQ_VAR_NUM) * ((BMQ_VAR_NUM) + 1)) / 2) + BMQ_VAR_NUM + 1)

#define BLIN_VAR_NUM 800
#define BLIN_EQ_NUM 707

#define RMQ_VAR_NUM 94
#define RMQ_EQ_NUM 38
#define RMQ_XVAR_NUM ((((RMQ_VAR_NUM) * ((RMQ_VAR_NUM) + 1)) / 2) + RMQ_VAR_NUM + 1)

#define IEQ_NUM 53

uint8_t mq_system_init[BMQ_EQ_NUM][BMQ_XVAR_NUM];
uint8_t mq_lin_system[BLIN_EQ_NUM][BLIN_VAR_NUM + 1];
uint8_t mq_lin_dep[BLIN_VAR_NUM][BLIN_VAR_NUM + 1];
uint32_t mq_fvar_idx[RMQ_VAR_NUM] = {0};
uint32_t reverse_mq_fvar_idx[BLIN_VAR_NUM];
uint8_t mq_system_result[RMQ_EQ_NUM][RMQ_XVAR_NUM];
uint8_t mq_iterative_system[IEQ_NUM][BLIN_VAR_NUM + 1];
uint8_t reduced_iterative_system[IEQ_NUM][RMQ_VAR_NUM + 1];

void reduceLinearSystem(uint8_t lin_system[BLIN_EQ_NUM][BLIN_VAR_NUM + 1],
                        uint32_t mq_fvar_idx[RMQ_VAR_NUM]) {
    uint8_t tmprow[BLIN_VAR_NUM];
    uint32_t mq_fvars[800] = {0};
    for (uint32_t i = 0; i < BLIN_EQ_NUM; i++) {
        uint32_t j = 0;
        uint32_t k = 0;
        for (j = i; j < BLIN_VAR_NUM; j++) {
            for (k = i; k < BLIN_EQ_NUM; k++) {
                if (lin_system[k][j])
                    break;
            }
            if (k != BLIN_EQ_NUM)
                break;
        }
        if (j == BLIN_VAR_NUM)
            continue;

        memcpy(tmprow, lin_system[k], BLIN_VAR_NUM + 1);
        memcpy(lin_system[k], lin_system[i], BLIN_VAR_NUM + 1);
        memcpy(lin_system[i], tmprow, BLIN_VAR_NUM + 1);

        // record linear dependency
        mq_fvars[j] = 1;
        memcpy(mq_lin_dep[j], lin_system[i], BLIN_VAR_NUM + 1);
        mq_lin_dep[j][j] = 0;

        for (k = 0; k < i; k++) {
            if (lin_system[k][j]) {
                for (uint32_t v = j; v < BLIN_VAR_NUM + 1; v++) {
                    lin_system[k][v] ^= lin_system[i][v];
                }
            }
        }
    }

    // record free variables
    uint32_t fvar_offset = 0;
    for (uint32_t i = 0; i < 800; i++) {
        if (mq_fvars[i] == 0) {
            mq_fvar_idx[fvar_offset] = i;
            reverse_mq_fvar_idx[i] = fvar_offset;
            fvar_offset++;
        } else {
            reverse_mq_fvar_idx[i] = 0xFFFFFFFF;
        }
    }
}

void reduceIterativeSystem() {
    uint32_t eq_idx, var_idx;
    for (eq_idx = 0; eq_idx < IEQ_NUM; eq_idx++) {
        for (var_idx = 0; var_idx < BLIN_VAR_NUM; var_idx++) {
            if (mq_iterative_system[eq_idx][var_idx]) {
                if (reverse_mq_fvar_idx[var_idx] == 0xFFFFFFFF) {
                    mq_iterative_system[eq_idx][var_idx] ^= 1;
                    for (uint32_t i = 0; i < BLIN_VAR_NUM + 1; i++) {
                        if (mq_lin_dep[var_idx][i])
                            mq_iterative_system[eq_idx][i] ^= 1;
                    }
                }
            }
        }
    }
    for (eq_idx = 0; eq_idx < IEQ_NUM; eq_idx++) {
        for (var_idx = 0; var_idx < RMQ_VAR_NUM; var_idx++) {
            if (mq_iterative_system[eq_idx][mq_fvar_idx[var_idx]]) {
                reduced_iterative_system[eq_idx][var_idx] ^= 1;
            }
        }
        if (mq_iterative_system[eq_idx][BLIN_VAR_NUM])
            reduced_iterative_system[eq_idx][RMQ_VAR_NUM] ^= 1;
    }
}

void reduceMQSystem() {
    uint32_t eq_idx;
    uint32_t var_idx;
    uint32_t i, j;
    for (eq_idx = 0; eq_idx < RMQ_EQ_NUM; eq_idx++) {
        memset(mq_system_result[eq_idx], 0, RMQ_XVAR_NUM * sizeof(uint8_t));
        uint32_t multix_1;
        uint32_t multix_2;
        for (multix_1 = 0; multix_1 < BLIN_VAR_NUM; multix_1++) {
            bool is_multix1_dep = (bool) (reverse_mq_fvar_idx[multix_1] == 0xFFFFFFFF);

            // modify quadratic part
            for (multix_2 = multix_1; multix_2 < BLIN_VAR_NUM; multix_2++) {
                var_idx = deg2midx2(multix_1, multix_2);
                if (mq_system_init[eq_idx][var_idx]) {
                    bool is_multix2_dep = (bool) (reverse_mq_fvar_idx[multix_2] == 0xFFFFFFFF);
                    if (is_multix1_dep && is_multix2_dep) {
                        // all multiplexers are dependent
                        for (i = 0; i < RMQ_VAR_NUM; i++)
                            for (j = 0; j < RMQ_VAR_NUM; j++)
                                if (mq_lin_dep[multix_1][i] && mq_lin_dep[multix_2][j])
                                    mq_system_result[eq_idx][deg2midx2(i, j)] ^= 1;

                        if (mq_lin_dep[multix_1][RMQ_VAR_NUM]) {
                            for (i = 0; i < RMQ_VAR_NUM; i++) {
                                if (mq_lin_dep[multix_2][i])
                                    mq_system_result[eq_idx][deg2midx1(RMQ_VAR_NUM, i)] ^= 1;
                            }
                        }

                        if (mq_lin_dep[multix_2][RMQ_VAR_NUM]) {
                            for (i = 0; i < RMQ_VAR_NUM; i++) {
                                if (mq_lin_dep[multix_1][i])
                                    mq_system_result[eq_idx][deg2midx1(RMQ_VAR_NUM, i)] ^= 1;
                            }
                        }

                        if (mq_lin_dep[multix_1][RMQ_VAR_NUM] && mq_lin_dep[multix_2][RMQ_VAR_NUM])
                            mq_system_result[eq_idx][RMQ_XVAR_NUM - 1] ^= 1;

                    } else if (is_multix1_dep) {
                        // only multiplexer 1 is dependent
                        uint32_t v = reverse_mq_fvar_idx[multix_2];
                        for (i = 0; i < RMQ_VAR_NUM; i++)
                            if (mq_lin_dep[multix_1][i])
                                mq_system_result[eq_idx][deg2midx2(i, v)] ^= 1;

                        if (mq_lin_dep[multix_1][RMQ_VAR_NUM])
                            mq_system_result[eq_idx][deg2midx1(RMQ_VAR_NUM, v)] ^= 1;
                    } else if (is_multix2_dep) {
                        // only multiplexer 2 is dependent
                        uint32_t v = reverse_mq_fvar_idx[multix_1];
                        for (i = 0; i < RMQ_VAR_NUM; i++)
                            if (mq_lin_dep[multix_2][i])
                                mq_system_result[eq_idx][deg2midx2(i, v)] ^= 1;

                        if (mq_lin_dep[multix_2][RMQ_VAR_NUM])
                            mq_system_result[eq_idx][deg2midx1(RMQ_VAR_NUM, v)] ^= 1;
                    } else {
                        // all multiplexers are free variables
                        uint32_t v1 = reverse_mq_fvar_idx[multix_1];
                        uint32_t v2 = reverse_mq_fvar_idx[multix_2];
                        if (v1 >= v2)
                            mq_system_result[eq_idx][deg2midx2(v2, v1)] ^= 1;
                        else
                            mq_system_result[eq_idx][deg2midx2(v1, v2)] ^= 1;
                    }
                }
            }
            var_idx = deg2midx1(BLIN_VAR_NUM, multix_1);
            if (mq_system_init[eq_idx][var_idx]) {
                if (is_multix1_dep) {
                    for (i = 0; i < RMQ_VAR_NUM; i++) {
                        if (mq_lin_dep[multix_1][i]) {
                            mq_system_result[eq_idx][deg2midx1(RMQ_VAR_NUM, i)] ^= 1;
                        }
                    }
                } else {
                    uint32_t v = reverse_mq_fvar_idx[multix_1];
                    mq_system_result[eq_idx][deg2midx1(RMQ_VAR_NUM, v)] ^= 1;
                }
            }
        }
        mq_system_result[eq_idx][RMQ_XVAR_NUM - 1] ^= mq_system_init[eq_idx][BMQ_XVAR_NUM - 1];
    }
}

int main() {
    FILE *fmq = fopen("../data/mq_analysis.dat", "r");
    char ch;
    uint32_t xv = 0;
    if (fmq == NULL) {
        printf("File is not available\n");
    } else {
        uint32_t k = 0;
        uint32_t eq_idx = 0;
        while ((ch = fgetc(fmq)) != EOF) {
            if (ch == '\n') {
                xv = k;
                k = 0;
                eq_idx++;
                continue;
            }
            mq_system_init[eq_idx][k] = (uint8_t) (ch - '0');
            k++;
        }
        fclose(fmq);
    }

//    printf("[+] read mq analysis from file\n"
//           "\tequation_number: %d\n"
//           "\txvar_number: %d\n"
//           "\tvar_number: %d\n",
//           BMQ_EQ_NUM, xv, BMQ_XVAR_NUM);

    FILE *flin = fopen("../data/lin_analysis.dat", "r");
    if (flin == NULL) {
        printf("File is not available\n");
    } else {
        for (uint32_t i = 0; i < BLIN_EQ_NUM; i++) {
            for (uint32_t j = 0; j < BLIN_VAR_NUM + 1; j++) {
                ch = fgetc(flin);
                mq_lin_system[i][j] = (uint8_t) (ch - '0');
            }
            fgetc(flin);
        }
        fclose(flin);
    }

    reduceLinearSystem(mq_lin_system, mq_fvar_idx);

//    for (uint32_t i = 0; i < BLIN_EQ_NUM; i++) {
//        for (uint32_t j = 0; j < BLIN_VAR_NUM + 1; j++)
//            printf("%d", mq_lin_system[i][j]);
//        printf("\n");
//    }

    FILE *fi = fopen("../data/variate_analysis.dat", "r");
    if (fi == NULL) {
        printf("File is not available\n");
    } else {
        for (uint32_t i = 0; i < IEQ_NUM; i++) {
            for (uint32_t j = 0; j < BLIN_VAR_NUM + 1; j++) {
                ch = fgetc(flin);
                mq_iterative_system[i][j] = (uint8_t) (ch - '0');
            }
            fgetc(fi);
        }
        fclose(fi);
    }
    reduceIterativeSystem();
//    for (uint32_t i = 0; i < 63; i++) {
//        for (uint32_t j = 0; j < BLIN_VAR_NUM + 1; j++)
//            printf("%d", mq_iterative_system[i][j]);
//        printf("\n");
//    }
//    for (uint32_t i = 0; i < IEQ_NUM; i++) {
//        for (uint32_t j = 0; j < RMQ_VAR_NUM + 1; j++)
//            printf("%d", reduced_iterative_system[i][j]);
//        printf("\n");
//    }
//    reduceMQSystem();

//    for (uint32_t i = 0; i < RMQ_VAR_NUM; i++)
//        printf("%d %d\n", i, mq_fvar_idx[i]);
//    for (uint32_t i = 0; i < RMQ_EQ_NUM; i++) {
//        for (uint32_t j = 0; j < RMQ_XVAR_NUM; j++)
//            printf("%d", mq_system_result[i][j]);
//        printf("\n");
//    }
//    for (uint32_t i = 0; i < BLIN_VAR_NUM; i++) {
//        bool isLinear = true;
//        for (uint32_t j = 0; j < RMQ_VAR_NUM; j++) {
//            if (mq_fvar_idx[j] == i) {
//                isLinear = false;
//                break;
//            }
//        }
//        if (!isLinear)
//            continue;
//        bool isConstance = true;
//        for (uint32_t j = 0; j < BLIN_VAR_NUM; j++) {
//            if (mq_lin_dep[i][j]) {
//                isConstance = false;
//                break;
//            }
//        }
//        if (isConstance)
//            printf("%d: %d\n", i, mq_lin_dep[i][BLIN_VAR_NUM]);
//    }
//    printf("[+] read constant linear analysis from file\n"
//           "\tequation_number: %d\n"
//           "\tvar_number: %d\n",
//           BLIN_EQ_NUM, BLIN_VAR_NUM);

    return 0;
}