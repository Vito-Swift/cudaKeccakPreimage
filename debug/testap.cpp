//
// Created by vitowu on 3/17/20.
//

#include <iostream>
#include <cstring>

#define AEQ_NUM 10
#define AVAR_NUM 800
#define AVAR_XVAR_NUM ((((AVAR_NUM) * ((AVAR_NUM) + 1)) / 2) + AVAR_NUM + 1)

#define LINEQ_NUM 707
#define ILINEQ_NUM 53
#define FVAR_NUM 41

#define kidx(i, j, k) \
    ((((5 * i) + j) * 32) + k)

#define cbinom2(n) \
    ( ((n) * ((n)+1)) / 2)

#define deg2midx1(vnum, var1_idx) \
    (cbinom2(vnum) + (var1_idx))

// pre-requirements: var2_idx > var1_idx
#define deg2midx2(var1_idx, var2_idx) \
    ((var2_idx) > (var1_idx) ? (cbinom2(var2_idx) + (var1_idx)) : (cbinom2(var1_idx) + (var2_idx)))

uint8_t append_system_origin[AEQ_NUM][AVAR_XVAR_NUM];
uint8_t lin_system[LINEQ_NUM][801];
uint8_t it_lin_system[ILINEQ_NUM][801];
uint8_t lin_dep[800][801];
uint32_t free_var_idx[FVAR_NUM];
uint32_t reverse_free_var_idx[800];
uint8_t reduced_append_system[AEQ_NUM][AVAR_XVAR_NUM];

void linearReduce(uint64_t guessingBits) {
    uint8_t it_tmp[LINEQ_NUM + ILINEQ_NUM][801];
    memcpy(it_tmp, lin_system, LINEQ_NUM * 801);
    memcpy(&it_tmp[LINEQ_NUM], it_lin_system, ILINEQ_NUM * 801);

    for (uint32_t i = 0; i < ILINEQ_NUM; i++) {
        uint8_t c = (guessingBits >> i) & 1;
        it_tmp[i + LINEQ_NUM][800] ^= c;
    }
    uint8_t lin_pivot[800] = {0};
    uint8_t tmprow[801];
    for (uint32_t i = 0; i < ILINEQ_NUM + LINEQ_NUM; i++) {
        uint32_t j = 0;
        uint32_t k = 0;
        for (j = 0; j < 800; j++) {
            for (k = i; k < ILINEQ_NUM + LINEQ_NUM; k++) {
                if (it_tmp[k][j])
                    break;
            }
            if (k != ILINEQ_NUM + LINEQ_NUM)
                break;
        }
        if (j == 800)
            continue;

        memcpy(tmprow, it_tmp[k], 801);
        memcpy(it_tmp[k], it_tmp[i], 801);
        memcpy(it_tmp[i], tmprow, 801);

        lin_pivot[j] = 1;
//        memcpy(lin_dep[j], it_tmp[i], 801);
//        lin_dep[j][j] = 0;

        for (k = 0; k < LINEQ_NUM + ILINEQ_NUM; k++) {
            if (it_tmp[k][j] && k != i)
                for (uint32_t v = j; v < 801; v++)
                    it_tmp[k][v] ^= it_tmp[i][v];
        }
    }

    for (uint32_t i = 0; i < LINEQ_NUM + ILINEQ_NUM; i++) {
        for (uint32_t j = 0; j < 800; j++) {
            if (it_tmp[i][j]) {
                memcpy(lin_dep[j], it_tmp[i], 801);
                lin_dep[j][j] = 0;
                break;
            }
        }
    }
    uint32_t fvar_offset = 0;
    for (uint32_t i = 0; i < 800; i++) {
        if (lin_pivot[i] == 0) {
//            printf("%d\n", i);
            free_var_idx[fvar_offset] = i;
            reverse_free_var_idx[i] = fvar_offset;
            fvar_offset++;
        } else {
            reverse_free_var_idx[i] = 0xFFFFFFFF;
        }
    }
    printf("dependent variable number: %d\n", fvar_offset);
}

void appendReduce() {
    uint32_t eq_idx, var_idx, i, j;

    for (eq_idx = 0; eq_idx < AEQ_NUM; eq_idx++) {
        memset(reduced_append_system[eq_idx], 0, AVAR_XVAR_NUM);
        uint32_t multix_1;
        uint32_t multix_2;
        for (multix_1 = 0; multix_1 < 800; multix_1++) {
            bool is_multix1_dep = (bool) (reverse_free_var_idx[multix_1] == 0xFFFFFFFF);

            // modify quadratic part
            for (multix_2 = multix_1; multix_2 < 800; multix_2++) {
                var_idx = deg2midx2(multix_1, multix_2);
                if (append_system_origin[eq_idx][var_idx]) {
                    bool is_multix2_dep = (bool) (reverse_free_var_idx[multix_2] == 0xFFFFFFFF);

                    if (is_multix1_dep && is_multix2_dep) {
                        for (i = 0; i < 800; i++) {
                            for (j = 0; j < 800; j++) {
                                if (lin_dep[multix_1][i] && lin_dep[multix_2][j]) {
                                    reduced_append_system[eq_idx][deg2midx2(j, i)] ^= 1;
                                }
                            }
                        }

                        if (lin_dep[multix_1][800]) {
                            for (i = 0; i < 800; i++) {
                                if (lin_dep[multix_2][i])
                                    reduced_append_system[eq_idx][deg2midx1(800, i)] ^= 1;
                            }
                        }

                        if (lin_dep[multix_2][800]) {
                            for (i = 0; i < 800; i++) {
                                if (lin_dep[multix_1][i])
                                    reduced_append_system[eq_idx][deg2midx1(800, i)] ^= 1;
                            }
                        }

                        if (lin_dep[multix_1][800] && lin_dep[multix_2][800])
                            reduced_append_system[eq_idx][AVAR_XVAR_NUM - 1] ^= 1;

                    } else if (is_multix1_dep) {
                        for (i = 0; i < 800; i++) {
                            if (lin_dep[multix_1][i]) {
                                reduced_append_system[eq_idx][deg2midx2(multix_2, i)] ^= 1;
                            }
                        }

                        if (lin_dep[multix_1][800])
                            reduced_append_system[eq_idx][deg2midx1(800, multix_2)] ^= 1;

                    } else if (is_multix2_dep) {
                        for (i = 0; i < 800; i++) {
                            if (lin_dep[multix_2][i]) {
                                reduced_append_system[eq_idx][deg2midx2(multix_1, i)] ^= 1;
                            }
                        }

                        if (lin_dep[multix_2][800])
                            reduced_append_system[eq_idx][deg2midx1(800, multix_1)] ^= 1;
                    } else {
                        reduced_append_system[eq_idx][deg2midx2(multix_1, multix_2)] ^= 1;
                    }
                }
            }
            if (append_system_origin[eq_idx][deg2midx1(800, multix_1)]) {
                if (is_multix1_dep) {
                    for (i = 0; i < 801; i++) {
                        if (lin_dep[multix_1][i])
                            reduced_append_system[eq_idx][deg2midx1(800, i)] ^= 1;
                    }
                } else {
                    reduced_append_system[eq_idx][deg2midx1(800, multix_1)] ^= 1;
                }
            }
        }
        reduced_append_system[eq_idx][AVAR_XVAR_NUM - 1] ^= append_system_origin[eq_idx][AVAR_XVAR_NUM - 1];
        for (i = 0; i < 800; i++)
            for (j = i; j < 800; j++)
                if (reduced_append_system[eq_idx][deg2midx2(i, j)])
                    printf("%d %d\n", i, j);
        exit(0);
    }
}

int main() {
    char ch;
    FILE *fmq = fopen("../data/append_analysis.dat", "r");
    if (fmq == NULL) {
        printf("File is not available\n");
    } else {
        for (uint32_t i = 0; i < AEQ_NUM; i++) {
            for (uint32_t j = 0; j < AVAR_XVAR_NUM; j++) {
                ch = fgetc(fmq);
                append_system_origin[i][j] = (uint8_t) (ch - '0');
            }
            fgetc(fmq);
        }
        fclose(fmq);
    }

    FILE *flin = fopen("../data/lin_analysis.dat", "r");
    if (flin == NULL) {
        printf("File is not available\n");
    } else {
        for (uint32_t i = 0; i < LINEQ_NUM; i++) {
            for (uint32_t j = 0; j < 801; j++) {
                ch = fgetc(flin);
                lin_system[i][j] = (uint8_t) (ch - '0');
            }
            fgetc(flin);
        }
        fclose(flin);
    }

    FILE *filin = fopen("../data/variate_analysis.dat", "r");
    if (filin == NULL) {
        printf("File is not available\n");
    } else {
        for (uint32_t i = 0; i < ILINEQ_NUM; i++) {
            for (uint32_t j = 0; j < 801; j++) {
                ch = fgetc(filin);
                it_lin_system[i][j] = (uint8_t) (ch - '0');
            }
            fgetc(filin);
        }
        fclose(filin);
    }

    linearReduce(0xFFFF);
//    for (uint32_t i = 0; i < 800; i++) {
//        if (reverse_free_var_idx[i] == 0xFFFFFFFF) {
//            printf("%d:", i);
//            for (uint32_t j = 0; j < 801; j++) {
//                if(lin_dep[i][j])
//                    printf(" %d ", j);
//            }
//            printf("\n");
//        }
//    }
    appendReduce();
    return 0;
}