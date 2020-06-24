//
// Created by vitowu on 2/19/20.
//

#include "status.h"

#include <unistd.h>

#define RC1 0x80000000
#define RC2 0x41010000
#define RC3 0x51010000
#define RC4 0x00010001
#define ROR32I(x, a) (((x)>>(a))|((x)<<(32-(a))))

void
reduce_lsub(bool *sys, uint64_t row_num, uint64_t col_num, uint64_t lcol_num) {
    assert(lcol_num <= col_num);
    const uint64_t bound = (row_num < lcol_num) ? row_num : lcol_num;
    uint64_t i, j, k;
    bool tmp_row[col_num];
    for (i = 0; i < bound; ++i) {
        // find the i-th pivot
        for (j = i; j < row_num; ++j) {
            if (true == sys[j * col_num + i]) { // sys[j][i]
                break;
            }
        }
        if (row_num == j) { // singular
            continue;
        }

        // swap i-th and j-th rows
        memcpy(tmp_row, sys + i * col_num, col_num);
        memcpy(sys + i * col_num, sys + j * col_num, col_num);
        memcpy(sys + j * col_num, tmp_row, col_num);

        // for all the rows below pivot
        for (j = i + 1; j < row_num; ++j) {
            if (true == sys[j * col_num + i]) { // sys[j][i]
                // subtract i-th row from the row
                for (k = 0; k < col_num; ++k) {
                    sys[j * col_num + k] ^= sys[i * col_num + k];
                }
            }
        }
    }
}

void dump_variate_constraints(std::vector<Keccak::Status> &constraints) {
    std::vector<uint32_t> status_idx = {5, 6, 7, 8, 9, 15, 16, 17, 18, 19};
    for (uint32_t i = 0; i < 10; i++) {
        for (uint32_t j = 0; j < 32; j++)
            if (constraints[i].var32[j].poly_ord != 0) {
//                std::cout << status_idx[i] << " " << j << " ";
                constraints[i].display_var_linpart(j);
            }
    }
}

int main() {
    std::vector<Keccak::Status> constrainedStatus;

//    std::cout << "Initialize" << std::endl;
    Keccak::Status A_1_0[5][5];
    Keccak::init_status(A_1_0);  // assign variable index

    Keccak::Status A_1_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_1_1[0][j] = A_1_0[0][j] ^ A_1_0[4][0] ^ A_1_0[4][1] ^ A_1_0[4][2] ^ A_1_0[4][3] ^ A_1_0[4][4]
                      ^ ROR32((A_1_0[1][0] ^ A_1_0[1][1] ^ A_1_0[1][2] ^ A_1_0[1][3] ^ A_1_0[1][4]), 1);

        A_1_1[1][j] = A_1_0[1][j] ^ A_1_0[0][0] ^ A_1_0[0][1] ^ A_1_0[0][2] ^ A_1_0[0][3] ^ A_1_0[0][4]
                      ^ ROR32((A_1_0[2][0] ^ A_1_0[2][1] ^ A_1_0[2][2] ^ A_1_0[2][3] ^ A_1_0[2][4]), 1);

        A_1_1[2][j] = A_1_0[2][j] ^ A_1_0[1][0] ^ A_1_0[1][1] ^ A_1_0[1][2] ^ A_1_0[1][3] ^ A_1_0[1][4]
                      ^ ROR32((A_1_0[3][0] ^ A_1_0[3][1] ^ A_1_0[3][2] ^ A_1_0[3][3] ^ A_1_0[3][4]), 1);

        A_1_1[3][j] = A_1_0[3][j] ^ A_1_0[2][0] ^ A_1_0[2][1] ^ A_1_0[2][2] ^ A_1_0[2][3] ^ A_1_0[2][4]
                      ^ ROR32((A_1_0[4][0] ^ A_1_0[4][1] ^ A_1_0[4][2] ^ A_1_0[4][3] ^ A_1_0[4][4]), 1);

        A_1_1[4][j] = A_1_0[4][j] ^ A_1_0[3][0] ^ A_1_0[3][1] ^ A_1_0[3][2] ^ A_1_0[3][3] ^ A_1_0[3][4]
                      ^ ROR32((A_1_0[0][0] ^ A_1_0[0][1] ^ A_1_0[0][2] ^ A_1_0[0][3] ^ A_1_0[0][4]), 1);
    }

    Keccak::Status A_1_2[5][5];
    A_1_2[0][0] = A_1_1[0][0];
    A_1_2[0][1] = ROR32(A_1_1[3][0], 28);
    A_1_2[0][2] = ROR32(A_1_1[1][0], 1);
    A_1_2[0][3] = ROR32(A_1_1[4][0], 27);
    A_1_2[0][4] = ROR32(A_1_1[2][0], 30);
    A_1_2[1][0] = ROR32(A_1_1[1][1], 12);
    A_1_2[1][1] = ROR32(A_1_1[4][1], 20);
    A_1_2[1][2] = ROR32(A_1_1[2][1], 6);
    A_1_2[1][3] = ROR32(A_1_1[0][1], 4);
    A_1_2[1][4] = ROR32(A_1_1[3][1], 23);
    A_1_2[2][0] = ROR32(A_1_1[2][2], 11);
    A_1_2[2][1] = ROR32(A_1_1[0][2], 3);
    A_1_2[2][2] = ROR32(A_1_1[3][2], 25);
    A_1_2[2][3] = ROR32(A_1_1[1][2], 10);
    A_1_2[2][4] = ROR32(A_1_1[4][2], 7);
    A_1_2[3][0] = ROR32(A_1_1[3][3], 21);
    A_1_2[3][1] = ROR32(A_1_1[1][3], 13);
    A_1_2[3][2] = ROR32(A_1_1[4][3], 8);
    A_1_2[3][3] = ROR32(A_1_1[2][3], 15);
    A_1_2[3][4] = ROR32(A_1_1[0][3], 9);
    A_1_2[4][0] = ROR32(A_1_1[4][4], 14);
    A_1_2[4][1] = ROR32(A_1_1[2][4], 29);
    A_1_2[4][2] = ROR32(A_1_1[0][4], 18);
    A_1_2[4][3] = ROR32(A_1_1[3][4], 24);
    A_1_2[4][4] = ROR32(A_1_1[1][4], 2);

    Keccak::Status A_2_0[5][5];
    A_2_0[0][0] = A_1_2[0][0] ^ RC1;
    Keccak::status_set_val(A_2_0[1][0], 0xFFFFFFFF);
    A_2_0[2][0] = A_1_2[0][0] ^ A_1_2[2][0];
    Keccak::status_set_val(A_2_0[3][0], 0);
    Keccak::status_set_val(A_2_0[4][0], 0xFFFFFFFF);
    A_2_0[0][1] = A_1_2[0][1];
    Keccak::status_set_val(A_2_0[1][1], 0xFFFFFFFF);
    A_2_0[2][1] = A_1_2[0][1] ^ A_1_2[2][1];
    Keccak::status_set_val(A_2_0[3][1], 0);
    Keccak::status_set_val(A_2_0[4][1], 0xFFFFFFFF);
    A_2_0[0][2] = A_1_2[0][2];
    Keccak::status_set_val(A_2_0[1][2], 0xFFFFFFFF);
    A_2_0[2][2] = A_1_2[0][2] ^ A_1_2[2][2];
    Keccak::status_set_val(A_2_0[3][2], 0);
    Keccak::status_set_val(A_2_0[4][2], 0xFFFFFFFF);
    A_2_0[0][3] = A_1_2[0][3];
    Keccak::status_set_val(A_2_0[1][3], 0xFFFFFFFF);
    A_2_0[2][3] = A_1_2[0][3] ^ A_1_2[2][3];
    Keccak::status_set_val(A_2_0[3][3], 0);
    Keccak::status_set_val(A_2_0[4][3], 0xFFFFFFFF);
    A_2_0[0][4] = A_1_2[0][4];
    Keccak::status_set_val(A_2_0[1][4], 0xFFFFFFFF);
    A_2_0[2][4] = A_1_2[0][4] ^ A_1_2[2][4];
    Keccak::status_set_val(A_2_0[3][4], 0);
    Keccak::status_set_val(A_2_0[4][4], 0xFFFFFFFF);

    // set 2 constraints
//    const uint32_t Alpha = 0x44e72;
//    const uint32_t Beta = 0xBAA20F;
    const uint32_t Alpha = 0x80000000;
    const uint32_t Beta = 0x00000000;

    Keccak::Status A_2_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_2_1[0][j] = A_2_0[0][j];
        Keccak::status_set_val(A_2_1[1][j], 0xFFFFFFFF ^ Alpha ^ ROR32I(Beta, 1));
        A_2_1[2][j] = A_2_0[2][j] ^ 0xFFFFFFFF;
        Keccak::status_set_val(A_2_1[3][j], Beta ^ 0xFFFFFFFF);
        Keccak::status_set_val(A_2_1[4][j], 0xFFFFFFFF ^ ROR32I(Alpha, 1));
    }

    Keccak::Status A_2_2[5][5];
    A_2_2[0][0] = A_2_1[0][0];
    A_2_2[0][1] = ROR32(A_2_1[3][0], 28);
    A_2_2[0][2] = ROR32(A_2_1[1][0], 1);
    A_2_2[0][3] = ROR32(A_2_1[4][0], 27);
    A_2_2[0][4] = ROR32(A_2_1[2][0], 30);
    A_2_2[1][0] = ROR32(A_2_1[1][1], 12);
    A_2_2[1][1] = ROR32(A_2_1[4][1], 20);
    A_2_2[1][2] = ROR32(A_2_1[2][1], 6);
    A_2_2[1][3] = ROR32(A_2_1[0][1], 4);
    A_2_2[1][4] = ROR32(A_2_1[3][1], 23);
    A_2_2[2][0] = ROR32(A_2_1[2][2], 11);
    A_2_2[2][1] = ROR32(A_2_1[0][2], 3);
    A_2_2[2][2] = ROR32(A_2_1[3][2], 25);
    A_2_2[2][3] = ROR32(A_2_1[1][2], 10);
    A_2_2[2][4] = ROR32(A_2_1[4][2], 7);
    A_2_2[3][0] = ROR32(A_2_1[3][3], 21);
    A_2_2[3][1] = ROR32(A_2_1[1][3], 13);
    A_2_2[3][2] = ROR32(A_2_1[4][3], 8);
    A_2_2[3][3] = ROR32(A_2_1[2][3], 15);
    A_2_2[3][4] = ROR32(A_2_1[0][3], 9);
    A_2_2[4][0] = ROR32(A_2_1[4][4], 14);
    A_2_2[4][1] = ROR32(A_2_1[2][4], 29);
    A_2_2[4][2] = ROR32(A_2_1[0][4], 18);
    A_2_2[4][3] = ROR32(A_2_1[3][4], 24);
    A_2_2[4][4] = ROR32(A_2_1[1][4], 2);

    Keccak::Status A_3_0[5][5];
    for (int j = 0; j < 5; j++) {
        A_3_0[0][j] = A_2_2[0][j] ^ ((A_2_2[1][j] ^ 0xFFFFFFFF) & A_2_2[2][j]);
        A_3_0[1][j] = A_2_2[1][j] ^ ((A_2_2[2][j] ^ 0xFFFFFFFF) & A_2_2[3][j]);
        A_3_0[2][j] = A_2_2[2][j] ^ ((A_2_2[3][j] ^ 0xFFFFFFFF) & A_2_2[4][j]);
        A_3_0[3][j] = A_2_2[3][j] ^ ((A_2_2[4][j] ^ 0xFFFFFFFF) & A_2_2[0][j]);
        A_3_0[4][j] = A_2_2[4][j] ^ ((A_2_2[0][j] ^ 0xFFFFFFFF) & A_2_2[1][j]);
    }
    A_3_0[0][0] = A_3_0[0][0] ^ RC2;

    Keccak::Status A_3_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_3_1[0][j] = A_3_0[0][j] ^ A_3_0[4][0] ^ A_3_0[4][1] ^ A_3_0[4][2] ^ A_3_0[4][3] ^ A_3_0[4][4]
                      ^ ROR32((A_3_0[1][0] ^ A_3_0[1][1] ^ A_3_0[1][2] ^ A_3_0[1][3] ^ A_3_0[1][4]), 1);
        A_3_1[1][j] = A_3_0[1][j] ^ A_3_0[0][0] ^ A_3_0[0][1] ^ A_3_0[0][2] ^ A_3_0[0][3] ^ A_3_0[0][4]
                      ^ ROR32((A_3_0[2][0] ^ A_3_0[2][1] ^ A_3_0[2][2] ^ A_3_0[2][3] ^ A_3_0[2][4]), 1);
        A_3_1[2][j] = A_3_0[2][j] ^ A_3_0[1][0] ^ A_3_0[1][1] ^ A_3_0[1][2] ^ A_3_0[1][3] ^ A_3_0[1][4]
                      ^ ROR32((A_3_0[3][0] ^ A_3_0[3][1] ^ A_3_0[3][2] ^ A_3_0[3][3] ^ A_3_0[3][4]), 1);
        A_3_1[3][j] = A_3_0[3][j] ^ A_3_0[2][0] ^ A_3_0[2][1] ^ A_3_0[2][2] ^ A_3_0[2][3] ^ A_3_0[2][4]
                      ^ ROR32((A_3_0[4][0] ^ A_3_0[4][1] ^ A_3_0[4][2] ^ A_3_0[4][3] ^ A_3_0[4][4]), 1);
        A_3_1[4][j] = A_3_0[4][j] ^ A_3_0[3][0] ^ A_3_0[3][1] ^ A_3_0[3][2] ^ A_3_0[3][3] ^ A_3_0[3][4]
                      ^ ROR32((A_3_0[0][0] ^ A_3_0[0][1] ^ A_3_0[0][2] ^ A_3_0[0][3] ^ A_3_0[0][4]), 1);
    }

    Keccak::Status A_3_2[5][5];
    A_3_2[0][0] = A_3_1[0][0];
    A_3_2[0][1] = ROR32(A_3_1[3][0], 28);
    A_3_2[0][2] = ROR32(A_3_1[1][0], 1);
    A_3_2[0][3] = ROR32(A_3_1[4][0], 27);
    A_3_2[0][4] = ROR32(A_3_1[2][0], 30);

    A_3_2[1][0] = ROR32(A_3_1[1][1], 12);
    A_3_2[1][1] = ROR32(A_3_1[4][1], 20);
    A_3_2[1][2] = ROR32(A_3_1[2][1], 6);
    A_3_2[1][3] = ROR32(A_3_1[0][1], 4);
    A_3_2[1][4] = ROR32(A_3_1[3][1], 23);

    A_3_2[2][0] = ROR32(A_3_1[2][2], 11);
    A_3_2[2][1] = ROR32(A_3_1[0][2], 3);
    A_3_2[2][2] = ROR32(A_3_1[3][2], 25);
    A_3_2[2][3] = ROR32(A_3_1[1][2], 10);
    A_3_2[2][4] = ROR32(A_3_1[4][2], 7);

    A_3_2[3][0] = ROR32(A_3_1[3][3], 21);
    A_3_2[3][1] = ROR32(A_3_1[1][3], 13);
    A_3_2[3][2] = ROR32(A_3_1[4][3], 8);
    A_3_2[3][3] = ROR32(A_3_1[2][3], 15);
    A_3_2[3][4] = ROR32(A_3_1[0][3], 9);

    A_3_2[4][0] = ROR32(A_3_1[4][4], 14);
    A_3_2[4][1] = ROR32(A_3_1[2][4], 29);
    A_3_2[4][2] = ROR32(A_3_1[0][4], 18);
    A_3_2[4][3] = ROR32(A_3_1[3][4], 24);
    A_3_2[4][4] = ROR32(A_3_1[1][4], 2);

//    std::vector<Keccak::Status> varc;
//    varc.push_back(A_3_2[1][0] & 0x00000746);
//    varc.push_back(A_3_2[1][1] & 0x00000746);
//    varc.push_back(A_3_2[1][2] & 0x00000746);
//    varc.push_back(A_3_2[1][3] & 0x00000746);
//    varc.push_back(A_3_2[1][4] & 0x00000746);
//    varc.push_back(A_3_2[3][0] & 0x00000E8C);
//    varc.push_back(A_3_2[3][1] & 0x00000FCC);
//    varc.push_back(A_3_2[3][2] & 0x00000E8C);
//    varc.push_back(A_3_2[3][3] & 0x00000E8C);
//    varc.push_back(A_3_2[3][4] & 0x00000E8C);
//    dump_variate_constraints(varc);

// IMPORTANT:
// --------- for test use round 3 constraints ----------
    std::vector<Keccak::Status> varc;
    varc.push_back(A_3_2[1][0] & 0xA0000026);
    varc.push_back(A_3_2[1][1] & 0xA0000026);
    varc.push_back(A_3_2[1][2] & 0xA0000026);
    varc.push_back(A_3_2[1][3] & 0xA0000026);
    varc.push_back(A_3_2[1][4] & 0xA0000026);
    varc.push_back(A_3_2[3][0] & 0x4000004D);
    varc.push_back(A_3_2[3][1] & 0xE000006F);
    varc.push_back(A_3_2[3][2] & 0x4000004D);
    varc.push_back(A_3_2[3][3] & 0x4000004D);
    varc.push_back(A_3_2[3][4] & 0x4000004D);
    dump_variate_constraints(varc);
// -----------------------------------------------------
}