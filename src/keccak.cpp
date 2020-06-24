//
// Created by vitowu on 3/18/20.
//

#include "keccak.h"

#define ALPHA 0x44E72
#define BETA 0xBAA20F
#define DIFF_TOLERANCE 4
#define RC1 0x80000000
#define RC2 0x41010000
#define RC3 0x51010000
#define RC4 0x00010001

bool cpu_VerifyKeccakResult(const uint32_t A[5][5], uint32_t *minDiff) {
#define _ROR32(x, a) ((((x)>>(a))|((x)<<(32-(a)))))
#define u32 uint32_t
    u32 A_1_0[5][5];
    for (uint32_t i = 0; i < 5; i++)
        memcpy(A_1_0[i], A[i], sizeof(uint32_t) * 5);

    u32 A_1_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_1_1[0][j] = A_1_0[0][j] ^ A_1_0[4][0] ^ A_1_0[4][1] ^ A_1_0[4][2] ^ A_1_0[4][3] ^ A_1_0[4][4]
                      ^ _ROR32((A_1_0[1][0] ^ A_1_0[1][1] ^ A_1_0[1][2] ^ A_1_0[1][3] ^ A_1_0[1][4]), 1);
        A_1_1[1][j] = A_1_0[1][j] ^ A_1_0[0][0] ^ A_1_0[0][1] ^ A_1_0[0][2] ^ A_1_0[0][3] ^ A_1_0[0][4]
                      ^ _ROR32((A_1_0[2][0] ^ A_1_0[2][1] ^ A_1_0[2][2] ^ A_1_0[2][3] ^ A_1_0[2][4]), 1);
        A_1_1[2][j] = A_1_0[2][j] ^ A_1_0[1][0] ^ A_1_0[1][1] ^ A_1_0[1][2] ^ A_1_0[1][3] ^ A_1_0[1][4]
                      ^ _ROR32((A_1_0[3][0] ^ A_1_0[3][1] ^ A_1_0[3][2] ^ A_1_0[3][3] ^ A_1_0[3][4]), 1);
        A_1_1[3][j] = A_1_0[3][j] ^ A_1_0[2][0] ^ A_1_0[2][1] ^ A_1_0[2][2] ^ A_1_0[2][3] ^ A_1_0[2][4]
                      ^ _ROR32((A_1_0[4][0] ^ A_1_0[4][1] ^ A_1_0[4][2] ^ A_1_0[4][3] ^ A_1_0[4][4]), 1);
        A_1_1[4][j] = A_1_0[4][j] ^ A_1_0[3][0] ^ A_1_0[3][1] ^ A_1_0[3][2] ^ A_1_0[3][3] ^ A_1_0[3][4]
                      ^ _ROR32((A_1_0[0][0] ^ A_1_0[0][1] ^ A_1_0[0][2] ^ A_1_0[0][3] ^ A_1_0[0][4]), 1);
    }

    u32 A_1_2[5][5];
    A_1_2[0][0] = A_1_1[0][0];
    A_1_2[0][1] = _ROR32(A_1_1[3][0], 28);
    A_1_2[0][2] = _ROR32(A_1_1[1][0], 1);
    A_1_2[0][3] = _ROR32(A_1_1[4][0], 27);
    A_1_2[0][4] = _ROR32(A_1_1[2][0], 30);
    A_1_2[1][0] = _ROR32(A_1_1[1][1], 12);
    A_1_2[1][1] = _ROR32(A_1_1[4][1], 20);
    A_1_2[1][2] = _ROR32(A_1_1[2][1], 6);
    A_1_2[1][3] = _ROR32(A_1_1[0][1], 4);
    A_1_2[1][4] = _ROR32(A_1_1[3][1], 23);
    A_1_2[2][0] = _ROR32(A_1_1[2][2], 11);
    A_1_2[2][1] = _ROR32(A_1_1[0][2], 3);
    A_1_2[2][2] = _ROR32(A_1_1[3][2], 25);
    A_1_2[2][3] = _ROR32(A_1_1[1][2], 10);
    A_1_2[2][4] = _ROR32(A_1_1[4][2], 7);
    A_1_2[3][0] = _ROR32(A_1_1[3][3], 21);
    A_1_2[3][1] = _ROR32(A_1_1[1][3], 13);
    A_1_2[3][2] = _ROR32(A_1_1[4][3], 8);
    A_1_2[3][3] = _ROR32(A_1_1[2][3], 15);
    A_1_2[3][4] = _ROR32(A_1_1[0][3], 9);
    A_1_2[4][0] = _ROR32(A_1_1[4][4], 14);
    A_1_2[4][1] = _ROR32(A_1_1[2][4], 29);
    A_1_2[4][2] = _ROR32(A_1_1[0][4], 18);
    A_1_2[4][3] = _ROR32(A_1_1[3][4], 24);
    A_1_2[4][4] = _ROR32(A_1_1[1][4], 2);

    u32 A_2_0[5][5];
    A_2_0[0][0] = A_1_2[0][0] ^ RC1;
    A_2_0[1][0] = 0xFFFFFFFF;
    A_2_0[2][0] = A_1_2[0][0] ^ A_1_2[2][0];
    A_2_0[3][0] = 0;
    A_2_0[4][0] = 0xFFFFFFFF;
    A_2_0[0][1] = A_1_2[0][1];
    A_2_0[1][1] = 0xFFFFFFFF;
    A_2_0[2][1] = A_1_2[0][1] ^ A_1_2[2][1];
    A_2_0[3][1] = 0;
    A_2_0[4][1] = 0xFFFFFFFF;
    A_2_0[0][2] = A_1_2[0][2];
    A_2_0[1][2] = 0xFFFFFFFF;
    A_2_0[2][2] = A_1_2[0][2] ^ A_1_2[2][2];
    A_2_0[3][2] = 0;
    A_2_0[4][2] = 0xFFFFFFFF;
    A_2_0[0][3] = A_1_2[0][3];
    A_2_0[1][3] = 0xFFFFFFFF;
    A_2_0[2][3] = A_1_2[0][3] ^ A_1_2[2][3];
    A_2_0[3][3] = 0;
    A_2_0[4][3] = 0xFFFFFFFF;
    A_2_0[0][4] = A_1_2[0][4];
    A_2_0[1][4] = 0xFFFFFFFF;
    A_2_0[2][4] = A_1_2[0][4] ^ A_1_2[2][4];
    A_2_0[3][4] = 0;
    A_2_0[4][4] = 0xFFFFFFFF;

#ifndef TEST_PRE
    if ((A_2_0[0][0] ^ A_2_0[0][1] ^ A_2_0[0][2] ^ A_2_0[0][3] ^ A_2_0[0][4]) != ALPHA)
        EXIT_WITH_MSG("round 2 constraint error\n");
    if ((A_2_0[2][0] ^ A_2_0[2][1] ^ A_2_0[2][2] ^ A_2_0[2][3] ^ A_2_0[2][4]) != BETA)
        EXIT_WITH_MSG("round 2 constraint error\n");
#endif

    u32 A_2_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_2_1[0][j] = A_2_0[0][j] ^ A_2_0[4][0] ^ A_2_0[4][1] ^ A_2_0[4][2] ^ A_2_0[4][3] ^ A_2_0[4][4]
                      ^ ROR32((A_2_0[1][0] ^ A_2_0[1][1] ^ A_2_0[1][2] ^ A_2_0[1][3] ^ A_2_0[1][4]), 1);
        A_2_1[1][j] = A_2_0[1][j] ^ A_2_0[0][0] ^ A_2_0[0][1] ^ A_2_0[0][2] ^ A_2_0[0][3] ^ A_2_0[0][4]
                      ^ ROR32((A_2_0[2][0] ^ A_2_0[2][1] ^ A_2_0[2][2] ^ A_2_0[2][3] ^ A_2_0[2][4]), 1);
        A_2_1[2][j] = A_2_0[2][j] ^ A_2_0[1][0] ^ A_2_0[1][1] ^ A_2_0[1][2] ^ A_2_0[1][3] ^ A_2_0[1][4]
                      ^ ROR32((A_2_0[3][0] ^ A_2_0[3][1] ^ A_2_0[3][2] ^ A_2_0[3][3] ^ A_2_0[3][4]), 1);
        A_2_1[3][j] = A_2_0[3][j] ^ A_2_0[2][0] ^ A_2_0[2][1] ^ A_2_0[2][2] ^ A_2_0[2][3] ^ A_2_0[2][4]
                      ^ ROR32((A_2_0[4][0] ^ A_2_0[4][1] ^ A_2_0[4][2] ^ A_2_0[4][3] ^ A_2_0[4][4]), 1);
        A_2_1[4][j] = A_2_0[4][j] ^ A_2_0[3][0] ^ A_2_0[3][1] ^ A_2_0[3][2] ^ A_2_0[3][3] ^ A_2_0[3][4]
                      ^ ROR32((A_2_0[0][0] ^ A_2_0[0][1] ^ A_2_0[0][2] ^ A_2_0[0][3] ^ A_2_0[0][4]), 1);
    }

    u32 A_2_2[5][5];
    A_2_2[0][0] = A_2_1[0][0];
    A_2_2[0][1] = _ROR32(A_2_1[3][0], 28);
    A_2_2[0][2] = _ROR32(A_2_1[1][0], 1);
    A_2_2[0][3] = _ROR32(A_2_1[4][0], 27);
    A_2_2[0][4] = _ROR32(A_2_1[2][0], 30);
    A_2_2[1][0] = _ROR32(A_2_1[1][1], 12);
    A_2_2[1][1] = _ROR32(A_2_1[4][1], 20);
    A_2_2[1][2] = _ROR32(A_2_1[2][1], 6);
    A_2_2[1][3] = _ROR32(A_2_1[0][1], 4);
    A_2_2[1][4] = _ROR32(A_2_1[3][1], 23);
    A_2_2[2][0] = _ROR32(A_2_1[2][2], 11);
    A_2_2[2][1] = _ROR32(A_2_1[0][2], 3);
    A_2_2[2][2] = _ROR32(A_2_1[3][2], 25);
    A_2_2[2][3] = _ROR32(A_2_1[1][2], 10);
    A_2_2[2][4] = _ROR32(A_2_1[4][2], 7);
    A_2_2[3][0] = _ROR32(A_2_1[3][3], 21);
    A_2_2[3][1] = _ROR32(A_2_1[1][3], 13);
    A_2_2[3][2] = _ROR32(A_2_1[4][3], 8);
    A_2_2[3][3] = _ROR32(A_2_1[2][3], 15);
    A_2_2[3][4] = _ROR32(A_2_1[0][3], 9);
    A_2_2[4][0] = _ROR32(A_2_1[4][4], 14);
    A_2_2[4][1] = _ROR32(A_2_1[2][4], 29);
    A_2_2[4][2] = _ROR32(A_2_1[0][4], 18);
    A_2_2[4][3] = _ROR32(A_2_1[3][4], 24);
    A_2_2[4][4] = _ROR32(A_2_1[1][4], 2);

    //3rd round
    u32 A_3_0[5][5];
    for (int j = 0; j < 5; j++) {
        A_3_0[0][j] = A_2_2[0][j] ^ ((A_2_2[1][j] ^ 0xFFFFFFFF) & A_2_2[2][j]);
        A_3_0[1][j] = A_2_2[1][j] ^ ((A_2_2[2][j] ^ 0xFFFFFFFF) & A_2_2[3][j]);
        A_3_0[2][j] = A_2_2[2][j] ^ ((A_2_2[3][j] ^ 0xFFFFFFFF) & A_2_2[4][j]);
        A_3_0[3][j] = A_2_2[3][j] ^ ((A_2_2[4][j] ^ 0xFFFFFFFF) & A_2_2[0][j]);
        A_3_0[4][j] = A_2_2[4][j] ^ ((A_2_2[0][j] ^ 0xFFFFFFFF) & A_2_2[1][j]);
    }
    A_3_0[0][0] ^= RC2;

    u32 A_3_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_3_1[0][j] = A_3_0[0][j] ^ A_3_0[4][0] ^ A_3_0[4][1] ^ A_3_0[4][2] ^ A_3_0[4][3] ^ A_3_0[4][4]
                      ^ _ROR32((A_3_0[1][0] ^ A_3_0[1][1] ^ A_3_0[1][2] ^ A_3_0[1][3] ^ A_3_0[1][4]), 1);
        A_3_1[1][j] = A_3_0[1][j] ^ A_3_0[0][0] ^ A_3_0[0][1] ^ A_3_0[0][2] ^ A_3_0[0][3] ^ A_3_0[0][4]
                      ^ _ROR32((A_3_0[2][0] ^ A_3_0[2][1] ^ A_3_0[2][2] ^ A_3_0[2][3] ^ A_3_0[2][4]), 1);
        A_3_1[2][j] = A_3_0[2][j] ^ A_3_0[1][0] ^ A_3_0[1][1] ^ A_3_0[1][2] ^ A_3_0[1][3] ^ A_3_0[1][4]
                      ^ _ROR32((A_3_0[3][0] ^ A_3_0[3][1] ^ A_3_0[3][2] ^ A_3_0[3][3] ^ A_3_0[3][4]), 1);
        A_3_1[3][j] = A_3_0[3][j] ^ A_3_0[2][0] ^ A_3_0[2][1] ^ A_3_0[2][2] ^ A_3_0[2][3] ^ A_3_0[2][4]
                      ^ _ROR32((A_3_0[4][0] ^ A_3_0[4][1] ^ A_3_0[4][2] ^ A_3_0[4][3] ^ A_3_0[4][4]), 1);
        A_3_1[4][j] = A_3_0[4][j] ^ A_3_0[3][0] ^ A_3_0[3][1] ^ A_3_0[3][2] ^ A_3_0[3][3] ^ A_3_0[3][4]
                      ^ _ROR32((A_3_0[0][0] ^ A_3_0[0][1] ^ A_3_0[0][2] ^ A_3_0[0][3] ^ A_3_0[0][4]), 1);
    }

    u32 A_3_2[5][5];
    A_3_2[0][0] = A_3_1[0][0];
    A_3_2[0][1] = _ROR32(A_3_1[3][0], 28);
    A_3_2[0][2] = _ROR32(A_3_1[1][0], 1);
    A_3_2[0][3] = _ROR32(A_3_1[4][0], 27);
    A_3_2[0][4] = _ROR32(A_3_1[2][0], 30);
    A_3_2[1][0] = _ROR32(A_3_1[1][1], 12);
    A_3_2[1][1] = _ROR32(A_3_1[4][1], 20);
    A_3_2[1][2] = _ROR32(A_3_1[2][1], 6);
    A_3_2[1][3] = _ROR32(A_3_1[0][1], 4);
    A_3_2[1][4] = _ROR32(A_3_1[3][1], 23);
    A_3_2[2][0] = _ROR32(A_3_1[2][2], 11);
    A_3_2[2][1] = _ROR32(A_3_1[0][2], 3);
    A_3_2[2][2] = _ROR32(A_3_1[3][2], 25);
    A_3_2[2][3] = _ROR32(A_3_1[1][2], 10);
    A_3_2[2][4] = _ROR32(A_3_1[4][2], 7);
    A_3_2[3][0] = _ROR32(A_3_1[3][3], 21);
    A_3_2[3][1] = _ROR32(A_3_1[1][3], 13);
    A_3_2[3][2] = _ROR32(A_3_1[4][3], 8);
    A_3_2[3][3] = _ROR32(A_3_1[2][3], 15);
    A_3_2[3][4] = _ROR32(A_3_1[0][3], 9);
    A_3_2[4][0] = _ROR32(A_3_1[4][4], 14);
    A_3_2[4][1] = _ROR32(A_3_1[2][4], 29);
    A_3_2[4][2] = _ROR32(A_3_1[0][4], 18);
    A_3_2[4][3] = _ROR32(A_3_1[3][4], 24);
    A_3_2[4][4] = _ROR32(A_3_1[1][4], 2);

    //4th round
    u32 A_4_0[5][5];
    for (int j = 0; j < 5; j++) {
        A_4_0[0][j] = A_3_2[0][j] ^ ((A_3_2[1][j] ^ 0xFFFFFFFF) & A_3_2[2][j]);
        A_4_0[1][j] = A_3_2[1][j] ^ ((A_3_2[2][j] ^ 0xFFFFFFFF) & A_3_2[3][j]);
        A_4_0[2][j] = A_3_2[2][j] ^ ((A_3_2[3][j] ^ 0xFFFFFFFF) & A_3_2[4][j]);
        A_4_0[3][j] = A_3_2[3][j] ^ ((A_3_2[4][j] ^ 0xFFFFFFFF) & A_3_2[0][j]);
        A_4_0[4][j] = A_3_2[4][j] ^ ((A_3_2[0][j] ^ 0xFFFFFFFF) & A_3_2[1][j]);
    }
    A_4_0[0][0] ^= RC3;

    u32 A_4_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_4_1[0][j] = A_4_0[0][j] ^ A_4_0[4][0] ^ A_4_0[4][1] ^ A_4_0[4][2] ^ A_4_0[4][3] ^ A_4_0[4][4]
                      ^ _ROR32((A_4_0[1][0] ^ A_4_0[1][1] ^ A_4_0[1][2] ^ A_4_0[1][3] ^ A_4_0[1][4]), 1);
        A_4_1[1][j] = A_4_0[1][j] ^ A_4_0[0][0] ^ A_4_0[0][1] ^ A_4_0[0][2] ^ A_4_0[0][3] ^ A_4_0[0][4]
                      ^ _ROR32((A_4_0[2][0] ^ A_4_0[2][1] ^ A_4_0[2][2] ^ A_4_0[2][3] ^ A_4_0[2][4]), 1);
        A_4_1[2][j] = A_4_0[2][j] ^ A_4_0[1][0] ^ A_4_0[1][1] ^ A_4_0[1][2] ^ A_4_0[1][3] ^ A_4_0[1][4]
                      ^ _ROR32((A_4_0[3][0] ^ A_4_0[3][1] ^ A_4_0[3][2] ^ A_4_0[3][3] ^ A_4_0[3][4]), 1);
        A_4_1[3][j] = A_4_0[3][j] ^ A_4_0[2][0] ^ A_4_0[2][1] ^ A_4_0[2][2] ^ A_4_0[2][3] ^ A_4_0[2][4]
                      ^ _ROR32((A_4_0[4][0] ^ A_4_0[4][1] ^ A_4_0[4][2] ^ A_4_0[4][3] ^ A_4_0[4][4]), 1);
        A_4_1[4][j] = A_4_0[4][j] ^ A_4_0[3][0] ^ A_4_0[3][1] ^ A_4_0[3][2] ^ A_4_0[3][3] ^ A_4_0[3][4]
                      ^ _ROR32((A_4_0[0][0] ^ A_4_0[0][1] ^ A_4_0[0][2] ^ A_4_0[0][3] ^ A_4_0[0][4]), 1);
    }

    u32 A_4_2[5][5];
    A_4_2[0][0] = A_4_1[0][0];
    A_4_2[0][1] = _ROR32(A_4_1[3][0], 28);
    A_4_2[0][2] = _ROR32(A_4_1[1][0], 1);
    A_4_2[0][3] = _ROR32(A_4_1[4][0], 27);
    A_4_2[0][4] = _ROR32(A_4_1[2][0], 30);
    A_4_2[1][0] = _ROR32(A_4_1[1][1], 12);
    A_4_2[1][1] = _ROR32(A_4_1[4][1], 20);
    A_4_2[1][2] = _ROR32(A_4_1[2][1], 6);
    A_4_2[1][3] = _ROR32(A_4_1[0][1], 4);
    A_4_2[1][4] = _ROR32(A_4_1[3][1], 23);
    A_4_2[2][0] = _ROR32(A_4_1[2][2], 11);
    A_4_2[2][1] = _ROR32(A_4_1[0][2], 3);
    A_4_2[2][2] = _ROR32(A_4_1[3][2], 25);
    A_4_2[2][3] = _ROR32(A_4_1[1][2], 10);
    A_4_2[2][4] = _ROR32(A_4_1[4][2], 7);
    A_4_2[3][0] = _ROR32(A_4_1[3][3], 21);
    A_4_2[3][1] = _ROR32(A_4_1[1][3], 13);
    A_4_2[3][2] = _ROR32(A_4_1[4][3], 8);
    A_4_2[3][3] = _ROR32(A_4_1[2][3], 15);
    A_4_2[3][4] = _ROR32(A_4_1[0][3], 9);
    A_4_2[4][0] = _ROR32(A_4_1[4][4], 14);
    A_4_2[4][1] = _ROR32(A_4_1[2][4], 29);
    A_4_2[4][2] = _ROR32(A_4_1[0][4], 18);
    A_4_2[4][3] = _ROR32(A_4_1[3][4], 24);
    A_4_2[4][4] = _ROR32(A_4_1[1][4], 2);

#ifndef TEST_PRE
    if ((A_4_2[0][0] & 0x00000744) != 0x00000004)
        EXIT_WITH_MSG("mq constraint error\n");
    if ((A_4_2[1][0] & 0x74400000) != 0x24000000)
        EXIT_WITH_MSG("mq constraint error\n");

    if ((A_4_2[0][0] & 0x27A98003) != 0x26090002)
        PRINTF_ERR("mq constraint error\n");
    if (((A_4_2[0][0] ^ A_4_2[2][0]) & 0xD85678B8) != 0x885068A0)
        PRINTF_ERR("mq constraint error\n");
    if ((A_4_2[1][0] & 0x8B040000) != 0x03000000)
        PRINTF_ERR("mq constraint error\n");
    if (((A_4_2[1][0] ^ A_4_2[3][0]) & 0x00BB0000) != 0x00A90000)
        PRINTF_ERR("mq constraint error\n");
#endif
    
    u32 A_5_0[5][5];
    for (int j = 0; j < 5; j++) {
        A_5_0[0][j] = A_4_2[0][j] ^ ((A_4_2[1][j] ^ 0xFFFFFFFF) & A_4_2[2][j]);
        A_5_0[1][j] = A_4_2[1][j] ^ ((A_4_2[2][j] ^ 0xFFFFFFFF) & A_4_2[3][j]);
        A_5_0[2][j] = A_4_2[2][j] ^ ((A_4_2[3][j] ^ 0xFFFFFFFF) & A_4_2[4][j]);
        A_5_0[3][j] = A_4_2[3][j] ^ ((A_4_2[4][j] ^ 0xFFFFFFFF) & A_4_2[0][j]);
        A_5_0[4][j] = A_4_2[4][j] ^ ((A_4_2[0][j] ^ 0xFFFFFFFF) & A_4_2[1][j]);
    }
    A_5_0[0][0] ^= RC4;

    u32 result[3];
    result[0] = A_5_0[0][0];
    result[1] = A_5_0[1][0];
    result[2] = A_5_0[2][0] & 0xFFFF0000;

    u32 bit_diff = 0;
#ifdef TEST_PRE
    bit_diff += __builtin_popcount(result[0] ^ 0xcc2dc21a);
    bit_diff += __builtin_popcount(result[1] ^ 0xe7d743bf);
    bit_diff += __builtin_popcount(result[2] ^ 0x02fe0000);
#else
    bit_diff += __builtin_popcount(result[0] ^ 0xAE5868A7);
    bit_diff += __builtin_popcount(result[1] ^ 0x27A98747);
    bit_diff += __builtin_popcount(result[2] ^ 0xFF440000);
#endif
    printHash(result);
    *minDiff = bit_diff;
    return bit_diff <= DIFF_TOLERANCE;
#undef u32
#undef _ROR32
}


void cpu_DumpKeccakResult(const uint32_t A[5][5], uint32_t hash[3]) {
#define _ROR32(x, a) ((((x)>>(a))|((x)<<(32-(a)))))
#define u32 uint32_t
    u32 A_1_0[5][5];
    for (uint32_t i = 0; i < 5; i++)
        memcpy(A_1_0[i], A[i], sizeof(uint32_t) * 5);

    u32 A_1_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_1_1[0][j] = A_1_0[0][j] ^ A_1_0[4][0] ^ A_1_0[4][1] ^ A_1_0[4][2] ^ A_1_0[4][3] ^ A_1_0[4][4]
                      ^ _ROR32((A_1_0[1][0] ^ A_1_0[1][1] ^ A_1_0[1][2] ^ A_1_0[1][3] ^ A_1_0[1][4]), 1);
        A_1_1[1][j] = A_1_0[1][j] ^ A_1_0[0][0] ^ A_1_0[0][1] ^ A_1_0[0][2] ^ A_1_0[0][3] ^ A_1_0[0][4]
                      ^ _ROR32((A_1_0[2][0] ^ A_1_0[2][1] ^ A_1_0[2][2] ^ A_1_0[2][3] ^ A_1_0[2][4]), 1);
        A_1_1[2][j] = A_1_0[2][j] ^ A_1_0[1][0] ^ A_1_0[1][1] ^ A_1_0[1][2] ^ A_1_0[1][3] ^ A_1_0[1][4]
                      ^ _ROR32((A_1_0[3][0] ^ A_1_0[3][1] ^ A_1_0[3][2] ^ A_1_0[3][3] ^ A_1_0[3][4]), 1);
        A_1_1[3][j] = A_1_0[3][j] ^ A_1_0[2][0] ^ A_1_0[2][1] ^ A_1_0[2][2] ^ A_1_0[2][3] ^ A_1_0[2][4]
                      ^ _ROR32((A_1_0[4][0] ^ A_1_0[4][1] ^ A_1_0[4][2] ^ A_1_0[4][3] ^ A_1_0[4][4]), 1);
        A_1_1[4][j] = A_1_0[4][j] ^ A_1_0[3][0] ^ A_1_0[3][1] ^ A_1_0[3][2] ^ A_1_0[3][3] ^ A_1_0[3][4]
                      ^ _ROR32((A_1_0[0][0] ^ A_1_0[0][1] ^ A_1_0[0][2] ^ A_1_0[0][3] ^ A_1_0[0][4]), 1);
    }

    u32 A_1_2[5][5];
    A_1_2[0][0] = A_1_1[0][0];
    A_1_2[0][1] = _ROR32(A_1_1[3][0], 28);
    A_1_2[0][2] = _ROR32(A_1_1[1][0], 1);
    A_1_2[0][3] = _ROR32(A_1_1[4][0], 27);
    A_1_2[0][4] = _ROR32(A_1_1[2][0], 30);
    A_1_2[1][0] = _ROR32(A_1_1[1][1], 12);
    A_1_2[1][1] = _ROR32(A_1_1[4][1], 20);
    A_1_2[1][2] = _ROR32(A_1_1[2][1], 6);
    A_1_2[1][3] = _ROR32(A_1_1[0][1], 4);
    A_1_2[1][4] = _ROR32(A_1_1[3][1], 23);
    A_1_2[2][0] = _ROR32(A_1_1[2][2], 11);
    A_1_2[2][1] = _ROR32(A_1_1[0][2], 3);
    A_1_2[2][2] = _ROR32(A_1_1[3][2], 25);
    A_1_2[2][3] = _ROR32(A_1_1[1][2], 10);
    A_1_2[2][4] = _ROR32(A_1_1[4][2], 7);
    A_1_2[3][0] = _ROR32(A_1_1[3][3], 21);
    A_1_2[3][1] = _ROR32(A_1_1[1][3], 13);
    A_1_2[3][2] = _ROR32(A_1_1[4][3], 8);
    A_1_2[3][3] = _ROR32(A_1_1[2][3], 15);
    A_1_2[3][4] = _ROR32(A_1_1[0][3], 9);
    A_1_2[4][0] = _ROR32(A_1_1[4][4], 14);
    A_1_2[4][1] = _ROR32(A_1_1[2][4], 29);
    A_1_2[4][2] = _ROR32(A_1_1[0][4], 18);
    A_1_2[4][3] = _ROR32(A_1_1[3][4], 24);
    A_1_2[4][4] = _ROR32(A_1_1[1][4], 2);

    u32 A_2_0[5][5];
    A_2_0[0][0] = A_1_2[0][0] ^ RC1;
    A_2_0[1][0] = 0xFFFFFFFF;
    A_2_0[2][0] = A_1_2[0][0] ^ A_1_2[2][0];
    A_2_0[3][0] = 0;
    A_2_0[4][0] = 0xFFFFFFFF;
    A_2_0[0][1] = A_1_2[0][1];
    A_2_0[1][1] = 0xFFFFFFFF;
    A_2_0[2][1] = A_1_2[0][1] ^ A_1_2[2][1];
    A_2_0[3][1] = 0;
    A_2_0[4][1] = 0xFFFFFFFF;
    A_2_0[0][2] = A_1_2[0][2];
    A_2_0[1][2] = 0xFFFFFFFF;
    A_2_0[2][2] = A_1_2[0][2] ^ A_1_2[2][2];
    A_2_0[3][2] = 0;
    A_2_0[4][2] = 0xFFFFFFFF;
    A_2_0[0][3] = A_1_2[0][3];
    A_2_0[1][3] = 0xFFFFFFFF;
    A_2_0[2][3] = A_1_2[0][3] ^ A_1_2[2][3];
    A_2_0[3][3] = 0;
    A_2_0[4][3] = 0xFFFFFFFF;
    A_2_0[0][4] = A_1_2[0][4];
    A_2_0[1][4] = 0xFFFFFFFF;
    A_2_0[2][4] = A_1_2[0][4] ^ A_1_2[2][4];
    A_2_0[3][4] = 0;
    A_2_0[4][4] = 0xFFFFFFFF;

    u32 A_2_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_2_1[0][j] = A_2_0[0][j];
        A_2_1[1][j] = 0xFFFFFFFF ^ ALPHA ^ _ROR32(BETA, 1);
        A_2_1[2][j] = A_2_0[2][j] ^ 0xFFFFFFFF;
        A_2_1[3][j] = BETA ^ 0xFFFFFFFF;
        A_2_1[4][j] = 0xFFFFFFFF ^ _ROR32(ALPHA, 1);
    }

    u32 A_2_2[5][5];
    A_2_2[0][0] = A_2_1[0][0];
    A_2_2[0][1] = _ROR32(A_2_1[3][0], 28);
    A_2_2[0][2] = _ROR32(A_2_1[1][0], 1);
    A_2_2[0][3] = _ROR32(A_2_1[4][0], 27);
    A_2_2[0][4] = _ROR32(A_2_1[2][0], 30);
    A_2_2[1][0] = _ROR32(A_2_1[1][1], 12);
    A_2_2[1][1] = _ROR32(A_2_1[4][1], 20);
    A_2_2[1][2] = _ROR32(A_2_1[2][1], 6);
    A_2_2[1][3] = _ROR32(A_2_1[0][1], 4);
    A_2_2[1][4] = _ROR32(A_2_1[3][1], 23);
    A_2_2[2][0] = _ROR32(A_2_1[2][2], 11);
    A_2_2[2][1] = _ROR32(A_2_1[0][2], 3);
    A_2_2[2][2] = _ROR32(A_2_1[3][2], 25);
    A_2_2[2][3] = _ROR32(A_2_1[1][2], 10);
    A_2_2[2][4] = _ROR32(A_2_1[4][2], 7);
    A_2_2[3][0] = _ROR32(A_2_1[3][3], 21);
    A_2_2[3][1] = _ROR32(A_2_1[1][3], 13);
    A_2_2[3][2] = _ROR32(A_2_1[4][3], 8);
    A_2_2[3][3] = _ROR32(A_2_1[2][3], 15);
    A_2_2[3][4] = _ROR32(A_2_1[0][3], 9);
    A_2_2[4][0] = _ROR32(A_2_1[4][4], 14);
    A_2_2[4][1] = _ROR32(A_2_1[2][4], 29);
    A_2_2[4][2] = _ROR32(A_2_1[0][4], 18);
    A_2_2[4][3] = _ROR32(A_2_1[3][4], 24);
    A_2_2[4][4] = _ROR32(A_2_1[1][4], 2);

    //3rd round
    u32 A_3_0[5][5];
    for (int j = 0; j < 5; j++) {
        A_3_0[0][j] = A_2_2[0][j] ^ ((A_2_2[1][j] ^ 0xFFFFFFFF) & A_2_2[2][j]);
        A_3_0[1][j] = A_2_2[1][j] ^ ((A_2_2[2][j] ^ 0xFFFFFFFF) & A_2_2[3][j]);
        A_3_0[2][j] = A_2_2[2][j] ^ ((A_2_2[3][j] ^ 0xFFFFFFFF) & A_2_2[4][j]);
        A_3_0[3][j] = A_2_2[3][j] ^ ((A_2_2[4][j] ^ 0xFFFFFFFF) & A_2_2[0][j]);
        A_3_0[4][j] = A_2_2[4][j] ^ ((A_2_2[0][j] ^ 0xFFFFFFFF) & A_2_2[1][j]);
    }
    A_3_0[0][0] ^= RC2;

    u32 A_3_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_3_1[0][j] = A_3_0[0][j] ^ A_3_0[4][0] ^ A_3_0[4][1] ^ A_3_0[4][2] ^ A_3_0[4][3] ^ A_3_0[4][4]
                      ^ _ROR32((A_3_0[1][0] ^ A_3_0[1][1] ^ A_3_0[1][2] ^ A_3_0[1][3] ^ A_3_0[1][4]), 1);
        A_3_1[1][j] = A_3_0[1][j] ^ A_3_0[0][0] ^ A_3_0[0][1] ^ A_3_0[0][2] ^ A_3_0[0][3] ^ A_3_0[0][4]
                      ^ _ROR32((A_3_0[2][0] ^ A_3_0[2][1] ^ A_3_0[2][2] ^ A_3_0[2][3] ^ A_3_0[2][4]), 1);
        A_3_1[2][j] = A_3_0[2][j] ^ A_3_0[1][0] ^ A_3_0[1][1] ^ A_3_0[1][2] ^ A_3_0[1][3] ^ A_3_0[1][4]
                      ^ _ROR32((A_3_0[3][0] ^ A_3_0[3][1] ^ A_3_0[3][2] ^ A_3_0[3][3] ^ A_3_0[3][4]), 1);
        A_3_1[3][j] = A_3_0[3][j] ^ A_3_0[2][0] ^ A_3_0[2][1] ^ A_3_0[2][2] ^ A_3_0[2][3] ^ A_3_0[2][4]
                      ^ _ROR32((A_3_0[4][0] ^ A_3_0[4][1] ^ A_3_0[4][2] ^ A_3_0[4][3] ^ A_3_0[4][4]), 1);
        A_3_1[4][j] = A_3_0[4][j] ^ A_3_0[3][0] ^ A_3_0[3][1] ^ A_3_0[3][2] ^ A_3_0[3][3] ^ A_3_0[3][4]
                      ^ _ROR32((A_3_0[0][0] ^ A_3_0[0][1] ^ A_3_0[0][2] ^ A_3_0[0][3] ^ A_3_0[0][4]), 1);
    }

    u32 A_3_2[5][5];
    A_3_2[0][0] = A_3_1[0][0];
    A_3_2[0][1] = _ROR32(A_3_1[3][0], 28);
    A_3_2[0][2] = _ROR32(A_3_1[1][0], 1);
    A_3_2[0][3] = _ROR32(A_3_1[4][0], 27);
    A_3_2[0][4] = _ROR32(A_3_1[2][0], 30);
    A_3_2[1][0] = _ROR32(A_3_1[1][1], 12);
    A_3_2[1][1] = _ROR32(A_3_1[4][1], 20);
    A_3_2[1][2] = _ROR32(A_3_1[2][1], 6);
    A_3_2[1][3] = _ROR32(A_3_1[0][1], 4);
    A_3_2[1][4] = _ROR32(A_3_1[3][1], 23);
    A_3_2[2][0] = _ROR32(A_3_1[2][2], 11);
    A_3_2[2][1] = _ROR32(A_3_1[0][2], 3);
    A_3_2[2][2] = _ROR32(A_3_1[3][2], 25);
    A_3_2[2][3] = _ROR32(A_3_1[1][2], 10);
    A_3_2[2][4] = _ROR32(A_3_1[4][2], 7);
    A_3_2[3][0] = _ROR32(A_3_1[3][3], 21);
    A_3_2[3][1] = _ROR32(A_3_1[1][3], 13);
    A_3_2[3][2] = _ROR32(A_3_1[4][3], 8);
    A_3_2[3][3] = _ROR32(A_3_1[2][3], 15);
    A_3_2[3][4] = _ROR32(A_3_1[0][3], 9);
    A_3_2[4][0] = _ROR32(A_3_1[4][4], 14);
    A_3_2[4][1] = _ROR32(A_3_1[2][4], 29);
    A_3_2[4][2] = _ROR32(A_3_1[0][4], 18);
    A_3_2[4][3] = _ROR32(A_3_1[3][4], 24);
    A_3_2[4][4] = _ROR32(A_3_1[1][4], 2);

    //4th round
    u32 A_4_0[5][5];
    for (int j = 0; j < 5; j++) {
        A_4_0[0][j] = A_3_2[0][j] ^ ((A_3_2[1][j] ^ 0xFFFFFFFF) & A_3_2[2][j]);
        A_4_0[1][j] = A_3_2[1][j] ^ ((A_3_2[2][j] ^ 0xFFFFFFFF) & A_3_2[3][j]);
        A_4_0[2][j] = A_3_2[2][j] ^ ((A_3_2[3][j] ^ 0xFFFFFFFF) & A_3_2[4][j]);
        A_4_0[3][j] = A_3_2[3][j] ^ ((A_3_2[4][j] ^ 0xFFFFFFFF) & A_3_2[0][j]);
        A_4_0[4][j] = A_3_2[4][j] ^ ((A_3_2[0][j] ^ 0xFFFFFFFF) & A_3_2[1][j]);
    }
    A_4_0[0][0] ^= RC3;

    u32 A_4_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_4_1[0][j] = A_4_0[0][j] ^ A_4_0[4][0] ^ A_4_0[4][1] ^ A_4_0[4][2] ^ A_4_0[4][3] ^ A_4_0[4][4]
                      ^ _ROR32((A_4_0[1][0] ^ A_4_0[1][1] ^ A_4_0[1][2] ^ A_4_0[1][3] ^ A_4_0[1][4]), 1);
        A_4_1[1][j] = A_4_0[1][j] ^ A_4_0[0][0] ^ A_4_0[0][1] ^ A_4_0[0][2] ^ A_4_0[0][3] ^ A_4_0[0][4]
                      ^ _ROR32((A_4_0[2][0] ^ A_4_0[2][1] ^ A_4_0[2][2] ^ A_4_0[2][3] ^ A_4_0[2][4]), 1);
        A_4_1[2][j] = A_4_0[2][j] ^ A_4_0[1][0] ^ A_4_0[1][1] ^ A_4_0[1][2] ^ A_4_0[1][3] ^ A_4_0[1][4]
                      ^ _ROR32((A_4_0[3][0] ^ A_4_0[3][1] ^ A_4_0[3][2] ^ A_4_0[3][3] ^ A_4_0[3][4]), 1);
        A_4_1[3][j] = A_4_0[3][j] ^ A_4_0[2][0] ^ A_4_0[2][1] ^ A_4_0[2][2] ^ A_4_0[2][3] ^ A_4_0[2][4]
                      ^ _ROR32((A_4_0[4][0] ^ A_4_0[4][1] ^ A_4_0[4][2] ^ A_4_0[4][3] ^ A_4_0[4][4]), 1);
        A_4_1[4][j] = A_4_0[4][j] ^ A_4_0[3][0] ^ A_4_0[3][1] ^ A_4_0[3][2] ^ A_4_0[3][3] ^ A_4_0[3][4]
                      ^ _ROR32((A_4_0[0][0] ^ A_4_0[0][1] ^ A_4_0[0][2] ^ A_4_0[0][3] ^ A_4_0[0][4]), 1);
    }

    u32 A_4_2[5][5];
    A_4_2[0][0] = A_4_1[0][0];
    A_4_2[0][1] = _ROR32(A_4_1[3][0], 28);
    A_4_2[0][2] = _ROR32(A_4_1[1][0], 1);
    A_4_2[0][3] = _ROR32(A_4_1[4][0], 27);
    A_4_2[0][4] = _ROR32(A_4_1[2][0], 30);
    A_4_2[1][0] = _ROR32(A_4_1[1][1], 12);
    A_4_2[1][1] = _ROR32(A_4_1[4][1], 20);
    A_4_2[1][2] = _ROR32(A_4_1[2][1], 6);
    A_4_2[1][3] = _ROR32(A_4_1[0][1], 4);
    A_4_2[1][4] = _ROR32(A_4_1[3][1], 23);
    A_4_2[2][0] = _ROR32(A_4_1[2][2], 11);
    A_4_2[2][1] = _ROR32(A_4_1[0][2], 3);
    A_4_2[2][2] = _ROR32(A_4_1[3][2], 25);
    A_4_2[2][3] = _ROR32(A_4_1[1][2], 10);
    A_4_2[2][4] = _ROR32(A_4_1[4][2], 7);
    A_4_2[3][0] = _ROR32(A_4_1[3][3], 21);
    A_4_2[3][1] = _ROR32(A_4_1[1][3], 13);
    A_4_2[3][2] = _ROR32(A_4_1[4][3], 8);
    A_4_2[3][3] = _ROR32(A_4_1[2][3], 15);
    A_4_2[3][4] = _ROR32(A_4_1[0][3], 9);
    A_4_2[4][0] = _ROR32(A_4_1[4][4], 14);
    A_4_2[4][1] = _ROR32(A_4_1[2][4], 29);
    A_4_2[4][2] = _ROR32(A_4_1[0][4], 18);
    A_4_2[4][3] = _ROR32(A_4_1[3][4], 24);
    A_4_2[4][4] = _ROR32(A_4_1[1][4], 2);


    u32 A_5_0[5][5];
    for (int j = 0; j < 5; j++) {
        A_5_0[0][j] = A_4_2[0][j] ^ ((A_4_2[1][j] ^ 0xFFFFFFFF) & A_4_2[2][j]);
        A_5_0[1][j] = A_4_2[1][j] ^ ((A_4_2[2][j] ^ 0xFFFFFFFF) & A_4_2[3][j]);
        A_5_0[2][j] = A_4_2[2][j] ^ ((A_4_2[3][j] ^ 0xFFFFFFFF) & A_4_2[4][j]);
        A_5_0[3][j] = A_4_2[3][j] ^ ((A_4_2[4][j] ^ 0xFFFFFFFF) & A_4_2[0][j]);
        A_5_0[4][j] = A_4_2[4][j] ^ ((A_4_2[0][j] ^ 0xFFFFFFFF) & A_4_2[1][j]);
    }
    A_5_0[0][0] ^= RC4;

    hash[0] = A_5_0[0][0];
    hash[1] = A_5_0[1][0];
    hash[2] = A_5_0[2][0] & 0xFFFF0000;
#undef u32
#undef _ROR32
}

void cpu_VerifyRound2(const uint32_t A[5][5]) {
#define _ROR32(x, a) ((((x)>>(a))|((x)<<(32-(a)))))
#define u32 uint32_t
    u32 A_1_0[5][5];
    for (uint32_t i = 0; i < 5; i++)
        memcpy(A_1_0[i], A[i], sizeof(uint32_t) * 5);

    u32 A_1_1[5][5];
    for (int j = 0; j < 5; j++) {
        A_1_1[0][j] = A_1_0[0][j] ^ A_1_0[4][0] ^ A_1_0[4][1] ^ A_1_0[4][2] ^ A_1_0[4][3] ^ A_1_0[4][4]
                      ^ _ROR32((A_1_0[1][0] ^ A_1_0[1][1] ^ A_1_0[1][2] ^ A_1_0[1][3] ^ A_1_0[1][4]), 1);
        A_1_1[1][j] = A_1_0[1][j] ^ A_1_0[0][0] ^ A_1_0[0][1] ^ A_1_0[0][2] ^ A_1_0[0][3] ^ A_1_0[0][4]
                      ^ _ROR32((A_1_0[2][0] ^ A_1_0[2][1] ^ A_1_0[2][2] ^ A_1_0[2][3] ^ A_1_0[2][4]), 1);
        A_1_1[2][j] = A_1_0[2][j] ^ A_1_0[1][0] ^ A_1_0[1][1] ^ A_1_0[1][2] ^ A_1_0[1][3] ^ A_1_0[1][4]
                      ^ _ROR32((A_1_0[3][0] ^ A_1_0[3][1] ^ A_1_0[3][2] ^ A_1_0[3][3] ^ A_1_0[3][4]), 1);
        A_1_1[3][j] = A_1_0[3][j] ^ A_1_0[2][0] ^ A_1_0[2][1] ^ A_1_0[2][2] ^ A_1_0[2][3] ^ A_1_0[2][4]
                      ^ _ROR32((A_1_0[4][0] ^ A_1_0[4][1] ^ A_1_0[4][2] ^ A_1_0[4][3] ^ A_1_0[4][4]), 1);
        A_1_1[4][j] = A_1_0[4][j] ^ A_1_0[3][0] ^ A_1_0[3][1] ^ A_1_0[3][2] ^ A_1_0[3][3] ^ A_1_0[3][4]
                      ^ _ROR32((A_1_0[0][0] ^ A_1_0[0][1] ^ A_1_0[0][2] ^ A_1_0[0][3] ^ A_1_0[0][4]), 1);
    }

    u32 A_1_2[5][5];
    A_1_2[0][0] = A_1_1[0][0];
    A_1_2[0][1] = _ROR32(A_1_1[3][0], 28);
    A_1_2[0][2] = _ROR32(A_1_1[1][0], 1);
    A_1_2[0][3] = _ROR32(A_1_1[4][0], 27);
    A_1_2[0][4] = _ROR32(A_1_1[2][0], 30);
    A_1_2[1][0] = _ROR32(A_1_1[1][1], 12);
    A_1_2[1][1] = _ROR32(A_1_1[4][1], 20);
    A_1_2[1][2] = _ROR32(A_1_1[2][1], 6);
    A_1_2[1][3] = _ROR32(A_1_1[0][1], 4);
    A_1_2[1][4] = _ROR32(A_1_1[3][1], 23);
    A_1_2[2][0] = _ROR32(A_1_1[2][2], 11);
    A_1_2[2][1] = _ROR32(A_1_1[0][2], 3);
    A_1_2[2][2] = _ROR32(A_1_1[3][2], 25);
    A_1_2[2][3] = _ROR32(A_1_1[1][2], 10);
    A_1_2[2][4] = _ROR32(A_1_1[4][2], 7);
    A_1_2[3][0] = _ROR32(A_1_1[3][3], 21);
    A_1_2[3][1] = _ROR32(A_1_1[1][3], 13);
    A_1_2[3][2] = _ROR32(A_1_1[4][3], 8);
    A_1_2[3][3] = _ROR32(A_1_1[2][3], 15);
    A_1_2[3][4] = _ROR32(A_1_1[0][3], 9);
    A_1_2[4][0] = _ROR32(A_1_1[4][4], 14);
    A_1_2[4][1] = _ROR32(A_1_1[2][4], 29);
    A_1_2[4][2] = _ROR32(A_1_1[0][4], 18);
    A_1_2[4][3] = _ROR32(A_1_1[3][4], 24);
    A_1_2[4][4] = _ROR32(A_1_1[1][4], 2);

#ifndef TEST_PRE
    if (A_1_2[1][0] != 0xFFFFFFFF) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[1][1] != 0xFFFFFFFF) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[1][2] != 0xFFFFFFFF) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[1][3] != 0xFFFFFFFF) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[1][4] != 0xFFFFFFFF) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[3][0] != 0) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[3][1] != 0) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[3][2] != 0) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[3][3] != 0) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[3][4] != 0) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[4][0] != A_1_2[0][0]) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[4][1] != A_1_2[0][1]) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[4][2] != A_1_2[0][2]) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[4][3] != A_1_2[0][3]) EXIT_WITH_MSG("round 1 constraint error\n");
    if (A_1_2[4][4] != A_1_2[0][4]) EXIT_WITH_MSG("round 1 constraint error\n");
#endif
    u32 A_2_0[5][5];
    A_2_0[0][0] = A_1_2[0][0] ^ RC1;
    A_2_0[1][0] = 0xFFFFFFFF;
    A_2_0[2][0] = A_1_2[0][0] ^ A_1_2[2][0];
    A_2_0[3][0] = 0;
    A_2_0[4][0] = 0xFFFFFFFF;
    A_2_0[0][1] = A_1_2[0][1];
    A_2_0[1][1] = 0xFFFFFFFF;
    A_2_0[2][1] = A_1_2[0][1] ^ A_1_2[2][1];
    A_2_0[3][1] = 0;
    A_2_0[4][1] = 0xFFFFFFFF;
    A_2_0[0][2] = A_1_2[0][2];
    A_2_0[1][2] = 0xFFFFFFFF;
    A_2_0[2][2] = A_1_2[0][2] ^ A_1_2[2][2];
    A_2_0[3][2] = 0;
    A_2_0[4][2] = 0xFFFFFFFF;
    A_2_0[0][3] = A_1_2[0][3];
    A_2_0[1][3] = 0xFFFFFFFF;
    A_2_0[2][3] = A_1_2[0][3] ^ A_1_2[2][3];
    A_2_0[3][3] = 0;
    A_2_0[4][3] = 0xFFFFFFFF;
    A_2_0[0][4] = A_1_2[0][4];
    A_2_0[1][4] = 0xFFFFFFFF;
    A_2_0[2][4] = A_1_2[0][4] ^ A_1_2[2][4];
    A_2_0[3][4] = 0;
    A_2_0[4][4] = 0xFFFFFFFF;

#ifndef TEST_PRE
    if ((A_2_0[0][0] ^ A_2_0[0][1] ^ A_2_0[0][2] ^ A_2_0[0][3] ^ A_2_0[0][4]) != ALPHA)
        EXIT_WITH_MSG("round 2 constraint error\n");
    if ((A_2_0[2][0] ^ A_2_0[2][1] ^ A_2_0[2][2] ^ A_2_0[2][3] ^ A_2_0[2][4]) != BETA)
        EXIT_WITH_MSG("round 2 constraint error\n");
#endif

#undef u32
#undef _ROR32
}
