//
// Created by vitowu on 3/13/20.
//

#include "../include/keccak.h"

// parameters relates to Keccak analysis
#define ALPHA 0x44E72
#define BETA 0xBAA20F
#define DIFF_TOLERANCE 0
#define RC1 0x80000000
#define RC2 0x41010000
#define RC3 0x51010000
#define RC4 0x00010001

__device__ __forceinline__
bool
gpu_VerifyKeccakResult(uint32_t A[25]) {
    uint32_t tmpA[25];
    uint32_t C[5];
    uint32_t D[5];

    // Round 1 start here
    // Unrolled THETA operation
    C[0] = A[5];
    C[0] = C[0] ^ A[6];
    C[0] = C[0] ^ A[7];
    C[0] = C[0] ^ A[8];
    C[0] = C[0] ^ A[9];
    D[0] = ROR32(C[0], 1);

    C[1] = A[10];
    C[1] = C[1] ^ A[11];
    C[1] = C[1] ^ A[12];
    C[1] = C[1] ^ A[13];
    C[1] = C[1] ^ A[14];
    D[1] = ROR32(C[1], 1);

    C[2] = A[15];
    C[2] = C[2] ^ A[16];
    C[2] = C[2] ^ A[17];
    C[2] = C[2] ^ A[18];
    C[2] = C[2] ^ A[19];
    D[2] = ROR32(C[2], 1);

    C[3] = A[20];
    C[3] = C[3] ^ A[21];
    C[3] = C[3] ^ A[22];
    C[3] = C[3] ^ A[23];
    C[3] = C[3] ^ A[24];
    D[3] = ROR32(C[3], 1);

    C[4] = A[0];
    C[4] = C[4] ^ A[1];
    C[4] = C[4] ^ A[2];
    C[4] = C[4] ^ A[3];
    C[4] = C[4] ^ A[4];
    D[4] = ROR32(C[4], 1);

    A[0] = A[0] ^ C[3] ^ D[0];
    A[1] = A[1] ^ C[3] ^ D[0];
    A[2] = A[2] ^ C[3] ^ D[0];
    A[3] = A[3] ^ C[3] ^ D[0];
    A[4] = A[4] ^ C[3] ^ D[0];

    A[5] = A[5] ^ C[4] ^ D[1];
    A[6] = A[6] ^ C[4] ^ D[1];
    A[7] = A[7] ^ C[4] ^ D[1];
    A[8] = A[8] ^ C[4] ^ D[1];
    A[9] = A[9] ^ C[4] ^ D[1];

    A[10] = A[10] ^ C[0] ^ D[2];
    A[11] = A[11] ^ C[0] ^ D[2];
    A[12] = A[12] ^ C[0] ^ D[2];
    A[13] = A[13] ^ C[0] ^ D[2];
    A[14] = A[14] ^ C[0] ^ D[2];

    A[15] = A[15] ^ C[1] ^ D[3];
    A[16] = A[16] ^ C[1] ^ D[3];
    A[17] = A[17] ^ C[1] ^ D[3];
    A[18] = A[18] ^ C[1] ^ D[3];
    A[19] = A[19] ^ C[1] ^ D[3];

    A[20] = A[20] ^ C[2] ^ D[4];
    A[21] = A[21] ^ C[2] ^ D[4];
    A[22] = A[22] ^ C[2] ^ D[4];
    A[23] = A[23] ^ C[2] ^ D[4];
    A[24] = A[24] ^ C[2] ^ D[4];

    // Unrolled RHO + PI operation
    // tmpA[0] = A[0];
    tmpA[1] = A[1];
    tmpA[2] = A[2];
    tmpA[3] = A[3];
    tmpA[4] = A[4];
    tmpA[5] = A[5];
    tmpA[6] = A[6];
    tmpA[7] = A[7];
    tmpA[8] = A[8];
    tmpA[9] = A[9];
    tmpA[10] = A[10];
    tmpA[11] = A[11];
    tmpA[12] = A[12];
    tmpA[13] = A[13];
    tmpA[14] = A[14];
    tmpA[15] = A[15];
    tmpA[16] = A[16];
    tmpA[17] = A[17];
    tmpA[18] = A[18];
    tmpA[19] = A[19];
    tmpA[20] = A[20];
    tmpA[21] = A[21];
    tmpA[22] = A[22];
    tmpA[23] = A[23];
    tmpA[24] = A[24];

    // A[0] = tmpA[0];
    A[1] = ROR32(tmpA[15], 28);
    A[2] = ROR32(tmpA[5], 1);
    A[3] = ROR32(tmpA[20], 27);
    A[4] = ROR32(tmpA[10], 30);

    A[5] = ROR32(tmpA[6], 12);
    A[6] = ROR32(tmpA[21], 20);
    A[7] = ROR32(tmpA[11], 6);
    A[8] = ROR32(tmpA[1], 4);
    A[9] = ROR32(tmpA[16], 23);

    A[10] = ROR32(tmpA[12], 11);
    A[11] = ROR32(tmpA[2], 3);
    A[12] = ROR32(tmpA[17], 25);
    A[13] = ROR32(tmpA[7], 10);
    A[14] = ROR32(tmpA[22], 7);

    A[15] = ROR32(tmpA[18], 21);
    A[16] = ROR32(tmpA[8], 13);
    A[17] = ROR32(tmpA[23], 8);
    A[18] = ROR32(tmpA[13], 15);
    A[19] = ROR32(tmpA[3], 9);

    A[20] = ROR32(tmpA[24], 14);
    A[21] = ROR32(tmpA[14], 29);
    A[22] = ROR32(tmpA[4], 18);
    A[23] = ROR32(tmpA[19], 24);
    A[24] = ROR32(tmpA[9], 2);

    // Unrolled CHI operation
    C[0] = A[0] ^ ((~A[5]) & (A[10]));
    C[1] = A[5] ^ ((~A[10]) & (A[15]));
    C[2] = A[10] ^ ((~A[15]) & (A[20]));
    C[3] = A[15] ^ ((~A[20]) & (A[0]));
    C[4] = A[20] ^ ((~A[0]) & (A[5]));
    A[0] = C[0];
    A[5] = C[1];
    A[10] = C[2];
    A[15] = C[3];
    A[20] = C[4];

    C[0] = A[1] ^ ((~A[6]) & (A[11]));
    C[1] = A[6] ^ ((~A[11]) & (A[16]));
    C[2] = A[11] ^ ((~A[16]) & (A[21]));
    C[3] = A[16] ^ ((~A[21]) & (A[1]));
    C[4] = A[21] ^ ((~A[1]) & (A[6]));
    A[1] = C[0];
    A[6] = C[1];
    A[11] = C[2];
    A[16] = C[3];
    A[21] = C[4];

    C[0] = A[2] ^ ((~A[7]) & (A[12]));
    C[1] = A[7] ^ ((~A[12]) & (A[17]));
    C[2] = A[12] ^ ((~A[17]) & (A[22]));
    C[3] = A[17] ^ ((~A[22]) & (A[2]));
    C[4] = A[22] ^ ((~A[2]) & (A[7]));
    A[2] = C[0];
    A[7] = C[1];
    A[12] = C[2];
    A[17] = C[3];
    A[22] = C[4];

    C[0] = A[3] ^ ((~A[8]) & (A[13]));
    C[1] = A[8] ^ ((~A[13]) & (A[18]));
    C[2] = A[13] ^ ((~A[18]) & (A[23]));
    C[3] = A[18] ^ ((~A[23]) & (A[3]));
    C[4] = A[23] ^ ((~A[3]) & (A[8]));
    A[3] = C[0];
    A[8] = C[1];
    A[13] = C[2];
    A[18] = C[3];
    A[23] = C[4];

    C[0] = A[4] ^ ((~A[9]) & (A[14]));
    C[1] = A[9] ^ ((~A[14]) & (A[19]));
    C[2] = A[14] ^ ((~A[19]) & (A[24]));
    C[3] = A[19] ^ ((~A[24]) & (A[4]));
    C[4] = A[24] ^ ((~A[4]) & (A[9]));
    A[4] = C[0];
    A[9] = C[1];
    A[14] = C[2];
    A[19] = C[3];
    A[24] = C[4];

    A[0] = A[0] ^ RC1;
    // --- Round 1 end here ---

    // --- Round 2 start here ---
    // Unrolled THETA operation
    C[0] = A[5];
    C[0] = C[0] ^ A[6];
    C[0] = C[0] ^ A[7];
    C[0] = C[0] ^ A[8];
    C[0] = C[0] ^ A[9];
    D[0] = ROR32(C[0], 1);

    C[1] = A[10];
    C[1] = C[1] ^ A[11];
    C[1] = C[1] ^ A[12];
    C[1] = C[1] ^ A[13];
    C[1] = C[1] ^ A[14];
    D[1] = ROR32(C[1], 1);

    C[2] = A[15];
    C[2] = C[2] ^ A[16];
    C[2] = C[2] ^ A[17];
    C[2] = C[2] ^ A[18];
    C[2] = C[2] ^ A[19];
    D[2] = ROR32(C[2], 1);

    C[3] = A[20];
    C[3] = C[3] ^ A[21];
    C[3] = C[3] ^ A[22];
    C[3] = C[3] ^ A[23];
    C[3] = C[3] ^ A[24];
    D[3] = ROR32(C[3], 1);

    C[4] = A[0];
    C[4] = C[4] ^ A[1];
    C[4] = C[4] ^ A[2];
    C[4] = C[4] ^ A[3];
    C[4] = C[4] ^ A[4];
    D[4] = ROR32(C[4], 1);

    A[0] = A[0] ^ C[3] ^ D[0];
    A[1] = A[1] ^ C[3] ^ D[0];
    A[2] = A[2] ^ C[3] ^ D[0];
    A[3] = A[3] ^ C[3] ^ D[0];
    A[4] = A[4] ^ C[3] ^ D[0];

    A[5] = A[5] ^ C[4] ^ D[1];
    A[6] = A[6] ^ C[4] ^ D[1];
    A[7] = A[7] ^ C[4] ^ D[1];
    A[8] = A[8] ^ C[4] ^ D[1];
    A[9] = A[9] ^ C[4] ^ D[1];

    A[10] = A[10] ^ C[0] ^ D[2];
    A[11] = A[11] ^ C[0] ^ D[2];
    A[12] = A[12] ^ C[0] ^ D[2];
    A[13] = A[13] ^ C[0] ^ D[2];
    A[14] = A[14] ^ C[0] ^ D[2];

    A[15] = A[15] ^ C[1] ^ D[3];
    A[16] = A[16] ^ C[1] ^ D[3];
    A[17] = A[17] ^ C[1] ^ D[3];
    A[18] = A[18] ^ C[1] ^ D[3];
    A[19] = A[19] ^ C[1] ^ D[3];

    A[20] = A[20] ^ C[2] ^ D[4];
    A[21] = A[21] ^ C[2] ^ D[4];
    A[22] = A[22] ^ C[2] ^ D[4];
    A[23] = A[23] ^ C[2] ^ D[4];
    A[24] = A[24] ^ C[2] ^ D[4];

    // Unrolled RHO + PI operation
    // tmpA[0] = A[0];
    tmpA[1] = A[1];
    tmpA[2] = A[2];
    tmpA[3] = A[3];
    tmpA[4] = A[4];
    tmpA[5] = A[5];
    tmpA[6] = A[6];
    tmpA[7] = A[7];
    tmpA[8] = A[8];
    tmpA[9] = A[9];
    tmpA[10] = A[10];
    tmpA[11] = A[11];
    tmpA[12] = A[12];
    tmpA[13] = A[13];
    tmpA[14] = A[14];
    tmpA[15] = A[15];
    tmpA[16] = A[16];
    tmpA[17] = A[17];
    tmpA[18] = A[18];
    tmpA[19] = A[19];
    tmpA[20] = A[20];
    tmpA[21] = A[21];
    tmpA[22] = A[22];
    tmpA[23] = A[23];
    tmpA[24] = A[24];

    // A[0] = tmpA[0];
    A[1] = ROR32(tmpA[15], 28);
    A[2] = ROR32(tmpA[5], 1);
    A[3] = ROR32(tmpA[20], 27);
    A[4] = ROR32(tmpA[10], 30);

    A[5] = ROR32(tmpA[6], 12);
    A[6] = ROR32(tmpA[21], 20);
    A[7] = ROR32(tmpA[11], 6);
    A[8] = ROR32(tmpA[1], 4);
    A[9] = ROR32(tmpA[16], 23);

    A[10] = ROR32(tmpA[12], 11);
    A[11] = ROR32(tmpA[2], 3);
    A[12] = ROR32(tmpA[17], 25);
    A[13] = ROR32(tmpA[7], 10);
    A[14] = ROR32(tmpA[22], 7);

    A[15] = ROR32(tmpA[18], 21);
    A[16] = ROR32(tmpA[8], 13);
    A[17] = ROR32(tmpA[23], 8);
    A[18] = ROR32(tmpA[13], 15);
    A[19] = ROR32(tmpA[3], 9);

    A[20] = ROR32(tmpA[24], 14);
    A[21] = ROR32(tmpA[14], 29);
    A[22] = ROR32(tmpA[4], 18);
    A[23] = ROR32(tmpA[19], 24);
    A[24] = ROR32(tmpA[9], 2);

    // Unrolled CHI operation
    C[0] = A[0] ^ ((~A[5]) & (A[10]));
    C[1] = A[5] ^ ((~A[10]) & (A[15]));
    C[2] = A[10] ^ ((~A[15]) & (A[20]));
    C[3] = A[15] ^ ((~A[20]) & (A[0]));
    C[4] = A[20] ^ ((~A[0]) & (A[5]));
    A[0] = C[0];
    A[5] = C[1];
    A[10] = C[2];
    A[15] = C[3];
    A[20] = C[4];

    C[0] = A[1] ^ ((~A[6]) & (A[11]));
    C[1] = A[6] ^ ((~A[11]) & (A[16]));
    C[2] = A[11] ^ ((~A[16]) & (A[21]));
    C[3] = A[16] ^ ((~A[21]) & (A[1]));
    C[4] = A[21] ^ ((~A[1]) & (A[6]));
    A[1] = C[0];
    A[6] = C[1];
    A[11] = C[2];
    A[16] = C[3];
    A[21] = C[4];

    C[0] = A[2] ^ ((~A[7]) & (A[12]));
    C[1] = A[7] ^ ((~A[12]) & (A[17]));
    C[2] = A[12] ^ ((~A[17]) & (A[22]));
    C[3] = A[17] ^ ((~A[22]) & (A[2]));
    C[4] = A[22] ^ ((~A[2]) & (A[7]));
    A[2] = C[0];
    A[7] = C[1];
    A[12] = C[2];
    A[17] = C[3];
    A[22] = C[4];

    C[0] = A[3] ^ ((~A[8]) & (A[13]));
    C[1] = A[8] ^ ((~A[13]) & (A[18]));
    C[2] = A[13] ^ ((~A[18]) & (A[23]));
    C[3] = A[18] ^ ((~A[23]) & (A[3]));
    C[4] = A[23] ^ ((~A[3]) & (A[8]));
    A[3] = C[0];
    A[8] = C[1];
    A[13] = C[2];
    A[18] = C[3];
    A[23] = C[4];

    C[0] = A[4] ^ ((~A[9]) & (A[14]));
    C[1] = A[9] ^ ((~A[14]) & (A[19]));
    C[2] = A[14] ^ ((~A[19]) & (A[24]));
    C[3] = A[19] ^ ((~A[24]) & (A[4]));
    C[4] = A[24] ^ ((~A[4]) & (A[9]));
    A[4] = C[0];
    A[9] = C[1];
    A[14] = C[2];
    A[19] = C[3];
    A[24] = C[4];

    A[0] = A[0] ^ RC2;
    // --- Round 2 end here ---

    // --- Round 3 start here ---
    // Unrolled THETA operation
    C[0] = A[5];
    C[0] = C[0] ^ A[6];
    C[0] = C[0] ^ A[7];
    C[0] = C[0] ^ A[8];
    C[0] = C[0] ^ A[9];
    D[0] = ROR32(C[0], 1);

    C[1] = A[10];
    C[1] = C[1] ^ A[11];
    C[1] = C[1] ^ A[12];
    C[1] = C[1] ^ A[13];
    C[1] = C[1] ^ A[14];
    D[1] = ROR32(C[1], 1);

    C[2] = A[15];
    C[2] = C[2] ^ A[16];
    C[2] = C[2] ^ A[17];
    C[2] = C[2] ^ A[18];
    C[2] = C[2] ^ A[19];
    D[2] = ROR32(C[2], 1);

    C[3] = A[20];
    C[3] = C[3] ^ A[21];
    C[3] = C[3] ^ A[22];
    C[3] = C[3] ^ A[23];
    C[3] = C[3] ^ A[24];
    D[3] = ROR32(C[3], 1);

    C[4] = A[0];
    C[4] = C[4] ^ A[1];
    C[4] = C[4] ^ A[2];
    C[4] = C[4] ^ A[3];
    C[4] = C[4] ^ A[4];
    D[4] = ROR32(C[4], 1);

    A[0] = A[0] ^ C[3] ^ D[0];
    A[1] = A[1] ^ C[3] ^ D[0];
    A[2] = A[2] ^ C[3] ^ D[0];
    A[3] = A[3] ^ C[3] ^ D[0];
    A[4] = A[4] ^ C[3] ^ D[0];

    A[5] = A[5] ^ C[4] ^ D[1];
    A[6] = A[6] ^ C[4] ^ D[1];
    A[7] = A[7] ^ C[4] ^ D[1];
    A[8] = A[8] ^ C[4] ^ D[1];
    A[9] = A[9] ^ C[4] ^ D[1];

    A[10] = A[10] ^ C[0] ^ D[2];
    A[11] = A[11] ^ C[0] ^ D[2];
    A[12] = A[12] ^ C[0] ^ D[2];
    A[13] = A[13] ^ C[0] ^ D[2];
    A[14] = A[14] ^ C[0] ^ D[2];

    A[15] = A[15] ^ C[1] ^ D[3];
    A[16] = A[16] ^ C[1] ^ D[3];
    A[17] = A[17] ^ C[1] ^ D[3];
    A[18] = A[18] ^ C[1] ^ D[3];
    A[19] = A[19] ^ C[1] ^ D[3];

    A[20] = A[20] ^ C[2] ^ D[4];
    A[21] = A[21] ^ C[2] ^ D[4];
    A[22] = A[22] ^ C[2] ^ D[4];
    A[23] = A[23] ^ C[2] ^ D[4];
    A[24] = A[24] ^ C[2] ^ D[4];

    // Unrolled RHO + PI operation
    // tmpA[0] = A[0];
    tmpA[1] = A[1];
    tmpA[2] = A[2];
    tmpA[3] = A[3];
    tmpA[4] = A[4];
    tmpA[5] = A[5];
    tmpA[6] = A[6];
    tmpA[7] = A[7];
    tmpA[8] = A[8];
    tmpA[9] = A[9];
    tmpA[10] = A[10];
    tmpA[11] = A[11];
    tmpA[12] = A[12];
    tmpA[13] = A[13];
    tmpA[14] = A[14];
    tmpA[15] = A[15];
    tmpA[16] = A[16];
    tmpA[17] = A[17];
    tmpA[18] = A[18];
    tmpA[19] = A[19];
    tmpA[20] = A[20];
    tmpA[21] = A[21];
    tmpA[22] = A[22];
    tmpA[23] = A[23];
    tmpA[24] = A[24];

    // A[0] = tmpA[0];
    A[1] = ROR32(tmpA[15], 28);
    A[2] = ROR32(tmpA[5], 1);
    A[3] = ROR32(tmpA[20], 27);
    A[4] = ROR32(tmpA[10], 30);

    A[5] = ROR32(tmpA[6], 12);
    A[6] = ROR32(tmpA[21], 20);
    A[7] = ROR32(tmpA[11], 6);
    A[8] = ROR32(tmpA[1], 4);
    A[9] = ROR32(tmpA[16], 23);

    A[10] = ROR32(tmpA[12], 11);
    A[11] = ROR32(tmpA[2], 3);
    A[12] = ROR32(tmpA[17], 25);
    A[13] = ROR32(tmpA[7], 10);
    A[14] = ROR32(tmpA[22], 7);

    A[15] = ROR32(tmpA[18], 21);
    A[16] = ROR32(tmpA[8], 13);
    A[17] = ROR32(tmpA[23], 8);
    A[18] = ROR32(tmpA[13], 15);
    A[19] = ROR32(tmpA[3], 9);

    A[20] = ROR32(tmpA[24], 14);
    A[21] = ROR32(tmpA[14], 29);
    A[22] = ROR32(tmpA[4], 18);
    A[23] = ROR32(tmpA[19], 24);
    A[24] = ROR32(tmpA[9], 2);

    // Unrolled CHI operation
    C[0] = A[0] ^ ((~A[5]) & (A[10]));
    C[1] = A[5] ^ ((~A[10]) & (A[15]));
    C[2] = A[10] ^ ((~A[15]) & (A[20]));
    C[3] = A[15] ^ ((~A[20]) & (A[0]));
    C[4] = A[20] ^ ((~A[0]) & (A[5]));
    A[0] = C[0];
    A[5] = C[1];
    A[10] = C[2];
    A[15] = C[3];
    A[20] = C[4];

    C[0] = A[1] ^ ((~A[6]) & (A[11]));
    C[1] = A[6] ^ ((~A[11]) & (A[16]));
    C[2] = A[11] ^ ((~A[16]) & (A[21]));
    C[3] = A[16] ^ ((~A[21]) & (A[1]));
    C[4] = A[21] ^ ((~A[1]) & (A[6]));
    A[1] = C[0];
    A[6] = C[1];
    A[11] = C[2];
    A[16] = C[3];
    A[21] = C[4];

    C[0] = A[2] ^ ((~A[7]) & (A[12]));
    C[1] = A[7] ^ ((~A[12]) & (A[17]));
    C[2] = A[12] ^ ((~A[17]) & (A[22]));
    C[3] = A[17] ^ ((~A[22]) & (A[2]));
    C[4] = A[22] ^ ((~A[2]) & (A[7]));
    A[2] = C[0];
    A[7] = C[1];
    A[12] = C[2];
    A[17] = C[3];
    A[22] = C[4];

    C[0] = A[3] ^ ((~A[8]) & (A[13]));
    C[1] = A[8] ^ ((~A[13]) & (A[18]));
    C[2] = A[13] ^ ((~A[18]) & (A[23]));
    C[3] = A[18] ^ ((~A[23]) & (A[3]));
    C[4] = A[23] ^ ((~A[3]) & (A[8]));
    A[3] = C[0];
    A[8] = C[1];
    A[13] = C[2];
    A[18] = C[3];
    A[23] = C[4];

    C[0] = A[4] ^ ((~A[9]) & (A[14]));
    C[1] = A[9] ^ ((~A[14]) & (A[19]));
    C[2] = A[14] ^ ((~A[19]) & (A[24]));
    C[3] = A[19] ^ ((~A[24]) & (A[4]));
    C[4] = A[24] ^ ((~A[4]) & (A[9]));
    A[4] = C[0];
    A[9] = C[1];
    A[14] = C[2];
    A[19] = C[3];
    A[24] = C[4];

    A[0] = A[0] ^ RC3;
    // --- Round 3 end here ---

    // --- Round 4 start here ---
    // Unrolled THETA operation
    C[0] = A[5];
    C[0] = C[0] ^ A[6];
    C[0] = C[0] ^ A[7];
    C[0] = C[0] ^ A[8];
    C[0] = C[0] ^ A[9];
    D[0] = ROR32(C[0], 1);

    C[1] = A[10];
    C[1] = C[1] ^ A[11];
    C[1] = C[1] ^ A[12];
    C[1] = C[1] ^ A[13];
    C[1] = C[1] ^ A[14];
    D[1] = ROR32(C[1], 1);

    C[2] = A[15];
    C[2] = C[2] ^ A[16];
    C[2] = C[2] ^ A[17];
    C[2] = C[2] ^ A[18];
    C[2] = C[2] ^ A[19];
    D[2] = ROR32(C[2], 1);

    C[3] = A[20];
    C[3] = C[3] ^ A[21];
    C[3] = C[3] ^ A[22];
    C[3] = C[3] ^ A[23];
    C[3] = C[3] ^ A[24];
    D[3] = ROR32(C[3], 1);

    C[4] = A[0];
    C[4] = C[4] ^ A[1];
    C[4] = C[4] ^ A[2];
    C[4] = C[4] ^ A[3];
    C[4] = C[4] ^ A[4];
    D[4] = ROR32(C[4], 1);

    A[0] = A[0] ^ C[3] ^ D[0];
    A[1] = A[1] ^ C[3] ^ D[0];
    A[2] = A[2] ^ C[3] ^ D[0];
    A[3] = A[3] ^ C[3] ^ D[0];
    A[4] = A[4] ^ C[3] ^ D[0];

    A[5] = A[5] ^ C[4] ^ D[1];
    A[6] = A[6] ^ C[4] ^ D[1];
    A[7] = A[7] ^ C[4] ^ D[1];
    A[8] = A[8] ^ C[4] ^ D[1];
    A[9] = A[9] ^ C[4] ^ D[1];

    A[10] = A[10] ^ C[0] ^ D[2];
    A[11] = A[11] ^ C[0] ^ D[2];
    A[12] = A[12] ^ C[0] ^ D[2];
    A[13] = A[13] ^ C[0] ^ D[2];
    A[14] = A[14] ^ C[0] ^ D[2];

    A[15] = A[15] ^ C[1] ^ D[3];
    A[16] = A[16] ^ C[1] ^ D[3];
    A[17] = A[17] ^ C[1] ^ D[3];
    A[18] = A[18] ^ C[1] ^ D[3];
    A[19] = A[19] ^ C[1] ^ D[3];

    A[20] = A[20] ^ C[2] ^ D[4];
    A[21] = A[21] ^ C[2] ^ D[4];
    A[22] = A[22] ^ C[2] ^ D[4];
    A[23] = A[23] ^ C[2] ^ D[4];
    A[24] = A[24] ^ C[2] ^ D[4];

    // Unrolled RHO + PI operation
    // tmpA[0] = A[0];
    tmpA[1] = A[1];
    tmpA[2] = A[2];
    tmpA[3] = A[3];
    tmpA[4] = A[4];
    tmpA[5] = A[5];
    tmpA[6] = A[6];
    tmpA[7] = A[7];
    tmpA[8] = A[8];
    tmpA[9] = A[9];
    tmpA[10] = A[10];
    tmpA[11] = A[11];
    tmpA[12] = A[12];
    tmpA[13] = A[13];
    tmpA[14] = A[14];
    tmpA[15] = A[15];
    tmpA[16] = A[16];
    tmpA[17] = A[17];
    tmpA[18] = A[18];
    tmpA[19] = A[19];
    tmpA[20] = A[20];
    tmpA[21] = A[21];
    tmpA[22] = A[22];
    tmpA[23] = A[23];
    tmpA[24] = A[24];

    // A[0] = tmpA[0];
    A[1] = ROR32(tmpA[15], 28);
    A[2] = ROR32(tmpA[5], 1);
    A[3] = ROR32(tmpA[20], 27);
    A[4] = ROR32(tmpA[10], 30);

    A[5] = ROR32(tmpA[6], 12);
    A[6] = ROR32(tmpA[21], 20);
    A[7] = ROR32(tmpA[11], 6);
    A[8] = ROR32(tmpA[1], 4);
    A[9] = ROR32(tmpA[16], 23);

    A[10] = ROR32(tmpA[12], 11);
    A[11] = ROR32(tmpA[2], 3);
    A[12] = ROR32(tmpA[17], 25);
    A[13] = ROR32(tmpA[7], 10);
    A[14] = ROR32(tmpA[22], 7);

    A[15] = ROR32(tmpA[18], 21);
    A[16] = ROR32(tmpA[8], 13);
    A[17] = ROR32(tmpA[23], 8);
    A[18] = ROR32(tmpA[13], 15);
    A[19] = ROR32(tmpA[3], 9);

    A[20] = ROR32(tmpA[24], 14);
    A[21] = ROR32(tmpA[14], 29);
    A[22] = ROR32(tmpA[4], 18);
    A[23] = ROR32(tmpA[19], 24);
    A[24] = ROR32(tmpA[9], 2);

    // Unrolled CHI operation
    C[0] = A[0] ^ ((~A[5]) & (A[10]));
    C[1] = A[5] ^ ((~A[10]) & (A[15]));
    C[2] = A[10] ^ ((~A[15]) & (A[20]));
    C[3] = A[15] ^ ((~A[20]) & (A[0]));
    C[4] = A[20] ^ ((~A[0]) & (A[5]));
    A[0] = C[0];
    A[5] = C[1];
    A[10] = C[2];
    A[15] = C[3];
    A[20] = C[4];

    C[0] = A[1] ^ ((~A[6]) & (A[11]));
    C[1] = A[6] ^ ((~A[11]) & (A[16]));
    C[2] = A[11] ^ ((~A[16]) & (A[21]));
    C[3] = A[16] ^ ((~A[21]) & (A[1]));
    C[4] = A[21] ^ ((~A[1]) & (A[6]));
    A[1] = C[0];
    A[6] = C[1];
    A[11] = C[2];
    A[16] = C[3];
    A[21] = C[4];

    C[0] = A[2] ^ ((~A[7]) & (A[12]));
    C[1] = A[7] ^ ((~A[12]) & (A[17]));
    C[2] = A[12] ^ ((~A[17]) & (A[22]));
    C[3] = A[17] ^ ((~A[22]) & (A[2]));
    C[4] = A[22] ^ ((~A[2]) & (A[7]));
    A[2] = C[0];
    A[7] = C[1];
    A[12] = C[2];
    A[17] = C[3];
    A[22] = C[4];

    C[0] = A[3] ^ ((~A[8]) & (A[13]));
    C[1] = A[8] ^ ((~A[13]) & (A[18]));
    C[2] = A[13] ^ ((~A[18]) & (A[23]));
    C[3] = A[18] ^ ((~A[23]) & (A[3]));
    C[4] = A[23] ^ ((~A[3]) & (A[8]));
    A[3] = C[0];
    A[8] = C[1];
    A[13] = C[2];
    A[18] = C[3];
    A[23] = C[4];

    C[0] = A[4] ^ ((~A[9]) & (A[14]));
    C[1] = A[9] ^ ((~A[14]) & (A[19]));
    C[2] = A[14] ^ ((~A[19]) & (A[24]));
    C[3] = A[19] ^ ((~A[24]) & (A[4]));
    C[4] = A[24] ^ ((~A[4]) & (A[9]));
    A[4] = C[0];
    A[9] = C[1];
    A[14] = C[2];
    A[19] = C[3];
    A[24] = C[4];

    A[0] = A[0] ^ RC4;
    // --- Round 4 end here ---

    // check correctness of hash
    uint32_t hash0 = A[0];
    uint32_t hash1 = A[5];
    uint32_t hash2 = A[10] & 0xFFFF0000;
    uint32_t bit_diff = 0;
    bit_diff += __popc(hash0 ^ 0x751A16E5);
    bit_diff += __popc(hash1 ^ 0xE495E1E2);
    bit_diff += __popc(hash2 ^ 0xFF220000);

    return bit_diff <= DIFF_TOLERANCE;
}
