//
// Created by vitowu on 6/13/20.
//

#include <iostream>

typedef unsigned int u32;
typedef unsigned long long u64;
#define ROR32(x, a) (((x)>>(a))|((x)<<(32-(a))))
#define RC1 0x80000000
#define RC2 0x41010000
#define RC3 0x51010000
#define RC4 0x00010001

int display() {
    //1st round
    u32 A_1_0[5][5] = {
        {0x600c0000, 0x9fffffff, 0x600c0006, 0x60000000, 0x00000000},
        {0x30003006, 0xcffffff9, 0xf0003006, 0x30000006, 0x00000000},
        {0xa3000000, 0x3fffffff, 0xc3000000, 0xc0000000, 0x00000000},
        {0x00000000, 0x9fffffff, 0x00000000, 0x60000000, 0x00000000},
        {0x00000000, 0xfffffffc, 0x00000000, 0x00000003, 0x00000000}
    };

    A_1_0[0][0] = 0x00000000;
    A_1_0[0][1] = 0x00000000;
    A_1_0[0][2] = 0xffffffff;
    A_1_0[0][3] = 0xffffffff;
    A_1_0[0][4] = 0x00000000;
    A_1_0[1][0] = 0x00000000;
    A_1_0[1][1] = 0x00000000;
    A_1_0[1][2] = 0xffffffff;
    A_1_0[1][3] = 0xffffffff;
    A_1_0[1][4] = 0x00000000;
    A_1_0[2][0] = 0x00000000;
    A_1_0[2][1] = 0x00000000;
    A_1_0[2][2] = 0x00000000;
    A_1_0[2][3] = 0xffffffff;
    A_1_0[2][4] = 0x00000000;
    A_1_0[3][0] = 0xffffffff;
    A_1_0[3][1] = 0xffffffff;
    A_1_0[3][2] = 0xffffffff;
    A_1_0[3][3] = 0x00000000;
    A_1_0[3][4] = 0x00000000;
    A_1_0[4][0] = 0xffffffff;
    A_1_0[4][1] = 0x00000000;
    A_1_0[4][2] = 0xffffffff;
    A_1_0[4][3] = 0xffffffff;
    A_1_0[4][4] = 0x00000000;

    u32 A_1_1[5][5];
    for (int j = 0; j < 5; j++)
    {
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

    u32 A_1_2[5][5];
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

    //2nd round
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

    u32 Alpha = A_2_0[0][0] ^A_2_0[0][1] ^A_2_0[0][2] ^A_2_0[0][3] ^A_2_0[0][4];
    u32 Beta = A_2_0[2][0] ^A_2_0[2][1] ^A_2_0[2][2] ^A_2_0[2][3] ^A_2_0[2][4];

    printf("Alpha: 0x%08x\tBeta: 0x%08x\n", Alpha, Beta);
    return 0;
}

int main() { return display(); }