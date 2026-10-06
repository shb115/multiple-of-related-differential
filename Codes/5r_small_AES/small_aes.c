/*
 * small_aes.c
 * Small-AES SR(n,4,4,4) of Cid, Murphy and Robshaw (FSE 2005) with its key schedule.
 * The whitening key k_0 is the 64-bit master key and round r adds the round key
 * k_r of the small-AES key schedule (an earlier version of this artifact added
 * the master key in every round).
 *
 * Key schedule (r = c = 4, e = 4, field GF(2^4) with x^4 + x + 1):
 *   k_0 = master key
 *   k_i[., 0] = k_{i-1}[., 0] + S(k_{i-1}[(row+1) mod 4, 3]) + (kappa_i, 0, 0, 0)
 *   k_i[., j] = k_{i-1}[., j] + k_i[., j-1],  j = 1, 2, 3
 *   kappa_i = x^{i-1}: 1, 2, 4, 8, 3, 6, C, B, 5, A
 */

#include "small_aes.h"
#include <string.h>

static const uint8_t sbox[16] = {
    0x6, 0xB, 0x5, 0x4, 0x2, 0xE, 0x7, 0xA,
    0x9, 0xD, 0xF, 0xC, 0x3, 0x1, 0x0, 0x8
};
static const uint8_t rcon[11] = { 0x0, 0x1, 0x2, 0x4, 0x8, 0x3, 0x6, 0xC, 0xB, 0x5, 0xA };

static const uint16_t T0[16] = {
    0xC66A, 0x5BBE, 0xA55F, 0x844C, 0x4226, 0xFEE1, 0xE779, 0x7AAD,
    0x1998, 0x9DD4, 0xDFF2, 0xBCC7, 0x6335, 0x2113, 0x0000, 0x388B
};
static const uint16_t T1[16] = {
    0xAC66, 0xE5BB, 0xFA55, 0xC844, 0x6422, 0x1FEE, 0x9E77, 0xD7AA,
    0x8199, 0x49DD, 0x2DFF, 0x7BCC, 0x5633, 0x3211, 0x0000, 0xB388
};
static const uint16_t T2[16] = {
    0x6AC6, 0xBE5B, 0x5FA5, 0x4C84, 0x2642, 0xE1FE, 0x79E7, 0xAD7A,
    0x9819, 0xD49D, 0xF2DF, 0xC7BC, 0x3563, 0x1321, 0x0000, 0x8B38
};
static const uint16_t T3[16] = {
    0x66AC, 0xBBE5, 0x55FA, 0x44C8, 0x2264, 0xEE1F, 0x779E, 0xAAD7,
    0x9981, 0xDD49, 0xFF2D, 0xCC7B, 0x3356, 0x1132, 0x0000, 0x88B3
};

static void bytesToState(const uint8_t *bytes, uint8_t state[4][4]) {
    for (int j = 0; j < 4; j++) {
        state[0][j] = (bytes[j * 2] >> 4) & 0xF;
        state[1][j] =  bytes[j * 2]       & 0xF;
        state[2][j] = (bytes[j * 2 + 1] >> 4) & 0xF;
        state[3][j] =  bytes[j * 2 + 1]       & 0xF;
    }
}
static void stateToBytes(uint8_t state[4][4], uint8_t *bytes) {
    for (int j = 0; j < 4; j++) {
        bytes[j * 2]     = (state[0][j] << 4) | state[1][j];
        bytes[j * 2 + 1] = (state[2][j] << 4) | state[3][j];
    }
}
static void addRoundKey(uint8_t state[4][4], uint8_t key[4][4]) {
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            state[i][j] ^= key[i][j];
}
static void subNibbles(uint8_t state[4][4]) {
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            state[i][j] = sbox[state[i][j]];
}
static void shiftRows(uint8_t state[4][4]) {
    uint8_t t;
    t = state[1][0]; state[1][0] = state[1][1]; state[1][1] = state[1][2]; state[1][2] = state[1][3]; state[1][3] = t;
    t = state[2][0]; state[2][0] = state[2][2]; state[2][2] = t;
    t = state[2][1]; state[2][1] = state[2][3]; state[2][3] = t;
    t = state[3][0]; state[3][0] = state[3][3]; state[3][3] = state[3][2]; state[3][2] = state[3][1]; state[3][1] = t;
}

/* round keys k_0 .. k_10, cached for the last master key */
static uint8_t cached_mk[8];
static int cache_valid = 0;
static uint8_t RK[11][4][4];

static void expand_key(const uint8_t *masterKey) {
    if (cache_valid && memcmp(cached_mk, masterKey, 8) == 0) return;
    bytesToState(masterKey, RK[0]);
    for (int i = 1; i <= 10; i++) {
        uint8_t t[4];
        for (int row = 0; row < 4; row++) t[row] = sbox[RK[i - 1][(row + 1) & 3][3]];
        t[0] ^= rcon[i];
        for (int row = 0; row < 4; row++) RK[i][row][0] = RK[i - 1][row][0] ^ t[row];
        for (int j = 1; j < 4; j++)
            for (int row = 0; row < 4; row++) RK[i][row][j] = RK[i - 1][row][j] ^ RK[i][row][j - 1];
    }
    memcpy(cached_mk, masterKey, 8);
    cache_valid = 1;
}

void small_aes_enc_ttable(uint8_t *ciphertext, const uint8_t *plaintext,
                          const uint8_t *masterKey, int num_rounds) {
    uint8_t state[4][4];
    expand_key(masterKey);
    bytesToState(plaintext, state);
    addRoundKey(state, RK[0]);
    for (int r = 0; r < num_rounds - 1; r++) {
        uint16_t next_cols[4];
        next_cols[0] = T0[state[0][0]] ^ T1[state[1][1]] ^ T2[state[2][2]] ^ T3[state[3][3]];
        next_cols[1] = T0[state[0][1]] ^ T1[state[1][2]] ^ T2[state[2][3]] ^ T3[state[3][0]];
        next_cols[2] = T0[state[0][2]] ^ T1[state[1][3]] ^ T2[state[2][0]] ^ T3[state[3][1]];
        next_cols[3] = T0[state[0][3]] ^ T1[state[1][0]] ^ T2[state[2][1]] ^ T3[state[3][2]];
        for (int j = 0; j < 4; j++) {
            state[0][j] = (next_cols[j] >> 12) & 0xF;
            state[1][j] = (next_cols[j] >>  8) & 0xF;
            state[2][j] = (next_cols[j] >>  4) & 0xF;
            state[3][j] =  next_cols[j]        & 0xF;
        }
        addRoundKey(state, RK[r + 1]);
    }
    subNibbles(state);
    shiftRows(state);
    addRoundKey(state, RK[num_rounds]);
    stateToBytes(state, ciphertext);
}
