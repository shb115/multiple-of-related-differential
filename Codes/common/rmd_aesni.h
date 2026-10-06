#ifndef RMD_AESNI_H
#define RMD_AESNI_H

#include <stdint.h>
#include <wmmintrin.h>

static inline __m128i rmd_aes128_key_expansion(__m128i key, __m128i keygened) {
    keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3, 3, 3, 3));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    return _mm_xor_si128(key, keygened);
}

#define RMD_AES128_KEY_EXP(k, rcon) \
    rmd_aes128_key_expansion((k), _mm_aeskeygenassist_si128((k), (rcon)))

static inline void rmd_aes128_generate_round_keys(const uint8_t *mk, __m128i *rk) {
    rk[0]  = _mm_loadu_si128((const __m128i *)mk);
    rk[1]  = RMD_AES128_KEY_EXP(rk[0], 0x01);
    rk[2]  = RMD_AES128_KEY_EXP(rk[1], 0x02);
    rk[3]  = RMD_AES128_KEY_EXP(rk[2], 0x04);
    rk[4]  = RMD_AES128_KEY_EXP(rk[3], 0x08);
    rk[5]  = RMD_AES128_KEY_EXP(rk[4], 0x10);
    rk[6]  = RMD_AES128_KEY_EXP(rk[5], 0x20);
    rk[7]  = RMD_AES128_KEY_EXP(rk[6], 0x40);
    rk[8]  = RMD_AES128_KEY_EXP(rk[7], 0x80);
    rk[9]  = RMD_AES128_KEY_EXP(rk[8], 0x1B);
    rk[10] = RMD_AES128_KEY_EXP(rk[9], 0x36);
}

static inline void rmd_aes128_encrypt(const uint8_t *pt, uint8_t *ct,
                                      const __m128i *round_keys, int rounds) {
    __m128i m = _mm_loadu_si128((const __m128i *)pt);
    m = _mm_xor_si128(m, round_keys[0]);
    for (int i = 1; i < rounds; i++)
        m = _mm_aesenc_si128(m, round_keys[i]);
    m = _mm_aesenclast_si128(m, round_keys[rounds]);
    _mm_storeu_si128((__m128i *)ct, m);
}

static const int RMD_AES_INV_DIAG_INDICES[4][4] = {
    { 0, 13, 10,  7},
    { 4,  1, 14, 11},
    { 8,  5,  2, 15},
    {12,  9,  6,  3}
};

static inline uint32_t rmd_aes_inv_diag_val(const uint8_t *ct, int d) {
    uint32_t v = 0;
    v |= (uint32_t)ct[RMD_AES_INV_DIAG_INDICES[d][0]] << 24;
    v |= (uint32_t)ct[RMD_AES_INV_DIAG_INDICES[d][1]] << 16;
    v |= (uint32_t)ct[RMD_AES_INV_DIAG_INDICES[d][2]] <<  8;
    v |= (uint32_t)ct[RMD_AES_INV_DIAG_INDICES[d][3]];
    return v;
}

#endif /* RMD_AESNI_H */
