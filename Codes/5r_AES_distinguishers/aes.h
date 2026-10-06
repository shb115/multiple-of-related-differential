#ifndef AES_H
#define AES_H
#include <stdint.h>
/* Standard AES-128. enc_rounds(out, in, rk, nr): nr-round AES (nr-1 full
 * rounds + 1 final round without MixColumns), rk = expanded key (44 words). */
void aes128_key_expansion(const uint8_t key[16], uint32_t rk[44]);
void aes_enc_rounds(uint8_t out[16], const uint8_t in[16], const uint32_t rk[44], int nr);
int aes_self_test(void);   /* returns 1 if FIPS-197 KAT passes */
#endif
