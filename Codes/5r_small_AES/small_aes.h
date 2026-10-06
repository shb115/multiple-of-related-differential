#ifndef SMALL_AES_H
#define SMALL_AES_H

#include <stdint.h>

/*
 * Small-AES encryption via T-tables.
 * 4-bit cells, 4x4 state (8-byte block), round keys from the small-AES
 * key schedule of Cid, Murphy and Robshaw (FSE 2005).
 *
 *  ciphertext [out] : 8-byte output buffer.
 *  plaintext  [in]  : 8-byte plaintext.
 *  masterKey  [in]  : 8-byte master key (whitening key k_0).
 *  num_rounds       : total number of rounds (1 to 10).
 */
void small_aes_enc_ttable(uint8_t *ciphertext,
                          const uint8_t *plaintext,
                          const uint8_t *masterKey,
                          int num_rounds);

#endif /* SMALL_AES_H */
