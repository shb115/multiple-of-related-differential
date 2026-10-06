#ifndef RMD_RANDOM_H
#define RMD_RANDOM_H

static inline unsigned int rmd_derive_trial_seed(unsigned int base, int t, int round) {
    unsigned int s = base * 0x9E3779B1u
                   + (unsigned int)t * 0xBF58476Du
                   + (unsigned int)round * 0x94D049BBu;
    s ^= s >> 16;
    s *= 0x7FEB352Du;
    s ^= s >> 15;
    return s | 1u;
}

static inline unsigned int rmd_derive_small_trial_seed(unsigned int base, int t) {
    unsigned int s = base * 0x9E3779B1u
                   + (unsigned int)t * 0xBF58476Du
                   + 0x94D049BBu;
    s ^= s >> 16;
    s *= 0x7FEB352Du;
    s ^= s >> 15;
    return s | 1u;
}

static inline int rmd_rand_byte(unsigned int *s) {
    *s = *s * 1103515245u + 12345u;
    return (int)((*s >> 16) & 0xFF);
}

#endif /* RMD_RANDOM_H */
