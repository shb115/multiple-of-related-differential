/* cp_nz2.c -- 5-round AES chosen-plaintext (CP) distinguisher, n_z=2 pairing
 * measured on one structure at data complexity N (Section 4.2; radix-sorted
 * global search).
 * 2D structure {0,10}=s0 x {5,15}=s1; ct[a*s1 + w].
 * base 4-active (a1,w1),(a2,w2) [byte0,byte10,byte5,byte15 all differ];
 * partner = swap {0,10}: (a2,w1),(a1,w2). Valid: both zero-diff on same inv diag; /2.
 * Usage: ./cp_nz2 s0 s1 rounds [trials] [seed]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "aes.h"

static const int INVD[4][4]={{0,13,10,7},{4,1,14,11},{8,5,2,15},{12,9,6,3}};
static inline uint32_t invd(const uint8_t*c,int d){
    return ((uint32_t)c[INVD[d][0]]<<24)|((uint32_t)c[INVD[d][1]]<<16)|((uint32_t)c[INVD[d][2]]<<8)|c[INVD[d][3]];
}
static uint64_t rng;
static uint64_t xr(void){ rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }

/* LSD radix sort of key[n] (uint64_t = inv_diag<<32 | idx) by the high 32 bits. O(n). */
static void radix32(uint64_t *key, uint64_t *tmp, long n){
    uint64_t *s=key,*d=tmp;
    for(int pass=0;pass<4;pass++){
        int sh=32+pass*8; long c[257]; memset(c,0,sizeof(c));
        for(long i=0;i<n;i++) c[((s[i]>>sh)&0xFF)+1]++;
        for(int b=1;b<=256;b++) c[b]+=c[b-1];
        for(long i=0;i<n;i++) d[c[(s[i]>>sh)&0xFF]++]=s[i];
        uint64_t*t=s;s=d;d=t;
    }
    if(s!=key) memcpy(key,s,(size_t)n*sizeof(uint64_t));
}

int main(int argc,char**argv){
    if(!aes_self_test()){ fprintf(stderr,"AES KAT FAILED\n"); return 2; }
    if(argc<4){ fprintf(stderr,"Usage: %s s0 s1 rounds [trials] [seed]\n",argv[0]); return 1; }
    int s0=atoi(argv[1]),s1=atoi(argv[2]),nr=atoi(argv[3]);
    int trials=(argc>=5)?atoi(argv[4]):100;
    if(nr<1||nr>10){ fprintf(stderr,"rounds must be 1..10\n"); return 1; }
    uint64_t seed0=(argc>=6)?strtoull(argv[5],0,10):42ull;
    long total=(long)s0*s1;

    uint8_t (*ct)[16]=malloc((size_t)total*16);
    uint64_t *key=malloc((size_t)total*8), *tmp=malloc((size_t)total*8);
    uint8_t *v0=malloc(s0),*v10=malloc(s0),*w5=malloc(s1),*w15=malloc(s1);
    if(!ct||!key||!tmp||!v0||!v10||!w5||!w15){ fprintf(stderr,"malloc fail (N=%ld ~%.1f GB)\n",total,total*32.0/1e9); return 1; }

    int hit=0; double sum=0;
    for(int t=0;t<trials;t++){
        rng = seed0 ^ (0x9E3779B97F4A7C15ull*(uint64_t)(t+1));
        uint8_t mk[16],bp[16]; for(int x=0;x<16;x++){mk[x]=xr();bp[x]=xr();}
        uint32_t rk[44]; aes128_key_expansion(mk,rk);
        { static uint16_t *pool=NULL; if(!pool) pool=malloc(sizeof(uint16_t)*(1<<16));
          for(int x=0;x<(1<<16);x++) pool[x]=(uint16_t)x;
          for(int x=0;x<s0;x++){ int r=x+(int)(xr()%((1u<<16)-x)); uint16_t tt=pool[x];pool[x]=pool[r];pool[r]=tt; v0[x]=pool[x]>>8; v10[x]=pool[x]&0xFF; } }
        { static uint16_t *pool=NULL; if(!pool) pool=malloc(sizeof(uint16_t)*(1<<16));
          for(int x=0;x<(1<<16);x++) pool[x]=(uint16_t)x;
          for(int x=0;x<s1;x++){ int r=x+(int)(xr()%((1u<<16)-x)); uint16_t tt=pool[x];pool[x]=pool[r];pool[r]=tt; w5[x]=pool[x]>>8; w15[x]=pool[x]&0xFF; } }

        for(int a=0;a<s0;a++)for(int w=0;w<s1;w++){
            uint8_t pt[16]; memcpy(pt,bp,16);
            pt[0]=v0[a]; pt[10]=v10[a]; pt[5]=w5[w]; pt[15]=w15[w];
            aes_enc_rounds(ct[(long)a*s1+w],pt,rk,nr);
        }
        uint64_t count=0;
        for(int d=0;d<4;d++){
            for(long idx=0;idx<total;idx++) key[idx]=((uint64_t)invd(ct[idx],d)<<32)|(uint64_t)idx;
            radix32(key,tmp,total);
            for(long p=0;p<total;){ long q=p+1; uint32_t v=(uint32_t)(key[p]>>32);
                while(q<total && (uint32_t)(key[q]>>32)==v) q++;
                for(long x=p;x<q;x++)for(long y=x+1;y<q;y++){
                    long i1=key[x]&0xFFFFFFFF,i2=key[y]&0xFFFFFFFF;
                    int a1=i1/s1,w1=i1%s1,a2=i2/s1,w2=i2%s1;
                    if(v0[a1]==v0[a2]||v10[a1]==v10[a2]||w5[w1]==w5[w2]||w15[w1]==w15[w2]) continue; /* 4-active */
                    if(invd(ct[(long)a2*s1+w1],d)==invd(ct[(long)a1*s1+w2],d)) count++;
                }
                p=q; }
        }
        count/=2;
        if(count>0)hit++;
        sum+=count;
        fprintf(stderr,"\r[%d/%d] n2=%llu",t+1,trials,(unsigned long long)count); fflush(stderr);
    }
    fprintf(stderr,"\n");
    printf("AES n_z=2: s0=%d s1=%d N=%ld(2^%.2f) rounds=%d trials=%d\n",s0,s1,total,
           total>0?log2((double)total):0.0,nr,trials);
    printf("  P_find=%d/%d=%.4f mean=%.4f\n",hit,trials,(double)hit/trials,sum/trials);
    free(ct);free(key);free(tmp);free(v0);free(v10);free(w5);free(w15); return 0;
}
