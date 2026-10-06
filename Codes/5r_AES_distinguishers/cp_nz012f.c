/* cp_nz012f.c -- Section 4.3, Table 13: speed-optimised version of cp_nz012.c with identical T lines
 * and footer (only the header line names cp_nz012f).
 * Same RNG stream, structure, pairings (a),(b),(c) and per-trial output line.
 * Differences (implementation only):
 *   - encryption with AES-NI (common/rmd_aesni.h), checked against aes.c at start-up;
 *   - one inverse diagonal at a time: the structure is re-encrypted for each d
 *     (4 encryptions per text) and only the word D[a*slice+w] of diagonal d is kept;
 *   - for pairing (a), D[d] is transposed per diagonal so the s0-column gather is
 *     contiguous (extra 4N bytes).
 * Usage: ./cp_nz012f s0 s1 s2 rounds trials seed [t_start]
 * Per-trial stdout line (same as cp_nz012):
 *   T t ca cb cc  rawA rawB rawC  baseA baseB  P2both  candC
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "aes.h"
#include "../common/rmd_aesni.h"

static const int INVD[4][4]={{0,13,10,7},{4,1,14,11},{8,5,2,15},{12,9,6,3}};
static inline uint32_t invd(const uint8_t*c,int d){
    return ((uint32_t)c[INVD[d][0]]<<24)|((uint32_t)c[INVD[d][1]]<<16)|((uint32_t)c[INVD[d][2]]<<8)|c[INVD[d][3]];
}
static uint64_t rng;
static uint64_t xr(void){ rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }

typedef struct { uint32_t v; int idx; } SV;
static void radix_sv(SV *a, SV *tmp, int n){
    SV *s=a,*d=tmp;
    for(int pass=0;pass<4;pass++){
        int sh=pass*8; int c[257]; memset(c,0,sizeof(c));
        for(int i=0;i<n;i++) c[((s[i].v>>sh)&0xFF)+1]++;
        for(int b=1;b<=256;b++) c[b]+=c[b-1];
        for(int i=0;i<n;i++) d[c[(s[i].v>>sh)&0xFF]++]=s[i];
        SV*t=s;s=d;d=t;
    }
    if(s!=a) memcpy(a,s,(size_t)n*sizeof(SV));
}

static int aesni_check(int nr){
    uint64_t save=rng; rng=0x123456789ABCDEFull;
    for(int it=0;it<2000;it++){
        uint8_t key[16],pt[16],c1[16],c2[16]; for(int x=0;x<16;x++){key[x]=xr();pt[x]=xr();}
        uint32_t rk[44]; aes128_key_expansion(key,rk); aes_enc_rounds(c1,pt,rk,nr);
        __m128i rkn[11]; rmd_aes128_generate_round_keys(key,rkn); rmd_aes128_encrypt(pt,c2,rkn,nr);
        if(memcmp(c1,c2,16)) { rng=save; return 0; }
    }
    rng=save; return 1;
}

int main(int argc,char**argv){
    if(!aes_self_test()){ fprintf(stderr,"AES KAT FAILED\n"); return 2; }
    if(argc<7){ fprintf(stderr,"Usage: %s s0 s1 s2 rounds trials seed [t_start]\n",argv[0]); return 1; }
    int s0=atoi(argv[1]),s1=atoi(argv[2]),s2=atoi(argv[3]),nr=atoi(argv[4]);
    int trials=atoi(argv[5]);
    uint64_t seed0=strtoull(argv[6],0,10);
    int tstart=(argc>=8)?atoi(argv[7]):0;
    if(nr<1||nr>10){ fprintf(stderr,"rounds must be 1..10\n"); return 1; }
    if(!aesni_check(nr)){ fprintf(stderr,"AES-NI != aes.c\n"); return 2; }
    int slice=s1*s2; long total=(long)s0*slice;

    uint32_t *D=malloc((size_t)total*4);            /* diagonal d only: D[a*slice + w] */
    uint32_t *TT=malloc((size_t)total*4);           /* transposed D[d]: TT[w*s0 + a] */
    SV *bufA=malloc(sizeof(SV)*slice), *bufAt=malloc(sizeof(SV)*slice);
    SV *bufB=malloc(sizeof(SV)*s0),    *bufBt=malloc(sizeof(SV)*s0);
    uint8_t *v0=malloc(s0),*v10=malloc(s0); int *nbothc=malloc(sizeof(int)*s0);
    if(!D||!TT||!bufA||!bufAt||!bufB||!bufBt||!v0||!v10||!nbothc){ fprintf(stderr,"malloc fail\n"); return 1; }

    printf("# cp_nz012f s0=%d s1=%d s2=%d N=%ld(2^%.4f) rounds=%d trials=%d seed=%llu t_start=%d\n",
           s0,s1,s2,total,log2((double)total),nr,trials,(unsigned long long)seed0,tstart);
    double W=(double)s1*(s1-1)*(double)s2*(s2-1)/2.0;
    printf("# Q_b=C(s0,2)*W=%.6e (2^%.4f)  W=%.0f ; Q_a=Q_c=P2both*W (per trial)\n",(double)s0*(s0-1)/2.0*W,log2((double)s0*(s0-1)/2.0*W),W);
    fflush(stdout);
    int hA=0,hB=0,hC=0,hAB=0,hAC=0,hABC=0,hBC=0;
    for(int tt=0;tt<trials;tt++){
        int t=tstart+tt;
        rng = seed0 ^ (0x9E3779B97F4A7C15ull*(uint64_t)(t+1));
        uint8_t key[16],bp[16]; for(int x=0;x<16;x++){key[x]=xr();bp[x]=xr();}
        __m128i rkn[11]; rmd_aes128_generate_round_keys(key,rkn);
        { static uint16_t *pool=NULL; if(!pool) pool=malloc(sizeof(uint16_t)*(1<<16));
          for(int x=0;x<(1<<16);x++) pool[x]=(uint16_t)x;
          for(int x=0;x<s0;x++){ int r=x+(int)(xr()%((1u<<16)-x)); uint16_t tq=pool[x];pool[x]=pool[r];pool[r]=tq; v0[x]=pool[x]>>8; v10[x]=pool[x]&0xFF; } }
        uint8_t n5[256],n15[256];
        { uint8_t p[256]; for(int x=0;x<256;x++)p[x]=x; for(int x=255;x>0;x--){int r=xr()%(x+1);uint8_t tq=p[x];p[x]=p[r];p[r]=tq;} for(int x=0;x<s1;x++)n5[x]=p[x]; }
        { uint8_t p[256]; for(int x=0;x<256;x++)p[x]=x; for(int x=255;x>0;x--){int r=xr()%(x+1);uint8_t tq=p[x];p[x]=p[r];p[r]=tq;} for(int x=0;x<s2;x++)n15[x]=p[x]; }

        uint64_t P2both;
        { uint64_t h0[256]={0},h10[256]={0}; for(int a=0;a<s0;a++){h0[v0[a]]++;h10[v10[a]]++;}
          uint64_t e0=0,e10=0; for(int x=0;x<256;x++){ if(h0[x]) e0+=h0[x]*(h0[x]-1)/2; if(h10[x]) e10+=h10[x]*(h10[x]-1)/2; }
          P2both=(uint64_t)s0*(s0-1)/2-e0-e10;
          for(int a=0;a<s0;a++) nbothc[a]=s0-(int)h0[v0[a]]-(int)h10[v10[a]]+1; }

        uint64_t rb=0, rc=0, baseB=0, candC=0, ra=0, baseA=0;
        for(int d=0;d<4;d++){
            for(int a=0;a<s0;a++)for(int j=0;j<s1;j++)for(int k=0;k<s2;k++){
                uint8_t pt[16],c[16]; memcpy(pt,bp,16);
                pt[0]=v0[a]; pt[10]=v10[a]; pt[5]=n5[j]; pt[15]=n15[k];
                rmd_aes128_encrypt(pt,c,rkn,nr);
                D[(long)a*slice+j*s2+k]=invd(c,d);
            }
            const uint32_t *Dd=D;
            /* ---- (b) and (c) ---- */
            for(int a=0;a<s0;a++){
                const uint32_t *sl=Dd+(long)a*slice;
                for(int idx=0;idx<slice;idx++){ bufA[idx].v=sl[idx]; bufA[idx].idx=idx; }
                radix_sv(bufA,bufAt,slice);
                for(int p=0;p<slice;){ int q=p+1; while(q<slice&&bufA[q].v==bufA[p].v)q++;
                    for(int x=p;x<q;x++)for(int y=x+1;y<q;y++){
                        int iP=bufA[x].idx,iQ=bufA[y].idx; int jP=iP/s2,kP=iP%s2,jQ=iQ/s2,kQ=iQ%s2;
                        if(jP==jQ||kP==kQ) continue;
                        baseB++; candC+=(uint64_t)nbothc[a];
                        int pa=jP*s2+kQ, pb=jQ*s2+kP;
                        for(int ap=0;ap<s0;ap++){ if(ap==a)continue;
                            const uint32_t *sl2=Dd+(long)ap*slice;
                            if(sl2[pa]==sl2[pb]) rb++;
                            if(v0[ap]!=v0[a] && v10[ap]!=v10[a] && sl2[iP]==sl2[iQ]) rc++;
                        }
                    }
                    p=q; }
            }
            /* ---- (a): transpose Dd (s0 x slice) -> TT (slice x s0), blocked ---- */
            const int BLK=64;
            for(int a0=0;a0<s0;a0+=BLK) for(int w0=0;w0<slice;w0+=BLK){
                int a1=a0+BLK<s0?a0+BLK:s0, w1=w0+BLK<slice?w0+BLK:slice;
                for(int a=a0;a<a1;a++){ const uint32_t *src=Dd+(long)a*slice; for(int w=w0;w<w1;w++) TT[(long)w*s0+a]=src[w]; }
            }
            for(int w=0;w<slice;w++){
                int jw=w/s2, kw=w%s2;
                const uint32_t *col=TT+(long)w*s0;
                for(int a=0;a<s0;a++){ bufB[a].v=col[a]; bufB[a].idx=a; }
                radix_sv(bufB,bufBt,s0);
                for(int p=0;p<s0;){ int q=p+1; while(q<s0&&bufB[q].v==bufB[p].v)q++;
                    for(int x=p;x<q;x++)for(int y=x+1;y<q;y++){
                        int a1=bufB[x].idx,a2=bufB[y].idx;
                        if(v0[a1]==v0[a2]||v10[a1]==v10[a2]) continue;
                        baseA++;
                        const uint32_t *r1=Dd+(long)a1*slice, *r2=Dd+(long)a2*slice;
                        for(int wp=0;wp<slice;wp++){ int jwp=wp/s2,kwp=wp%s2; if(jwp==jw||kwp==kw)continue;
                            if(r1[wp]==r2[wp]) ra++; }
                    }
                    p=q; }
            }
        }
        uint64_t ca=ra/2, cb=rb/2, cc=rc/2;
        int fa=ca>0, fb=cb>0, fc=cc>0;
        hA+=fa; hB+=fb; hC+=fc; hAB+=(fa|fb); hAC+=(fa|fc); hBC+=(fb|fc); hABC+=(fa|fb|fc);
        printf("T %d %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu\n",t,
               (unsigned long long)ca,(unsigned long long)cb,(unsigned long long)cc,
               (unsigned long long)ra,(unsigned long long)rb,(unsigned long long)rc,
               (unsigned long long)baseA,(unsigned long long)baseB,(unsigned long long)P2both,
               (unsigned long long)candC);
        fflush(stdout);
    }
    printf("# P_find a=%d b=%d c=%d a|b=%d a|c=%d b|c=%d a|b|c=%d of %d\n",hA,hB,hC,hAB,hAC,hBC,hABC,trials);
    free(D);free(TT);free(bufA);free(bufAt);free(bufB);free(bufBt);free(v0);free(v10);free(nbothc); return 0;
}
