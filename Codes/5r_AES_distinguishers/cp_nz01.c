/* cp_nz01.c -- 5-round AES chosen-plaintext (CP) distinguisher: n_z=0 and n_z=1
 * pairings measured on the SAME structure at one data complexity N (Section 4.2).
 * Slice-wise radix sort (O(N)) for the collision search.
 * Structure (3D) {0,10}=s0 x {5}=s1 x {15}=s2; ct[a*slice + j*s2 + k], slice=s1*s2.
 * n_z=0 (base 2-active {5,15}, same {0,10}): base (a,jP,kP),(a,jQ,kQ);
 *        partner (a',jP,kQ),(a',jQ,kP) for a'!=a (byte15 cross-swap, other {0,10}).
 * n_z=1 (base 2-active {0,10}, same {5,15}, Yan-Tan-Qi pairing):
 *        base (a1,w),(a2,w); partner (a1,w'),(a2,w') for w' a 2-active {5,15}-shift of w
 *        (both byte5 and byte15 differ; matches count_quartets_nz1.c offsets).
 * Valid: base and partner both zero-diff on the SAME inverse diagonal; /2.
 * Usage: ./cp_nz01 s0 s1 s2 rounds [trials] [seed]
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

typedef struct { uint32_t v; int idx; } SV;
/* LSD radix sort of SV[n] by the 32-bit key v. O(n). */
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

int main(int argc,char**argv){
    if(!aes_self_test()){ fprintf(stderr,"AES KAT FAILED\n"); return 2; }
    if(argc<5){ fprintf(stderr,"Usage: %s s0 s1 s2 rounds [trials] [seed]\n",argv[0]); return 1; }
    int s0=atoi(argv[1]),s1=atoi(argv[2]),s2=atoi(argv[3]),nr=atoi(argv[4]);
    int trials=(argc>=6)?atoi(argv[5]):200;
    if(nr<1||nr>10){ fprintf(stderr,"rounds must be 1..10\n"); return 1; }
    uint64_t seed0=(argc>=7)?strtoull(argv[6],0,10):42ull;
    int slice=s1*s2; long total=(long)s0*slice;

    uint8_t (*ct)[16]=malloc((size_t)total*16);
    SV *bufA=malloc(sizeof(SV)*slice), *bufAt=malloc(sizeof(SV)*slice);
    SV *bufB=malloc(sizeof(SV)*s0),    *bufBt=malloc(sizeof(SV)*s0);
    uint8_t *v0=malloc(s0),*v10=malloc(s0);
    if(!ct||!bufA||!bufAt||!bufB||!bufBt||!v0||!v10){ fprintf(stderr,"malloc fail (N=%ld ~%.1f GB)\n",total,total*16.0/1e9); return 1; }

    int hit0=0,hit1=0; double sum0=0,sum1=0;
    uint64_t B0t=0,C0t=0,RV0t=0, B1t=0,C1t=0,RV1t=0; /* diag: base collisions, partner checks, raw valids */
    for(int t=0;t<trials;t++){
        rng = seed0 ^ (0x9E3779B97F4A7C15ull*(uint64_t)(t+1));
        uint8_t key[16],bp[16]; for(int x=0;x<16;x++){key[x]=xr();bp[x]=xr();}
        uint32_t rk[44]; aes128_key_expansion(key,rk);
        { static uint16_t *pool=NULL; if(!pool) pool=malloc(sizeof(uint16_t)*(1<<16));
          for(int x=0;x<(1<<16);x++) pool[x]=(uint16_t)x;
          for(int x=0;x<s0;x++){ int r=x+(int)(xr()%((1u<<16)-x)); uint16_t tt=pool[x];pool[x]=pool[r];pool[r]=tt; v0[x]=pool[x]>>8; v10[x]=pool[x]&0xFF; } }
        uint8_t n5[256],n15[256];
        { uint8_t p[256]; for(int x=0;x<256;x++)p[x]=x; for(int x=255;x>0;x--){int r=xr()%(x+1);uint8_t tt=p[x];p[x]=p[r];p[r]=tt;} for(int x=0;x<s1;x++)n5[x]=p[x]; }
        { uint8_t p[256]; for(int x=0;x<256;x++)p[x]=x; for(int x=255;x>0;x--){int r=xr()%(x+1);uint8_t tt=p[x];p[x]=p[r];p[r]=tt;} for(int x=0;x<s2;x++)n15[x]=p[x]; }

        for(int a=0;a<s0;a++)for(int j=0;j<s1;j++)for(int k=0;k<s2;k++){
            uint8_t pt[16]; memcpy(pt,bp,16);
            pt[0]=v0[a]; pt[10]=v10[a]; pt[5]=n5[j]; pt[15]=n15[k];
            aes_enc_rounds(ct[(long)a*slice+j*s2+k],pt,rk,nr);
        }

        /* ---- n_z=0: base 2-active {5,15} within an {0,10}-slice ---- */
        uint64_t c0=0;
        for(int d=0;d<4;d++) for(int a=0;a<s0;a++){
            uint8_t (*sl)[16]=&ct[(long)a*slice];
            for(int idx=0;idx<slice;idx++){ bufA[idx].v=invd(sl[idx],d); bufA[idx].idx=idx; }
            radix_sv(bufA,bufAt,slice);
            for(int p=0;p<slice;){ int q=p+1; while(q<slice&&bufA[q].v==bufA[p].v)q++;
                for(int x=p;x<q;x++)for(int y=x+1;y<q;y++){
                    int iP=bufA[x].idx,iQ=bufA[y].idx; int jP=iP/s2,kP=iP%s2,jQ=iQ/s2,kQ=iQ%s2;
                    if(jP==jQ||kP==kQ) continue;
                    B0t++; C0t+=(uint64_t)(s0-1);
                    int pa=jP*s2+kQ, pb=jQ*s2+kP;
                    for(int ap=0;ap<s0;ap++){ if(ap==a)continue;
                        uint8_t (*sl2)[16]=&ct[(long)ap*slice];
                        if(invd(sl2[pa],d)==invd(sl2[pb],d)) c0++; }
                }
                p=q; }
        }
        RV0t+=c0;
        c0/=2;

        /* ---- n_z=1 on {0,10}: base 2-active {0,10} within a {5,15}-slice (Yan et al.) ---- */
        uint64_t c1=0;
        for(int d=0;d<4;d++) for(int w=0;w<slice;w++){
            int jw=w/s2, kw=w%s2;
            for(int a=0;a<s0;a++){ bufB[a].v=invd(ct[(long)a*slice+w],d); bufB[a].idx=a; }
            radix_sv(bufB,bufBt,s0);
            for(int p=0;p<s0;){ int q=p+1; while(q<s0&&bufB[q].v==bufB[p].v)q++;
                for(int x=p;x<q;x++)for(int y=x+1;y<q;y++){
                    int a1=bufB[x].idx,a2=bufB[y].idx;
                    if(v0[a1]==v0[a2]||v10[a1]==v10[a2]) continue;
                    B1t++; C1t+=(uint64_t)(s1-1)*(uint64_t)(s2-1);
                    for(int wp=0;wp<slice;wp++){ int jwp=wp/s2,kwp=wp%s2; if(jwp==jw||kwp==kw)continue; /* 2-active {5,15} shift only */
                        if(invd(ct[(long)a1*slice+wp],d)==invd(ct[(long)a2*slice+wp],d)) c1++; }
                }
                p=q; }
        }
        RV1t+=c1;
        c1/=2;

        if(c0>0)hit0++;
        sum0+=c0;
        if(c1>0)hit1++;
        sum1+=c1;
        fprintf(stderr,"\r[%d/%d] n0=%llu n1=%llu",t+1,trials,(unsigned long long)c0,(unsigned long long)c1); fflush(stderr);
    }
    fprintf(stderr,"\n");
    printf("AES cp_nz01: s0=%d s1=%d s2=%d N=%ld(2^%.2f) rounds=%d trials=%d\n",s0,s1,s2,total,
           total>0?log2((double)total):0.0,nr,trials);
    printf("  n_z=0: P_find=%d/%d=%.4f mean=%.4f\n",hit0,trials,(double)hit0/trials,sum0/trials);
    printf("  n_z=1: P_find=%d/%d=%.4f mean=%.4f\n",hit1,trials,(double)hit1/trials,sum1/trials);
    printf("  [diag] base-coll/trial : n0=%.1f  n1=%.1f\n",(double)B0t/trials,(double)B1t/trials);
    printf("  [diag] checks/trial    : n0=%.3e n1=%.3e\n",(double)C0t/trials,(double)C1t/trials);
    printf("  [diag] partner rate RV/checks : n0=%.3e (2^%.2f)  n1=%.3e (2^%.2f)\n",
           C0t?(double)RV0t/(double)C0t:0.0, (C0t&&RV0t)?log2((double)RV0t/(double)C0t):0.0,
           C1t?(double)RV1t/(double)C1t:0.0, (C1t&&RV1t)?log2((double)RV1t/(double)C1t):0.0);
    free(ct);free(bufA);free(bufAt);free(bufB);free(bufBt);free(v0);free(v10); return 0;
}
