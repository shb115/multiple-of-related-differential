/* acp_all.c -- 5-round AES adaptive-chosen-plaintext (ACP) distinguisher, one
 * pairing (n_z = 0, 1, or 2) per run, all three reachable via the nz argument
 * (Section 5.2).  Success rate is measured at a FIXED data budget D so the
 * three pairings are compared at the SAME data complexity.  Partners of each
 * base collision are encrypted ON DEMAND (adaptive), as in Yan-Tan-Qi.
 *
 *  nz=1 (YTQ): base structures vary {0,10} (<=2^16/struct); partner shifts {5,15}
 *              by Dx' over all (2^8-1)^2 values: (P0+Dx', P1+Dx').
 *  nz=0:       base structures vary {5,15} (<=2^16/struct, {0,10}=a const); partner
 *              moves {0,10} to a'!=a and cross-swaps byte15:
 *              base (a,jP,kP),(a,jQ,kQ) -> partner (a',jP,kQ),(a',jQ,kP), a' over 2^16-1.
 *  nz=2:       base structure = M sampled 4-active plaintexts (vary {0,5,10,15});
 *              partner = single swap of {0,10}: base (a1,w1),(a2,w2) -> (a2,w1),(a1,w2).
 *
 * Valid quartet: base pair AND a partner pair both zero-diff on the SAME inverse diagonal.
 * Usage: ./acp_all nz Dexp rounds [trials] [seed]      (D = 2^Dexp total queries)
 *        defaults: trials = 200, seed = 42.
 * Committed results (Results/5r_AES_distinguishers/acp_all_2p2332.txt: Table 14 and
 * the 2^22 block): 300 trials, seed 2024, for nz = 1, 0, 2 at Dexp = 22 and 23.32.
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
    if(argc<4){ fprintf(stderr,"Usage: %s nz Dexp rounds [trials] [seed]\n",argv[0]); return 1; }
    char *nz_end; long nzl=strtol(argv[1],&nz_end,10);
    if(argv[1][0]=='\0'||*nz_end||nzl<0||nzl>2){ fprintf(stderr,"nz must be 0, 1 or 2\n"); return 1; }
    int nz=(int)nzl; double Dexp=atof(argv[2]); int nr=atoi(argv[3]);
    int trials=(argc>=5)?atoi(argv[4]):200;
    if(nr<1||nr>10){ fprintf(stderr,"rounds must be 1..10\n"); return 1; }
    uint64_t seed0=(argc>=6)?strtoull(argv[5],0,10):42ull;
    uint64_t D=(uint64_t)(pow(2.0,Dexp)+0.5);
    const int S=1<<16;

    uint8_t *ct=NULL; SV *sv=NULL,*svt=NULL;
    uint8_t (*ct2)[16]=NULL; uint64_t *key=NULL,*tmp=NULL; uint8_t *bb0=NULL,*bb5=NULL,*bb10=NULL,*bb15=NULL;
    long Mmax=0;
    if(nz==0||nz==1){
        ct=malloc((size_t)S*16); sv=malloc(sizeof(SV)*S); svt=malloc(sizeof(SV)*S);
        if(!ct||!sv||!svt){ fprintf(stderr,"malloc fail\n"); return 1; }
    } else {
        Mmax=(long)D+16;
        ct2=malloc((size_t)Mmax*16); key=malloc((size_t)Mmax*8); tmp=malloc((size_t)Mmax*8);
        bb0=malloc(Mmax); bb5=malloc(Mmax); bb10=malloc(Mmax); bb15=malloc(Mmax);
        if(!ct2||!key||!tmp||!bb0||!bb5||!bb10||!bb15){ fprintf(stderr,"malloc fail (M=%ld ~%.1f GB)\n",Mmax,Mmax*40.0/1e9); return 1; }
    }
    uint8_t (*c16)[16]=(uint8_t(*)[16])ct;

    int hit=0; double sum=0,datasum=0; uint64_t hist[64]={0};
    uint64_t Bt=0,Ct=0,RVt=0; /* diag: base collisions, partner pairs queried (N_p), raw valids */
    for(int t=0;t<trials;t++){
        rng = seed0 ^ (0x9E3779B97F4A7C15ull*(uint64_t)(t+1));
        uint8_t mk[16]; for(int x=0;x<16;x++) mk[x]=xr();
        uint32_t rk[44]; aes128_key_expansion(mk,rk);
        uint64_t valid=0, data=0;

        if(nz==1){
            while(data<D){
                uint8_t base[16]; for(int x=0;x<16;x++) base[x]=xr();
                uint64_t rem=D-data; int lim=(rem<(uint64_t)S)?(int)rem:S;
                for(int i=0;i<lim;i++){ uint8_t pt[16]; memcpy(pt,base,16);
                    pt[0]=i>>8; pt[10]=i&0xFF; aes_enc_rounds(c16[i],pt,rk,nr); }
                data+=lim; if(data>=D) break;
                uint8_t c5=base[5],c15=base[15];
                for(int d=0;d<4;d++){
                    for(int i=0;i<S;i++){ sv[i].v=invd(c16[i],d); sv[i].idx=i; }
                    radix_sv(sv,svt,S);
                    for(int p=0;p<S;){ int q=p+1; while(q<S&&sv[q].v==sv[p].v)q++;
                        for(int x=p;x<q;x++)for(int y=x+1;y<q;y++){
                            int i1=sv[x].idx,i2=sv[y].idx;
                            int b0_1=i1>>8,b10_1=i1&0xFF,b0_2=i2>>8,b10_2=i2&0xFF;
                            if(b0_1==b0_2||b10_1==b10_2) continue; Bt++;
                            for(int dd=1;dd<S;dd++){
                                if((dd>>8)==0||(dd&0xFF)==0) continue; /* clean n_z=1: 2-active {5,15} shift only */
                                if(data>=D) goto done1;
                                uint8_t P2[16],P3[16],C2[16],C3[16];
                                memcpy(P2,base,16); memcpy(P3,base,16);
                                P2[0]=b0_1;P2[10]=b10_1; P3[0]=b0_2;P3[10]=b10_2;
                                P2[5]=P3[5]=c5^(dd>>8); P2[15]=P3[15]=c15^(dd&0xFF);
                                aes_enc_rounds(C2,P2,rk,nr); aes_enc_rounds(C3,P3,rk,nr); data+=2; Ct++;
                                if(invd(C2,d)==invd(C3,d)) valid++;
                            }
                        }
                        p=q; }
                }
            }
            done1:;
        } else if(nz==0){
            while(data<D){
                uint8_t base[16]; for(int x=0;x<16;x++) base[x]=xr();
                uint8_t a0=base[0],a10=base[10];
                uint64_t rem=D-data; int lim=(rem<(uint64_t)S)?(int)rem:S;
                for(int i=0;i<lim;i++){ uint8_t pt[16]; memcpy(pt,base,16);
                    pt[5]=i>>8; pt[15]=i&0xFF; aes_enc_rounds(c16[i],pt,rk,nr); }
                data+=lim; if(data>=D) break;
                for(int d=0;d<4;d++){
                    for(int i=0;i<S;i++){ sv[i].v=invd(c16[i],d); sv[i].idx=i; }
                    radix_sv(sv,svt,S);
                    for(int p=0;p<S;){ int q=p+1; while(q<S&&sv[q].v==sv[p].v)q++;
                        for(int x=p;x<q;x++)for(int y=x+1;y<q;y++){
                            int i1=sv[x].idx,i2=sv[y].idx;
                            int jP=i1>>8,kP=i1&0xFF,jQ=i2>>8,kQ=i2&0xFF;
                            if(jP==jQ||kP==kQ) continue; Bt++;
                            for(int ap=0;ap<S;ap++){
                                int ap0=ap>>8,ap10=ap&0xFF; if(ap0==a0&&ap10==a10) continue;
                                if(data>=D) goto done0;
                                uint8_t P2[16],P3[16],C2[16],C3[16];
                                memcpy(P2,base,16); memcpy(P3,base,16);
                                P2[0]=P3[0]=ap0; P2[10]=P3[10]=ap10;
                                P2[5]=jP; P2[15]=kQ;  P3[5]=jQ; P3[15]=kP;
                                aes_enc_rounds(C2,P2,rk,nr); aes_enc_rounds(C3,P3,rk,nr); data+=2; Ct++;
                                if(invd(C2,d)==invd(C3,d)) valid++;
                            }
                        }
                        p=q; }
                }
            }
            done0:;
        } else { /* nz==2 : base-dominated; reserve a little budget for the single swaps */
            uint8_t base[16]; for(int x=0;x<16;x++) base[x]=xr();
            uint64_t reserve=(uint64_t)((double)D*(double)D/1073741824.0);
            long M=(long)((D>reserve)?(D-reserve):(D/2)); if(M>Mmax)M=Mmax;
            for(long i=0;i<M;i++){ bb0[i]=xr();bb5[i]=xr();bb10[i]=xr();bb15[i]=xr();
                uint8_t pt[16]; memcpy(pt,base,16);
                pt[0]=bb0[i];pt[5]=bb5[i];pt[10]=bb10[i];pt[15]=bb15[i];
                aes_enc_rounds(ct2[i],pt,rk,nr); }
            data+=M;
            for(int d=0;d<4;d++){
                for(long i=0;i<M;i++) key[i]=((uint64_t)invd(ct2[i],d)<<32)|(uint64_t)i;
                radix32(key,tmp,M);
                for(long p=0;p<M;){ long q=p+1; uint32_t v=(uint32_t)(key[p]>>32);
                    while(q<M&&(uint32_t)(key[q]>>32)==v)q++;
                    for(long x=p;x<q;x++)for(long y=x+1;y<q;y++){
                        long i1=key[x]&0xFFFFFFFF,i2=key[y]&0xFFFFFFFF;
                        if(bb0[i1]==bb0[i2]||bb10[i1]==bb10[i2]||bb5[i1]==bb5[i2]||bb15[i1]==bb15[i2]) continue;
                        Bt++;
                        if(data>=D) goto done2;
                        uint8_t P2[16],P3[16],C2[16],C3[16];
                        memcpy(P2,base,16); memcpy(P3,base,16);
                        P2[0]=bb0[i2];P2[5]=bb5[i1];P2[10]=bb10[i2];P2[15]=bb15[i1];
                        P3[0]=bb0[i1];P3[5]=bb5[i2];P3[10]=bb10[i1];P3[15]=bb15[i2];
                        aes_enc_rounds(C2,P2,rk,nr); aes_enc_rounds(C3,P3,rk,nr); data+=2; Ct++;
                        if(invd(C2,d)==invd(C3,d)) valid++;
                    }
                    p=q; }
            }
            done2:;
        }

        if(valid>0)hit++; sum+=valid; datasum+=(double)data; RVt+=valid; hist[valid<63?(int)valid:63]++;
        fprintf(stderr,"\r[%d/%d] valid=%llu data=2^%.2f",t+1,trials,(unsigned long long)valid,
                data>0?log((double)data)/log(2.0):0); fflush(stderr);
    }
    fprintf(stderr,"\n");
    printf("ACP nz=%d D=2^%.2f rounds=%d trials=%d : P_find=%d/%d=%.4f mean=%.4f data~2^%.3f\n",
           nz,Dexp,nr,trials,hit,trials,(double)hit/trials,sum/trials,
           datasum>0?log(datasum/trials)/log(2.0):0);
    printf("  [diag] base-coll/trial=%.1f  N_p/trial=%.3e  rate=RV/N_p=%.3e (2^%.2f)\n",
           (double)Bt/trials,(double)Ct/trials,(Ct&&RVt)?(double)RVt/(double)Ct:0.0,
           (Ct&&RVt)?log((double)RVt/(double)Ct)/log(2.0):0.0);
    printf("  [hist] trials by valid count 0..11 (12..62 not shown):");
    for(int h=0;h<12;h++) printf(" %llu",(unsigned long long)hist[h]);
    printf("   [>=63 bucket:%llu]\n",(unsigned long long)hist[63]);
    return 0;
}
