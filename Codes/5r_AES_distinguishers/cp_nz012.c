/* cp_nz012.c -- Section 4.3, Table 13: the pairings n_z = 1 on {0, 10} (a, Yan et al.), n_z = 0 (b)
 * and n_z = 1 on {5, 15} (c) on one product structure, and their union.
 * Extension of Codes/5r_AES_distinguishers/cp_nz01.c.
 * Same structure, same per-trial key/plaintext generation (identical RNG stream),
 * so pairings (a) and (b) reproduce cp_nz01 trial by trial.
 *
 * Structure (3D) {0,10}=s0 x {5}=s1 x {15}=s2; ct[a*slice + j*s2 + k], slice=s1*s2.
 * (a) n_z=1 on {0,10} (Yan et al.): base (a1,w),(a2,w) [bytes 0 AND 10 differ, 5,15 equal];
 *     partner (a1,w'),(a2,w') with w' differing from w at byte 5 AND byte 15.
 * (b) n_z=0: base (a,jP,kP),(a,jQ,kQ) [bytes 5 AND 15 differ, 0,10 equal];
 *     partner (a',jP,kQ),(a',jQ,kP), a'!=a (byte-15 swap + {0,10} shift).
 * (c) n_z=1 on {5,15}: same base as (b) [base pairs differ at bytes 5 and 15];
 *     partner (a',jP,kP),(a',jQ,kQ) with byte 0 AND byte 10 of a' different
 *     from those of a (common nonzero shift of bytes 0 and 10 applied to both
 *     base texts).
 * Valid: base and partner both zero difference on the SAME inverse diagonal.
 * Each quartet is found twice (base<->partner), raw counts are halved.
 *
 * Usage: ./cp_nz012 s0 s1 s2 rounds trials seed [t_start] [verify]
 *   runs trials t = t_start .. t_start+trials-1, with the cp_nz01 seed rule
 *   rng = seed ^ (0x9E3779B97F4A7C15 * (t+1)).
 *   verify=1: additionally recount (a),(b),(c) by an independent brute force
 *   written from the plaintext-byte definitions (small grids only), and check the
 *   round-1 relation Delta_R1(P,Q)=Delta_R1(P',Q') on every candidate quartet.
 * Per-trial stdout line:
 *   T t ca cb cc  rawA rawB rawC  baseA baseB  P2both  candC
 *   (baseA/baseB = base collisions, i.e. colliding base pairs whose two base-active
 *    bytes both differ, summed over the four inverse diagonals; candC = sum over (b)/(c) base
 *    collisions of the number of admissible (c) partners a')
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "aes.h"

static const int INVD[4][4]={{0,13,10,7},{4,1,14,11},{8,5,2,15},{12,9,6,3}};
#ifndef INVD_MASK            /* test builds only: -DINVD_MASK=0xFF000000u compares 8 of the 32 bits */
#define INVD_MASK 0xFFFFFFFFu  /* production: full 32-bit inverse-diagonal collision */
#endif
static inline uint32_t invd(const uint8_t*c,int d){
    return INVD_MASK & (((uint32_t)c[INVD[d][0]]<<24)|((uint32_t)c[INVD[d][1]]<<16)|((uint32_t)c[INVD[d][2]]<<8)|c[INVD[d][3]]);
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

/* ---------- independent brute force (verify mode, small grids) ---------- */
static long NB; static int BS1,BS2; static uint8_t (*PT)[16]; static uint8_t (*CT)[16]; static uint8_t (*R1)[16];
static int invA[65536], invJ[256], invK[256];
static long find_pt(const uint8_t *p){ /* plaintext -> index via product-grid inverse tables */
    int a=invA[(p[0]<<8)|p[10]], j=invJ[p[5]], k=invK[p[15]];
    if(a<0||j<0||k<0) return -1;
    return (long)a*BS1*BS2+(long)j*BS2+k;
}
static int coll(long i,long j,int d){ return invd(CT[i],d)==invd(CT[j],d); }
static int r1rel(long p,long q,long pp,long qq){ for(int x=0;x<16;x++) if((R1[p][x]^R1[q][x])!=(R1[pp][x]^R1[qq][x])) return 0; return 1; }
/* raw counts out[3] normalised to "each quartet twice" (as the fast code);
 * r1bad[3]: candidate quartets violating Delta_R1(base)=Delta_R1(partner); ncand[3]: #quartets */
static void brute(uint64_t out[3], uint64_t r1bad[3], uint64_t ncand[3]){
    for(int i=0;i<3;i++){out[i]=0;r1bad[i]=0;ncand[i]=0;}
    for(long P=0;P<NB;P++) for(long Q=0;Q<NB;Q++){ if(P==Q) continue;
        const uint8_t *p=PT[P],*q=PT[Q];
        int eq0=p[0]==q[0], eq10=p[10]==q[10], eq5=p[5]==q[5], eq15=p[15]==q[15];
        /* (a): base differs at 0 AND 10, equal at 5,15; shift delta at 5,15 both nonzero */
        if(!eq0&&!eq10&&eq5&&eq15){
            for(long Pp=0;Pp<NB;Pp++){ const uint8_t *pp=PT[Pp];
                if(pp[0]!=p[0]||pp[10]!=p[10]||pp[5]==p[5]||pp[15]==p[15]) continue;
                uint8_t qq[16]; memcpy(qq,q,16); qq[5]^=pp[5]^p[5]; qq[15]^=pp[15]^p[15];
                long Qp=find_pt(qq); if(Qp<0) continue;
                ncand[0]++; if(!r1rel(P,Q,Pp,Qp)) r1bad[0]++;
                for(int d=0;d<4;d++) if(coll(P,Q,d)&&coll(Pp,Qp,d)) out[0]++;
            }
        }
        if(eq0&&eq10&&!eq5&&!eq15){
            for(long Pp=0;Pp<NB;Pp++){ const uint8_t *pp=PT[Pp];
                /* (b): partner P^delta, Q^delta, delta=(d0,d10,0,p15^q15), (d0,d10)!=(0,0) */
                if(pp[5]==p[5]&&pp[15]==q[15]&&!(pp[0]==p[0]&&pp[10]==p[10])){
                    uint8_t qq[16]; memcpy(qq,q,16); qq[0]^=pp[0]^p[0]; qq[10]^=pp[10]^p[10]; qq[15]^=p[15]^q[15];
                    long Qp=find_pt(qq); if(Qp>=0){ ncand[1]++; if(!r1rel(P,Q,Pp,Qp)) r1bad[1]++;
                        for(int d=0;d<4;d++) if(coll(P,Q,d)&&coll(Pp,Qp,d)) out[1]++; }
                }
                /* (c): partner P^delta, Q^delta, delta at bytes 0 AND 10 both nonzero */
                if(pp[5]==p[5]&&pp[15]==p[15]&&pp[0]!=p[0]&&pp[10]!=p[10]){
                    uint8_t qq[16]; memcpy(qq,q,16); qq[0]^=pp[0]^p[0]; qq[10]^=pp[10]^p[10];
                    long Qp=find_pt(qq); if(Qp>=0){ ncand[2]++; if(!r1rel(P,Q,Pp,Qp)) r1bad[2]++;
                        for(int d=0;d<4;d++) if(coll(P,Q,d)&&coll(Pp,Qp,d)) out[2]++; }
                }
            }
        }
    }
    /* ordered (P,Q) visits each unordered base pair twice, and every quartet is reached
       from both of its pairs -> raw = 4 x (#quartet,d); fast code reports 2 x */
    for(int i=0;i<3;i++){ out[i]/=2; ncand[i]/=4; r1bad[i]/=4; }
}

int main(int argc,char**argv){
    if(!aes_self_test()){ fprintf(stderr,"AES KAT FAILED\n"); return 2; }
    if(argc<7){ fprintf(stderr,"Usage: %s s0 s1 s2 rounds trials seed [t_start] [verify]\n",argv[0]); return 1; }
    int s0=atoi(argv[1]),s1=atoi(argv[2]),s2=atoi(argv[3]),nr=atoi(argv[4]);
    int trials=atoi(argv[5]);
    uint64_t seed0=strtoull(argv[6],0,10);
    int tstart=(argc>=8)?atoi(argv[7]):0;
    int verify=(argc>=9)?atoi(argv[8]):0;
    if(nr<1||nr>10){ fprintf(stderr,"rounds must be 1..10\n"); return 1; }
    int slice=s1*s2; long total=(long)s0*slice;

    uint8_t (*ct)[16]=malloc((size_t)total*16);
    SV *bufA=malloc(sizeof(SV)*slice), *bufAt=malloc(sizeof(SV)*slice);
    SV *bufB=malloc(sizeof(SV)*s0),    *bufBt=malloc(sizeof(SV)*s0);
    uint8_t *v0=malloc(s0),*v10=malloc(s0); int *nbothc=malloc(sizeof(int)*s0);
    if(!ct||!bufA||!bufAt||!bufB||!bufBt||!v0||!v10||!nbothc){ fprintf(stderr,"malloc fail\n"); return 1; }
    if(verify){ NB=total; BS1=s1; BS2=s2; PT=malloc((size_t)total*16); CT=ct; R1=malloc((size_t)total*16); }

    printf("# cp_nz012 s0=%d s1=%d s2=%d N=%ld(2^%.4f) rounds=%d trials=%d seed=%llu t_start=%d\n",
           s0,s1,s2,total,log2((double)total),nr,trials,(unsigned long long)seed0,tstart);
    double W=(double)s1*(s1-1)*(double)s2*(s2-1)/2.0; /* unordered {w,w'} differing in 5 AND 15 */
    printf("# Q_b=C(s0,2)*W=%.6e (2^%.4f)  W=%.0f ; Q_a=Q_c=P2both*W (per trial)\n",(double)s0*(s0-1)/2.0*W,log2((double)s0*(s0-1)/2.0*W),W);
    fflush(stdout);
    int hA=0,hB=0,hC=0,hAB=0,hAC=0,hABC=0,hBC=0;
    int vfail=0;
    for(int tt=0;tt<trials;tt++){
        int t=tstart+tt;
        rng = seed0 ^ (0x9E3779B97F4A7C15ull*(uint64_t)(t+1));
        uint8_t key[16],bp[16]; for(int x=0;x<16;x++){key[x]=xr();bp[x]=xr();}
        uint32_t rk[44]; aes128_key_expansion(key,rk);
        { static uint16_t *pool=NULL; if(!pool) pool=malloc(sizeof(uint16_t)*(1<<16));
          for(int x=0;x<(1<<16);x++) pool[x]=(uint16_t)x;
          for(int x=0;x<s0;x++){ int r=x+(int)(xr()%((1u<<16)-x)); uint16_t tq=pool[x];pool[x]=pool[r];pool[r]=tq; v0[x]=pool[x]>>8; v10[x]=pool[x]&0xFF; } }
        uint8_t n5[256],n15[256];
        { uint8_t p[256]; for(int x=0;x<256;x++)p[x]=x; for(int x=255;x>0;x--){int r=xr()%(x+1);uint8_t tq=p[x];p[x]=p[r];p[r]=tq;} for(int x=0;x<s1;x++)n5[x]=p[x]; }
        { uint8_t p[256]; for(int x=0;x<256;x++)p[x]=x; for(int x=255;x>0;x--){int r=xr()%(x+1);uint8_t tq=p[x];p[x]=p[r];p[r]=tq;} for(int x=0;x<s2;x++)n15[x]=p[x]; }

        if(verify){
            for(int x=0;x<65536;x++) invA[x]=-1;
            for(int x=0;x<256;x++){invJ[x]=-1;invK[x]=-1;}
            for(int a=0;a<s0;a++) invA[(v0[a]<<8)|v10[a]]=a;
            for(int j=0;j<s1;j++) invJ[n5[j]]=j;
            for(int k=0;k<s2;k++) invK[n15[k]]=k;
        }
        for(int a=0;a<s0;a++)for(int j=0;j<s1;j++)for(int k=0;k<s2;k++){
            uint8_t pt[16]; memcpy(pt,bp,16);
            pt[0]=v0[a]; pt[10]=v10[a]; pt[5]=n5[j]; pt[15]=n15[k];
            long id=(long)a*slice+j*s2+k;
            aes_enc_rounds(ct[id],pt,rk,nr);
            if(verify){ memcpy(PT[id],pt,16); aes_enc_rounds(R1[id],pt,rk,1); /* SR(SB(p^k0))^k1 */ }
        }
        /* P2both: unordered pairs a1<a2 with byte0 AND byte10 different;
           nbothc[a] = #a' with byte0 AND byte10 different from a */
        uint64_t P2both;
        { uint64_t h0[256]={0},h10[256]={0}; for(int a=0;a<s0;a++){h0[v0[a]]++;h10[v10[a]]++;}
          uint64_t e0=0,e10=0; for(int x=0;x<256;x++){ if(h0[x]) e0+=h0[x]*(h0[x]-1)/2; if(h10[x]) e10+=h10[x]*(h10[x]-1)/2; }
          P2both=(uint64_t)s0*(s0-1)/2-e0-e10; /* (v0,v10) distinct, so no pair equal in both */
          for(int a=0;a<s0;a++) nbothc[a]=s0-(int)h0[v0[a]]-(int)h10[v10[a]]+1; }

        /* ---- (b) n_z=0 and (c) n_z=1 on {5,15}: base 2-active {5,15} within an {0,10}-slice ---- */
        uint64_t rb=0, rc=0, baseB=0, candC=0;
        for(int d=0;d<4;d++) for(int a=0;a<s0;a++){
            uint8_t (*sl)[16]=&ct[(long)a*slice];
            for(int idx=0;idx<slice;idx++){ bufA[idx].v=invd(sl[idx],d); bufA[idx].idx=idx; }
            radix_sv(bufA,bufAt,slice);
            for(int p=0;p<slice;){ int q=p+1; while(q<slice&&bufA[q].v==bufA[p].v)q++;
                for(int x=p;x<q;x++)for(int y=x+1;y<q;y++){
                    int iP=bufA[x].idx,iQ=bufA[y].idx; int jP=iP/s2,kP=iP%s2,jQ=iQ/s2,kQ=iQ%s2;
                    if(jP==jQ||kP==kQ) continue;
                    baseB++; candC+=(uint64_t)nbothc[a];
                    int pa=jP*s2+kQ, pb=jQ*s2+kP;
                    for(int ap=0;ap<s0;ap++){ if(ap==a)continue;
                        uint8_t (*sl2)[16]=&ct[(long)ap*slice];
                        if(invd(sl2[pa],d)==invd(sl2[pb],d)) rb++;                       /* (b) */
                        if(v0[ap]!=v0[a] && v10[ap]!=v10[a] &&
                           invd(sl2[iP],d)==invd(sl2[iQ],d)) rc++;                        /* (c) */
                    }
                }
                p=q; }
        }

        /* ---- (a) n_z=1 on {0,10}: base 2-active {0,10} within a {5,15}-slice (Yan et al.) ---- */
        uint64_t ra=0, baseA=0;
        for(int d=0;d<4;d++) for(int w=0;w<slice;w++){
            int jw=w/s2, kw=w%s2;
            for(int a=0;a<s0;a++){ bufB[a].v=invd(ct[(long)a*slice+w],d); bufB[a].idx=a; }
            radix_sv(bufB,bufBt,s0);
            for(int p=0;p<s0;){ int q=p+1; while(q<s0&&bufB[q].v==bufB[p].v)q++;
                for(int x=p;x<q;x++)for(int y=x+1;y<q;y++){
                    int a1=bufB[x].idx,a2=bufB[y].idx;
                    if(v0[a1]==v0[a2]||v10[a1]==v10[a2]) continue;
                    baseA++;
                    for(int wp=0;wp<slice;wp++){ int jwp=wp/s2,kwp=wp%s2; if(jwp==jw||kwp==kw)continue;
                        if(invd(ct[(long)a1*slice+wp],d)==invd(ct[(long)a2*slice+wp],d)) ra++; }
                }
                p=q; }
        }
        uint64_t ca=ra/2, cb=rb/2, cc=rc/2;
        int fa=ca>0, fb=cb>0, fc=cc>0;
        hA+=fa; hB+=fb; hC+=fc; hAB+=(fa|fb); hAC+=(fa|fc); hBC+=(fb|fc); hABC+=(fa|fb|fc);
        printf("T %d %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu\n",t,
               (unsigned long long)ca,(unsigned long long)cb,(unsigned long long)cc,
               (unsigned long long)ra,(unsigned long long)rb,(unsigned long long)rc,
               (unsigned long long)baseA,(unsigned long long)baseB,(unsigned long long)P2both,
               (unsigned long long)candC);
        if(verify){
            uint64_t bo[3],r1b[3],nc[3]; brute(bo,r1b,nc);
            double Qa=(double)P2both*W, Qb=(double)s0*(s0-1)/2.0*W;
            int ok=(bo[0]==ra&&bo[1]==rb&&bo[2]==rc&&r1b[0]==0&&r1b[1]==0&&r1b[2]==0&&
                    (double)nc[0]==Qa&&(double)nc[2]==Qa&&(double)nc[1]==Qb);
            if(!ok) vfail++;
            printf("V %d fast_raw(a,b,c)=(%llu,%llu,%llu) brute_raw=(%llu,%llu,%llu) R1viol(a,b,c)=(%llu,%llu,%llu) ncand_brute=(%llu,%llu,%llu) ncand_formula=(%.0f,%.0f,%.0f) id3: baseA*(s1-1)(s2-1)=%llu candC=%llu %s\n",t,
               (unsigned long long)ra,(unsigned long long)rb,(unsigned long long)rc,
               (unsigned long long)bo[0],(unsigned long long)bo[1],(unsigned long long)bo[2],
               (unsigned long long)r1b[0],(unsigned long long)r1b[1],(unsigned long long)r1b[2],
               (unsigned long long)nc[0],(unsigned long long)nc[1],(unsigned long long)nc[2],Qa,Qb,Qa,
               (unsigned long long)(baseA*(uint64_t)(s1-1)*(uint64_t)(s2-1)),(unsigned long long)candC,
               ok?"OK":"MISMATCH");
        }
        fflush(stdout);
    }
    printf("# P_find a=%d b=%d c=%d a|b=%d a|c=%d b|c=%d a|b|c=%d of %d%s\n",hA,hB,hC,hAB,hAC,hBC,hABC,trials,
           verify?(vfail?"  VERIFY_FAIL":"  VERIFY_ALL_OK"):"");
    free(ct);free(bufA);free(bufAt);free(bufB);free(bufBt);free(v0);free(v10);free(nbothc); return 0;
}
