#ifndef PARAMS_H
#define PARAMS_H 

#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>

static inline
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtscp" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

#define TIC printf("\n"); uint64_t cl = rdtsc();
#define TOC(A) printf("%s cycles = %lu \n",#A ,rdtsc() - cl); cl = rdtsc();

#define PRIME 31
#define PRIME_BITS 5
#define M 60  
#define N 62   
#define O 6   
#define K 10  

#define KC2 (K*(K+1)/2)
#define MONOMIALS (N*(N+1)/2)
#define P1MONOMIALS ((N-O)*(N-O+1)/2 + (N-O)*O)
#define P2MONOMIALS (O*(O+1)/2)
#define SEED_BYTES 16
#define HASH_BYTES 32

#define P1_BYTES (M*P1MONOMIALS)
#define P2_BYTES (M*P2MONOMIALS)
#define OIL_SPACE_BYTES (O*(N-O)) 

#define PK_SEED(pk) (pk)
#define PK_P2(pk) (pk + SEED_BYTES)
#define PK_BYTES (PK_P2(0) + M*O*(O+1)/2)

#define PK_EXP_BYTES P1_BYTES

#define SK_EXP_P1(sk_exp) (sk_exp)
#define SK_EXP_OIL(sk_exp) (sk_exp + P1_BYTES)
#define SK_EXP_BILINEAR(sk_exp) (SK_EXP_OIL(sk_exp) + OIL_SPACE_BYTES)
#define SK_EXP_BYTES (SK_EXP_BILINEAR(0) + (M*(N-O)*O))

#define SK_PUBLIC_SEED(sk) (sk)
#define SK_PRIVATE_SEED(sk) (sk + SEED_BYTES)
#define SK_BYTES (2*SEED_BYTES)

#define SIG_SALT(sig) (sig)
#define SIG_INPUTS(sig) (sig + SEED_BYTES)
#define SIG_BYTES (SEED_BYTES + (K*N)) 

#include "libXKCP.a.headers/SimpleFIPS202.h"
#define HASH(data,len,out) SHAKE128(out, HASH_BYTES, data, len);
#define EXPAND(data,len,out,outlen) SHAKE128(out, outlen, data, len);

#endif
