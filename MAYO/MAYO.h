#ifndef MAYO_H
#define MAYO_H 
#include "MAYO_params.h"
#include "arithmetic.h"
#include <openssl/rand.h>
#include <assert.h>

void computeP2(const unsigned char* oil_space, const unsigned char* P1, unsigned char* P2);

void sample_oil_space(const unsigned char *seed, unsigned char *oil_space);
void compute_bilinear_part(const unsigned char *P1, const unsigned char *oil_space, unsigned char *bilinear);
void add_oil(unsigned char *inputs, const unsigned char *oil, const unsigned char *oil_space);

int keygen(unsigned char* pk, unsigned char* sk);

int sign(const unsigned char* m , long long m_len, const unsigned char* sk, unsigned char* sig);
int verify(const unsigned char* m, long long m_len, const unsigned char* pk, unsigned char* sig);

// fast functions
void expand_pk(const unsigned char *pk, unsigned char *pk_exp);
void expand_sk(const unsigned char *sk, unsigned char *sk_exp);
int sign_fast(const unsigned char* m , long long m_len, const unsigned char* sk, const unsigned char* sk_exp, unsigned char* sig);
int verify_fast(const unsigned char* m, long long m_len, const unsigned char* pk, const unsigned char* pk_exp, unsigned char* sig);

#endif
