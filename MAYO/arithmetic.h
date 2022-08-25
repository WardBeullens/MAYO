#ifndef ARITHMETIC_H
#define ARITHMETIC_H 
#include "MAYO_params.h"
#include <stdint.h>
#include <stdio.h>
#include <openssl/rand.h>

void negate(unsigned char* v, int len);
uint32_t mod_inverse(uint32_t a);


void add_vectors(const unsigned char *v1, const unsigned char *v2, unsigned char *out);
void _linear_combination_avx(const unsigned char* vecs, const unsigned char* coeffs, int len, unsigned char* out);
void _linear_combination_avx_new(const unsigned char* vecs, const unsigned char* coeffs, int len, unsigned char* out);
void _linear_combination(const unsigned char* vecs, const unsigned char* coeffs, int len, unsigned char* out);
void evaluateP(const unsigned char *inputs, const unsigned char *P1, const unsigned char *P2, unsigned char *out);
void evaluateP_vinegar(const unsigned char *inputs, const unsigned char *P1, unsigned char *out);
void reduce(unsigned char* v, int len);

#if M <= 64
	#define linear_combination _linear_combination_avx
#else
	#define linear_combination _linear_combination
#endif


int sample_oil(const unsigned char *rhs,const unsigned char *linear, unsigned char *solution);

void zero_out(unsigned char *out);

#define print_vec(vec) \
	printf("%s :\n", #vec); \
	for (int print_vec_counter = 0; print_vec_counter < M; ++print_vec_counter) \
	{ \
		printf("%3d ", vec[print_vec_counter]); \
		if(print_vec_counter == 25) { \
			printf("\n"); \
		}\
	} \
	printf("\n");


#endif