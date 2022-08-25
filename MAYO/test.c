#include <stdio.h>
#include <stdint.h>
#include "MAYO.h"
#include "MAYO_params.h"
#include <time.h>
#include <stdlib.h>

#define TRIALS 2000

#define MESSAGE_LENGTH 100


#define TEST(fn) \
printf("%30s: ", #fn); \
if(fn()){ \
	printf("Failed \n"); \
	fails ++; \
} else { \
	printf("Success \n"); \
} 

unsigned char test_randomness_seed[1] = {42};

#include "eval_tests.c"
#include "key_gen_tests.c"
#include "arithmetic_tests.c"
#include "signing_tests.c"

int tests(){
	srand(time(NULL));

	int fails = 0;

	printf("Do unit tests:\n");
	TEST(test_unit_vector_eval);
	TEST(test_double_unit_vector_eval);
	TEST(test_double_eval);
	TEST(test_zero_oil_space);
	TEST(test_almost_zero_oil_space);
	TEST(test_eval_zero_oil_vector);
	TEST(test_eval_oil_vector);

	TEST(test_modular_inverse);
	TEST(test_sample_solution);

	TEST(test_bilinear_zero_oilspace);
	TEST(test_bilinear);
	TEST(test_bilinear2);

	TEST(test_linear_combination_avx);

	printf("\n");

	return fails;
}


int main(){

	if( tests()){
		return -1;
	}

	printf("Benchmarking:\n");

	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};

	printf("pk bytes : %d\n", PK_BYTES );
	printf("sk bytes : %d\n", SK_BYTES );
	printf("sig bytes : %ld\n", (uint64_t) SIG_BYTES );

	printf("big pk bytes : %d\n", PK_EXP_BYTES );
	printf("big sk bytes : %d\n", SK_EXP_BYTES );

	unsigned char message[MESSAGE_LENGTH];
	message[0] = 42;
	unsigned char sig[SIG_BYTES];
	unsigned char pk_exp[PK_EXP_BYTES];
	unsigned char sk_exp[SK_EXP_BYTES];

	uint64_t keygenTime = 0;
	uint64_t signTime = 0;
	uint64_t verifyTime = 0;
	uint64_t signFastTime = 0;
	uint64_t verifyFastTime = 0;
	uint64_t t;

	printf("Doing %d trials...\n", TRIALS);

	for(int i=0 ; i<TRIALS; i++){
		t = rdtsc();
		keygen(pk,sk);
		keygenTime += rdtsc()-t;

		t = rdtsc();
		sign(message, MESSAGE_LENGTH, sk, sig);
		signTime += rdtsc()-t;

		t = rdtsc();
		int ver = verify(message, MESSAGE_LENGTH, pk, sig);;
		verifyTime += rdtsc()-t;

		if(ver < 0){
			printf("Signature invalid! \n");
		}

		// test fast functions
		expand_pk(pk, pk_exp);
		expand_sk(sk, sk_exp);

		t = rdtsc();
		sign_fast(message, MESSAGE_LENGTH, sk, sk_exp, sig);
		signFastTime += rdtsc()-t;

		t = rdtsc();
		ver = verify_fast(message, MESSAGE_LENGTH, pk, pk_exp, sig);
		verifyFastTime += rdtsc()-t;

		if(ver < 0){
			printf("Signature invalid! \n");
		}

	}

	printf("keygen cycles :            %lu \n", keygenTime/TRIALS );
	printf("signing cycles :           %lu \n", signTime/TRIALS );
	printf("verification cycles :      %lu \n", verifyTime/TRIALS );
	printf("fast signing cycles :      %lu \n", signFastTime/TRIALS );
	printf("fast verification cycles : %lu \n", verifyFastTime/TRIALS );

	return 0;
}
