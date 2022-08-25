#include "arithmetic.h"
#include "immintrin.h"
#include <stdalign.h>


void add_vectors(const unsigned char *v1, const unsigned char *v2, unsigned char *out){
	for (int i = 0; i < M; ++i)
	{
		out[i] = (((uint16_t) v1[i]) + ((uint16_t) v2[i])) % PRIME;
	}
}

void negate(unsigned char* v, int len){
	for (int i = 0; i < len; ++i)
	{
		v[i] = (PRIME - v[i]) % PRIME;
	}
}

void reduce(unsigned char* v, int len){
	for (int i = 0; i < len; ++i)
	{
		v[i] = v[i] % PRIME;
	}
}

void scalar_multiply(unsigned char* v, unsigned char a, int len){
	for (int i = 0; i < len; ++i)
	{
		v[i] = (((uint16_t)v[i])*a) % PRIME;
	}
}

#define REDUCTION_CT _mm256_set_epi16(2114, 2114, 2114, 2114, 2114, 2114, 2114, 2114, 2114, 2114, 2114, 2114, 2114, 2114, 2114, 2114)
#define PRIMEX16 _mm256_set_epi16(PRIME, PRIME, PRIME, PRIME, PRIME, PRIME, PRIME, PRIME, PRIME, PRIME, PRIME, PRIME, PRIME, PRIME, PRIME, PRIME)

void _linear_combination_avx(const unsigned char* vecs, const unsigned char* coeffs, int len, unsigned char *out){
	__m256i accumulators[4] = {0};
	//const __m256i zero = {0};

	for (int i = 0; i < len; i+=2)
	{
		unsigned char coef1 = coeffs[i]%PRIME;
		unsigned char coef2 = 0;
		if(i+1<len){
			coef2 = coeffs[i+1]%PRIME;
		}

		__m256i coef  = _mm256_unpacklo_epi8(_mm256_set1_epi8(coef1),_mm256_set1_epi8(coef2));
		
		// do first 32 entries of vectors
		__m256i v11 = _mm256_loadu_si256( (__m256i*)(vecs +     i*M));
		__m256i v21 = _mm256_loadu_si256( (__m256i*)(vecs +     i*M + 32));
		__m256i v12 = _mm256_loadu_si256( (__m256i*)(vecs + (i+1)*M));
		__m256i v22 = _mm256_loadu_si256( (__m256i*)(vecs + (i+1)*M + 32));

		__m256i v1lo = _mm256_unpacklo_epi8(v11, v12); // 0-7,  16-23
		__m256i v1hi = _mm256_unpackhi_epi8(v11, v12); // 8-15, 24-31
		__m256i v2lo = _mm256_unpacklo_epi8(v21, v22); // 32-39, 48-55
		__m256i v2hi = _mm256_unpackhi_epi8(v21, v22); // 40-47, 56-63

		__m256i p1lo = _mm256_maddubs_epi16(v1lo,coef); // 0-7,  16-23
		__m256i p1hi = _mm256_maddubs_epi16(v1hi,coef); // 8-15, 24-31
		__m256i p2lo = _mm256_maddubs_epi16(v2lo,coef); // 32-39, 48-55
		__m256i p2hi = _mm256_maddubs_epi16(v2hi,coef); // 40-47, 56-63

		accumulators[0] = _mm256_add_epi16(accumulators[0], p1lo); 
		accumulators[1] = _mm256_add_epi16(accumulators[1], p1hi); 
		accumulators[2] = _mm256_add_epi32(accumulators[2], p2lo); 
		accumulators[3] = _mm256_add_epi32(accumulators[3], p2hi); 

		if (i%70 == 68){
			//reduce accumulators

			__m256i m0 = _mm256_mulhi_epu16(REDUCTION_CT, accumulators[0]);
			__m256i m1 = _mm256_mulhi_epu16(REDUCTION_CT, accumulators[1]);
			__m256i m2 = _mm256_mulhi_epu16(REDUCTION_CT, accumulators[2]);
			__m256i m3 = _mm256_mulhi_epu16(REDUCTION_CT, accumulators[3]);

			accumulators[0] = _mm256_sub_epi16(accumulators[0],_mm256_mullo_epi16(m0, PRIMEX16));
			accumulators[1] = _mm256_sub_epi16(accumulators[1],_mm256_mullo_epi16(m1, PRIMEX16));
			accumulators[2] = _mm256_sub_epi16(accumulators[2],_mm256_mullo_epi16(m2, PRIMEX16));
			accumulators[3] = _mm256_sub_epi16(accumulators[3],_mm256_mullo_epi16(m3, PRIMEX16));
		}
	}

	unsigned char pre_out[64];

	uint16_t *buf = (uint16_t*) accumulators;

	for (int j = 0; j < 8; ++j)
	{
		pre_out[   j] = (unsigned char) (buf[   j] % PRIME);
		pre_out[8 +j] = (unsigned char) (buf[16+j] % PRIME);
		pre_out[16+j] = (unsigned char) (buf[8 +j] % PRIME);
		pre_out[24+j] = (unsigned char) (buf[24+j] % PRIME);
		pre_out[32+j] = (unsigned char) (buf[32+j] % PRIME);
		pre_out[40+j] = (unsigned char) (buf[48+j] % PRIME);
		pre_out[48+j] = (unsigned char) (buf[40+j] % PRIME);
		pre_out[56+j] = (unsigned char) (buf[56+j] % PRIME);
	}

	memcpy(out,pre_out,M);
}

void _linear_combination(const unsigned char* vecs, const unsigned char* coeffs, int len, unsigned char *out){
	uint32_t accumulators[M] = {0};
	for (int i = 0; i < len; ++i)
	{
		for(int j=0; j< M; j++){
			accumulators[j] += ((uint32_t) vecs[i*M + j]) * ((uint32_t) coeffs[i]);
		}
	}

	for (int i = 0; i < M; ++i)
	{
		out[i] = (unsigned char) (accumulators[i] % PRIME);
	}
}

void zero_out(unsigned char *out){
	for (int i = 0; i < M; ++i)
	{
		out[i] = 0;
	}
}

void evaluateP(const unsigned char *input, const unsigned char *P1, const unsigned char *P2, unsigned char *output){
	unsigned char products[MONOMIALS];
//TIC
	int counter = 0;
	// vinegar x vinegar
	for (int i = 0; i < N-O; ++i)
	{
		for (int j = i; j < N-O; ++j)
		{
			products[counter++] = (((uint32_t)input[i])*((uint32_t) input[j])) % PRIME;		
		}
	}

	// vinegar x oil
	for (int i = 0; i < N-O; ++i)
	{
		for (int j = N-O; j < N; ++j)
		{
			products[counter++] = (((uint32_t)input[i])*((uint32_t) input[j])) % PRIME;		
		}
	}

	// oil x oil
	for (int i = N-O; i < N; ++i)
	{
		for (int j = i; j < N; ++j)
		{
			products[counter++] = (((uint32_t)input[i])*((uint32_t) input[j])) % PRIME;		
		}
	}
	unsigned char part1[M] = {0};
	unsigned char part2[M] = {0};
//TOC(products)
	linear_combination(P1,products              , P1MONOMIALS, part1);
	linear_combination(P2,products + P1MONOMIALS, P2MONOMIALS, part2);
//TOC(lin_comb)
	add_vectors(part1,part2,output);
}

void evaluateP_vinegar(const unsigned char *input, const unsigned char *P1, unsigned char *out){

	unsigned char products[(N-O)*(N-O+1)/2];

	int counter = 0;
	for (int i = 0; i < N-O; ++i)
	{
		for (int j = i; j < N-O; ++j)
		{
			products[counter++] = (((uint32_t)input[i])*((uint32_t) input[j])) % PRIME;	
		}
	}

	linear_combination(P1,products, (N-O)*(N-O+1)/2, out);
}

void print_system(uint32_t *matrix){
	printf("\n");
	for (int i = 0; i < M; ++i)
	{
		printf("%5d %5d", matrix[i*64], matrix[i*64 + 1]);
		for (int j = 1; j < O*K; ++j)
		{
			printf("%s", matrix[i*64 + j] == 0? "_" : "X" );
			// printf("%3d ", matrix[i*64 + j]);
		}
		printf("| %3d \n", matrix[i*64 + O*K]);
	}
	printf("\n");
}

void swap_row(uint32_t *matrix, int a, int b){
	__m256i *rowa = (__m256i*) (matrix + 64*a);
	__m256i *rowb = (__m256i*) (matrix + 64*b);
	for (int i = 0; i < 8; ++i)
	{
		__m256i tmp = rowa[i];
		rowa[i] = rowb[i];
		rowb[i] = tmp;
	}
}

// TODO: make less retarded version if needed.
uint32_t mod_inverse(uint32_t a){
	uint16_t c = 1;
	for (int i = 0; i < PRIME-2; ++i)
	{
		c = (c*a) % PRIME;
	}
	return c;
}

void scale(uint32_t *matrix, int row, uint32_t a){
	for (int i = 0; i < (K*O)+1; ++i)
	{
		matrix[row*64 + i] = (((matrix[row*64 + i]%PRIME) * a) %PRIME);
	}
}

void row_op(uint32_t *matrix, int s, int d, int coef){
	for (int i = s; i < K*O+1; ++i)
	{
		matrix[d*64 + i] += (matrix[s*64 + i] * coef);
	}
}

void gauss_reduction(uint32_t *matrix){
	int row = 0;
	int col = 0;

	while (1) {
		//find next nonzero entry
		int find_row = row;
		while (matrix[find_row*64 + col]%PRIME == 0){
			matrix[find_row*64 + col] = 0;
			find_row++;
			if(find_row == M){
				col ++;
				find_row = row;
				if(col == K*O){
					return;
				}
			}
		}

		if (find_row != row){
			swap_row(matrix, row, find_row);
		}

		scale(matrix, row, mod_inverse(matrix[row*64 + col] % PRIME));

		for (int i = find_row+1; i < M; ++i)
		{
			row_op(matrix, row, i, (PRIME - (matrix[i*64 + col]%PRIME)) % PRIME );
		}

		row ++;
		col ++;
		if(row == M){
			return;
		}
	}
}

unsigned char random_element(){
	return 1;
	unsigned char r = PRIME;
	while( r >= PRIME){
		RAND_bytes(&r, 1);
	}
	return r;
}

int sample_oil(const unsigned char *rhs,const unsigned char *linear, unsigned char *solution){
//TIC
	#if K*O > 8*8 - 1
		Error: K*O > 8*8 -1 not supported
	#endif
	__m256i aug_matrix_avx[M*8];
	uint32_t *aug_matrix = (uint32_t*) aug_matrix_avx;
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < K*O; ++j)
		{
			aug_matrix[i*64 + j] = linear[j*M + i];
		}
		aug_matrix[i*64 + K*O] = rhs[i];
	}
//TOC(copy)
	gauss_reduction(aug_matrix);
//TOC(gauss_reduction)
	int col = (K*O);
	int row = M-1;
	while (row >= 0){
		int col2 = 0;
		while (aug_matrix[row*64 + col2]%PRIME == 0){
			col2 ++;
			if(col2 == K*O+1){
				break;
			}
		}

		if(col2 == K*O+1){
			row--;
			continue;
		}

		if(col2 == K*O){
			return -1;
		}

		while(col > col2 +1){
			col--;
			// choose solution entry at random
			solution[col] = random_element();

			for (int i = 0; i < M; ++i)
			{
				aug_matrix[i*64 + O*K] = (((uint16_t) aug_matrix[i*64 + O*K]) + ((uint16_t)(PRIME-solution[col]))*aug_matrix[i*64 + col]) % PRIME;
			}
		}

		col--;
		solution[col] = aug_matrix[row*64 + K*O];
		for (int i = 0; i < M; ++i)
		{
			aug_matrix[i*64 + O*K] = (((uint16_t) aug_matrix[i*64 + O*K]) + ((uint16_t)(PRIME-solution[col]))*aug_matrix[i*64 + col]) % PRIME;
		}

		row --;
	}
//TOC(substitution)
	return 0;
}
