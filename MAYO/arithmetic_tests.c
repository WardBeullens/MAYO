
int test_linear_combination_avx(){
	#define TEST_LEN 10000
	for (int i = 0; i < 10; ++i)
	{
		unsigned char matrix[M*TEST_LEN];
		RAND_bytes(matrix,M*TEST_LEN);
		reduce(matrix,M*TEST_LEN);

		unsigned char coefs[TEST_LEN];
		RAND_bytes(coefs,TEST_LEN);
		reduce(coefs,TEST_LEN);

		unsigned char outcome1[M];
		unsigned char outcome2[M];

		_linear_combination(matrix,coefs,TEST_LEN,outcome1);
		_linear_combination_avx(matrix,coefs,TEST_LEN,outcome2);

		if(memcmp(outcome1,outcome2,M) != 0){
			print_vec(outcome1);
			print_vec(outcome2);
			return -1;
		}
	}

	return 0;
}

int test_modular_inverse(){
	for (int i = 1; i < PRIME; ++i)
	{
		if( i*mod_inverse(i) % PRIME != 1 ){
			return -1;
		}
	}
	return 0;
}


int test_sample_solution(){
	unsigned char linear[M*K*O];
	unsigned char rhs[M];
	unsigned char solution[K*O];

	for (int i = 0; i < 1000; ++i)
	{
		while(1){
			RAND_bytes(linear, M*K*O);
			RAND_bytes(rhs, M);
			reduce(linear, M*K*O);
			reduce(rhs, M);

			if( sample_oil(rhs,linear,solution) == 0){
				break;
			}
		}

		unsigned char rhs2[M];
		linear_combination(linear,solution,K*O,rhs2);

		if(memcmp(rhs,rhs2,M) != 0){
			print_vec(rhs);
			print_vec(rhs2);
			return -1;
		} 
	}

	return 0;
}