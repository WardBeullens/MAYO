
int test_unit_vector_eval(){
	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};
	keygen(pk,sk);

	unsigned char P1[P1_BYTES];
	EXPAND(test_randomness_seed,1,P1,P1_BYTES);

	for (int i = 0; i < N; ++i)
	{
		unsigned char input[N] = {0};
		input[i] = 1;

		unsigned char output[M];
		evaluateP(input, P1, PK_P2(pk), output);

		unsigned char col[M];

		if (i < N-O){
			int col_index = (N-O)*(N-O+1)/2 - (N-O-i)*(N-O-i+1)/2;
			memcpy(col, P1 + M*col_index, M);
		} 
		else {
			int col_index = O*(O+1)/2 - (N-i)*(N-i+1)/2;
			memcpy(col, PK_P2(pk) + M*col_index, M);
		}

		reduce(col,M);

		if( memcmp(output, col ,M) != 0){
			printf("i: %d \n", i);
			print_vec(output);
			print_vec(col);
			return -1;
		}
	}
	return 0;
}

int test_double_unit_vector_eval(){
	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};
	keygen(pk,sk);

	unsigned char P1[P1_BYTES];
	EXPAND(test_randomness_seed,1,P1,P1_BYTES);

	for (int i = 0; i < N; ++i)
	{
		unsigned char input[N] = {0};
		input[i] = 2;

		unsigned char output[M];
		evaluateP(input, P1, PK_P2(pk), output);

		unsigned char col[M];

		if (i < N-O){
			int col_index = (N-O)*(N-O+1)/2 - (N-O-i)*(N-O-i+1)/2;
			memcpy(col, P1 + M*col_index, M);
		} 
		else {
			int col_index = O*(O+1)/2 - (N-i)*(N-i+1)/2;
			memcpy(col, PK_P2(pk) + M*col_index, M);
		}

		// multiply column by 4;
		add_vectors(col,col,col);
		add_vectors(col,col,col);

		if( memcmp(output, col ,M) != 0){
			printf("i: %d \n", i);
			print_vec(output);
			print_vec(col);
			return -1;
		}
	}
	return 0;
}

// test P(2x) = 4P(x) for random x
int test_double_eval(){

	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};
	keygen(pk,sk);

	unsigned char P1[P1_BYTES];
	expand_pk(pk, P1);
	//xorshift_expand(PK_SEED(pk), P1, P1_BYTES);

	for (int i = 0; i < 10; ++i)
	{
		unsigned char in[N] = {0};
		RAND_bytes(in, N);
		reduce(in,N);

		unsigned char in2[N];
		
		for (int j = 0; j < N; ++j)
		{
			in2[j] = (in[j]*2) % PRIME; 
		}

		unsigned char output[M];
		unsigned char output2[M];

		evaluateP(in,  P1, PK_P2(pk), output);
		evaluateP(in2, P1, PK_P2(pk), output2);

		// multiply output by 4
		add_vectors(output,output,output);
		add_vectors(output,output,output);

		if( memcmp(output, output2 ,M) != 0){
			print_vec(output);
			print_vec(output2);
			return -1;
		}
	}
	return 0;
}