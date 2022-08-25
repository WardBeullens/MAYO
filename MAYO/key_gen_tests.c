

int test_zero_oil_space(){

	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};
	keygen(pk,sk);

	unsigned char P1[P1_BYTES];
	EXPAND(test_randomness_seed, 1, P1, P1_BYTES);

	unsigned char oil_space[OIL_SPACE_BYTES] = {0};
	unsigned char P2[P2_BYTES];

	computeP2(oil_space, P1, P2);

	int non_zero = 0;
	for (int i = 0; i < P2_BYTES; ++i)
	{
		if(P2[i] != 0){
			non_zero = -1;
		}
	}
	return non_zero;
}

int test_almost_zero_oil_space(){

	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};
	keygen(pk,sk);

	unsigned char P1[P1_BYTES];
	EXPAND(test_randomness_seed, 1, P1, P1_BYTES);

	unsigned char oil_space[OIL_SPACE_BYTES] = {0};
	oil_space[0] = 1;
	unsigned char P2[P2_BYTES];

	computeP2(oil_space, P1, P2);

	negate(P2, P2_BYTES);
	reduce(P1, P1_BYTES);
	add_vectors(P1 + (N-O)*(N-O+1)/2*M, P1, P1 + (N-O)*(N-O+1)/2*M);

	if( memcmp(P2 , P1 + (N-O)*(N-O+1)/2*M, O*M) != 0){
		for (int i = 0; i < O*M; ++i)
		{
			printf("%3d: %3d %3d\n", i,  P2[i], P1[(N-O)*M + i]);
		}
		return -1;
	}

	return 0;
}

int test_eval_zero_oil_vector(){

	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};
	keygen(pk,sk);

	unsigned char P1[P1_BYTES];
	EXPAND(test_randomness_seed, 1, P1, P1_BYTES);

	unsigned char oil_space[OIL_SPACE_BYTES] = {0};
	unsigned char P2[P2_BYTES];

	computeP2(oil_space, P1, P2);

	unsigned char input[N] = {0};
	memcpy(input, oil_space, N-O);
	
	for (int i = 0; i < O; ++i)
	{
		input[N-O+i] = rand()%PRIME;
	}

	unsigned char output[M];
	evaluateP(input, P1, P2, output);

	for (int i = 0; i < M; ++i)
	{
		if(output[i] != 0){
			print_vec(output);
			return -1;
		}
	}

	return 0;
}

int test_eval_oil_vector(){

	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};
	keygen(pk,sk);

	unsigned char P1[P1_BYTES];
	expand_pk(PK_SEED(pk),P1);

	unsigned char oil_space[OIL_SPACE_BYTES];
	sample_oil_space(SK_PRIVATE_SEED(sk),oil_space);

	unsigned char input[N] = {0};
	
	memcpy(input, oil_space, N-O);
	input[N-O] = 1;

	unsigned char output[M];
	evaluateP(input, P1, PK_P2(pk), output);

	for (int i = 0; i < M; ++i)
	{
		if(output[i] != 0){
			print_vec(output);
			return -1;
		}
	}

	return 0;
}