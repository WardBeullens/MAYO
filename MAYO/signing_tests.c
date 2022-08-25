int test_bilinear_zero_oilspace(){
	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};
	keygen(pk,sk);

	unsigned char P1[P1_BYTES];
	EXPAND(test_randomness_seed, 1, P1, P1_BYTES);
	reduce(P1, P1_BYTES);

	unsigned char oil_space[OIL_SPACE_BYTES] = {0};

	unsigned char P2[P2_BYTES];
	computeP2(oil_space, P1, P2);

	unsigned char bilinear[M*(N-O)*O];
	compute_bilinear_part(P1, oil_space, bilinear);

	for (int i = 0; i < N-O; ++i)
	{
		for (int j = 0; j < O; ++j)
		{
			if(memcmp(P1 + (M*(N-O)*(N-O+1)/2) + M*(i*O + j), bilinear + M*(j*(N-O) + i) , M) ){
				printf("%d %d \n", i, j);
				return -1;
			}
		}
	}

	return 0;
}

int test_bilinear(){
	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};
	keygen(pk,sk);

	unsigned char P1[P1_BYTES] = {0};
	EXPAND(test_randomness_seed, 1, P1, P1_BYTES);
	reduce(P1, P1_BYTES);

	unsigned char oil_space[OIL_SPACE_BYTES] = {0};
	//oil_space[0] = 1;
	sample_oil_space(SK_PRIVATE_SEED(sk),oil_space);

	unsigned char P2[P2_BYTES];
	computeP2(oil_space, P1, P2);

	unsigned char bilinear[M*(N-O)*O];
	compute_bilinear_part(P1, oil_space, bilinear);

	for (int i = 0; i < N-O; ++i)
	{
		for (int j = 0; j < O; ++j)
		{
			unsigned char input[N*K] = {0};
			input[i] = 1;

			unsigned char output1[M];
			evaluateP(input, P1, P2, output1);

			unsigned char oil[O*K] = {0};
			oil[j] = 1;

			add_oil(input, oil, oil_space);

			unsigned char output2[M];
			evaluateP(input, P1, P2, output2);

			negate(output1,M);
			add_vectors(output2,output1,output2);

			if(memcmp(output2, bilinear + M*(j*(N-O)+i), M) != 0){
				printf("%d %d\n", i,j);
				print_vec(output2);
				print_vec((bilinear + M*(j*(N-O)+i)));
				return -1;
			}
		}
	}

	return 0;
}

int test_bilinear2(){
	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};
	keygen(pk,sk);

	unsigned char P1[P1_BYTES] = {0};
	EXPAND(test_randomness_seed, 1, P1, P1_BYTES);
	reduce(P1, P1_BYTES);

	unsigned char oil_space[OIL_SPACE_BYTES] = {0};
	sample_oil_space(SK_PRIVATE_SEED(sk),oil_space);

	unsigned char P2[P2_BYTES];
	computeP2(oil_space, P1, P2);

	unsigned char bilinear[M*(N-O)*O];
	compute_bilinear_part(P1, oil_space, bilinear);

	for (int T = 0; T < 1000; ++T)
	{

		//pick vinegar and oil
		unsigned char inputs[N*K] = {0};
		RAND_bytes(inputs, N-O);
		unsigned char oil[O*K] = {0};
		RAND_bytes(oil,O);

		unsigned char output1[M];
		evaluateP(inputs,P1,P2,output1);

		unsigned char output3[M];

		unsigned char products[(N-O)*O];
		for (int i = 0; i < O; ++i)
		{
			for (int j = 0; j < N-O; ++j)
			{
				products[i*(N-O) + j] = ((uint32_t) oil[i])*((uint32_t) inputs[j]) % PRIME;
			}
		}
		linear_combination(bilinear,products,(N-O)*O,output3);

		// add oil vector
		add_oil(inputs, oil, oil_space);

		unsigned char output2[M];
		evaluateP(inputs,P1,P2,output2);

		negate(output1,M);
		add_vectors(output2,output1,output2);


		if(memcmp(output2,output3,M) != 0){
			print_vec(output2);
			print_vec(output3);
			return -1;
		}

	}

	return 0;
}
