QUAD_SIEVE() = {
	read("jacobi.gp");
	read("tonelli.gp");
	read("trial.gp");
	local(M, sieving_interval, B, S, tp, i, j, r, row_sum);

	/* Number to be factored */
	\\n = 135291768612457;
	n = 4999486012441;

	/* Choose sieving interval [-M, M] */
	M = 100;
	sieving_interval = [-M, M];

	/* Choose factor base B = {-1, 2} union {p <= B: (n/p) = 1} */
	B = concat( 2, BUILD_FACTOR_BASE(n) );
	B = concat( -1, B );

/* Sieve */

	\\ Create zero matrix of size 2M x B
	S = matrix(2*M, #B);

	\\ Use tonelli's algorithm to solve
	\\ tp^2 = n mod p for all primes in
	\\ the factor base
	tp = listcreate(2*#B);
	listput(tp, [0, 0]);
	listput(tp, [1, -1]);

	for(i=3, #B,
		listput(tp, TONELLI(n, B[i]));
	);
	print(tp);

	\\ For each row check r = +-tp mod p for r=floor(sqrt(n))-M+i
	\\ for all primes in B
	for(i=1, 2*M,
		r = floor(sqrt(n)) - M + i;
		for(j=2, #B,
			if( (r % B[j]) == tp[j][1] || (r % B[j]) == tp[j][2],
				\\ add extra value
				S[i, j] = S[i, j] + floor( 0.5 + log(B[j]) );
			);
		);
	);
	print(S);

/* Threshold and trial division */

	\\ Choose threshold T a bit less than 2
	T = 1.5;

	\\ Attempt trial division of Q(r) for which rows
	\\ of S sum to at least a certain value
	for(i=1, 2*M,
		row_sum = 0;

		for(j=1, #B,
			row_sum = row_sum + S[i, j];
		);

		if(row_sum >= (0.5*log(n) + log(M) - T*log(1229)),
			print(floor(sqrt(n) - M + i), " ", TRIAL(floor(sqrt(n)) - M + i, 10000));
		);
	);
	

	/* Gaussian elimination */

	/* Kraitchik test */
}

BUILD_FACTOR_BASE(n) =  {
	local(factor_base, prime_list, i);

	\\ use the sieve of eratosthenes to
	\\ find all primes under B=1229
	prime_list = SIEVE_ERATOSTHENES(1229);

	\\ build a factor base using only those primes
	\\ in prime_list that satisfy (n/p)=1
	factor_base = listcreate(#prime_list);
	for(i=2, #prime_list,
		if(JACOBI(n, prime_list[i]) == 1,
			listput(factor_base, prime_list[i]);
		);
	);

	return(factor_base);
}

