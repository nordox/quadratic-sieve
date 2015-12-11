read("jacobi.gp");
read("tonelli.gp");
read("trial.gp");
read("fast-exp.gp");
read("gcd.gp");
read("gauss.gp");
read("kraitchik.gp");
read("hensel.gp");

QUAD_SIEVE() = {
	local(B, b, S, tonelli_results, r, trial_div_results,
		T, trial_ret, E, E_copy, E_gauss_reduced, kr, x, y, gcd_val);

	n = 135291768612457;	\\ number to be factored (global)

	M = 1900;	\\ sieving interval
	b  = 60;	\\ bound for factor base

/* Choose factor base B = {-1, 2} union {p <= b: (n/p) = 1} */

	B = concat( 2, BUILD_FACTOR_BASE(n, b) );
	B = concat( -1, B );
	\\print("There are ", #B, " primes in the factor base.");

/* Sieve */

	\\ Create zero matrix of size 2M x B
	S = matrix(2*M, #B);

	tonelli_results = List();
	listput(tonelli_results, 0);	\\ for -1
	listput(tonelli_results, 1);	\\ for 2

	\\ Use tonelli's algorithm to solvetp^2 = n mod p
	\\ for all primes in the factor base
	for(i=3, #B,
		listput(tonelli_results, TONELLI(n, B[i]));
	);

	\\ For each row check r = +-tp mod p for all primes in B
	for(i=1, 2*M,
		r = floor(sqrt(n)) - M + i;

		\\ case for 2
		if( r % B[2] == 1,
			S[i, 2] = S[i, 2] + floor( 0.5 + log(B[2]) );
		);

		for(j=3, #B,
			\\ add extra value if the entry satisfies the congruence
			if( (r % B[j]) == tonelli_results[j][1] || (r % B[j]) == tonelli_results[j][2],
				S[i, j] = S[i, j] + floor( 0.5 + log(B[j]) );
			);
		);
	);

/* Threshold and trial division */

	T = 0.7;	\\ threshold a bit less than 2
	trial_div_results = List();

	\\ Attempt trial division of Q(r) for which rows
	\\ of S sum to at least a certain value
	for(i=1, 2*M,
		if(ROW_SUM(S[i,]) >= (0.5*log(n) + log(M) - T*log(b)),
			trial_ret = TRIAL(abs((floor(sqrt(n)) - M + i)^2 - n), 397);

			\\ Store index along with result from trial division. Only keep
			\\ return values that factor into primes within the factor base
			if(trial_ret[1][#trial_ret[1]] <= b,
				listput(trial_div_results, [i, trial_ret]);
			);
		);
	);

/* Gaussian elimination */

	\\ Record full factorizations
	E = matrix(#trial_div_results, #B, i, j, CREATE_EXPONENT_MATRIX(i, j, B, trial_div_results));

	\\ Create identity matrix
	E_id = matid(#trial_div_results);

	\\ Adjoin identity matrix to the exponent matrix
	E_copy = concat(E, E_id);

	\\ Find perfect squares
	E_gauss_reduced = GAUSS_ELIM(E_copy % 2, #E);

	if(#E_gauss_reduced == 0,
		return("Can't find any factors");
	);

/* Kraitchik test */

	\\ go through the zero rows and perform Kraitchik test on the
	\\ second part of the row to try to find a factor
	for(k=1, #E_gauss_reduced,
		\\print(E_gauss_reduced[k]);
		kr = KRAITCHIK(E_gauss_reduced[k], trial_div_results);
		x = kr[1];
		y = kr[2];

		\\print("x: ", x);
		\\print("y: ", y);

		gcd_val = GCD(abs(x - y), n);

		if(gcd_val != 1 && gcd_val != n,
			return([gcd_val, n/gcd_val]);
		);
	);
}

BUILD_FACTOR_BASE(n, b) =  {
	local(factor_base, prime_list, i);

	\\ use the sieve of eratosthenes to
	\\ find all primes under b
	prime_list = SIEVE_ERATOSTHENES(b);

	\\ build a factor base using only those primes
	\\ in prime_list that satisfy (n/p)=1
	factor_base = List();

	for(i=2, #prime_list,
		if(JACOBI(n, prime_list[i]) == 1,
			listput(factor_base, prime_list[i]);
		);
	);

	return(factor_base);
}

\\ sum the columns of the input row
ROW_SUM(cur_row) =  {
	local(row_sum);

	row_sum = 0;

	for(j=1, #cur_row,
		row_sum += cur_row[j];
	);

	return(row_sum);
}

\\ for each row in the full factorization matrix, returns whether
\\ the prime at location (i, j) is a factor in corresponding column in td
CREATE_EXPONENT_MATRIX(row, col, factor_base, td) =  {
	\\ -1 column
	if(col == 1,
		return(1);
	);

	for(i=1, #td[row][2][1],
		if(factor_base[col] == td[row][2][1][i],
			return(td[row][2][2][i]);
		);
	);

	return(0);
}

