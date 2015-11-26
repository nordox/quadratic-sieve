read("jacobi.gp");
read("tonelli.gp");
read("trial.gp");
read("fast-exp.gp");

QUAD_SIEVE() = {
	local(sieving_interval, B, b, S, tp, i, j, r, row_sum, td);

	/* Number to be factored */
	n = 135291768612457;
	\\n = 4999486012441;

	/* Choose sieving interval [-M, M] */
	M = 1000;
	b  = 1229;
	sieving_interval = [-M, M];

	/* Choose factor base B = {-1, 2} union {p <= B: (n/p) = 1} */
	B = concat( 2, BUILD_FACTOR_BASE(n, b) );
	B = concat( -1, B );
	print("There are ", #B, " primes in the factor base.");

/* Sieve */

	\\ Create zero matrix of size 2M x B
	S = matrix(2*M, #B);

	\\ Use tonelli's algorithm to solve
	\\ tp^2 = n mod p for all primes in
	\\ the factor base
	tp = listcreate(2*#B);
	listput(tp, 0);
	listput(tp, 1);

	for(i=3, #B,
		listput(tp, TONELLI(n, B[i]));
	);

	\\ For each row check r = +-tp mod p for r=floor(sqrt(n))-M+i
	\\ for all primes in B
	for(i=1, 2*M,
		r = floor(sqrt(n)) - M + i;

		\\ case for j = 2
		if( r % B[2] == 1,
			S[i, 2] = S[i, 2] + floor( 0.5 + log(B[2]) );
		);

		for(j=3, #B,
			if( (r % B[j]) == tp[j][1] || (r % B[j]) == tp[j][2],
				\\ add extra value
				S[i, j] = S[i, j] + floor( 0.5 + log(B[j]) );
			);
		);
	);

/* Threshold and trial division */

	\\ Choose threshold T a bit less than 2
	T = 1.5;
	td_results = List();

	\\ Attempt trial division of Q(r) for which rows
	\\ of S sum to at least a certain value
	for(i=1, 2*M,
		if(ROW_SUM(i, #B) >= (0.5*log(n) + log(M) - T*log(b)),
			retval = TRIAL(abs((floor(sqrt(n)) - M + i)^2 - n), 10000);

			\\ only output return values that factor
			\\ into primes within the factor base
			if(retval[1][length(retval[1])] <= b,
				listput(td_results, retval);
			);
		);
	);

/* Gaussian elimination */

	\\ Record full factorizations
	E = matrix(length(td_results), length(B), i, j, PR(i, j, B, td_results));

	\\ Make copy of full factorization matrix
	E_copy = E;

	\\ Create identity matrix of length B x B
	E_id = matid(length(td_results));

	\\ Adjoin identity matrix
	E_gauss = concat(E_copy, E_id);

	\\ Find perfect square
	E_gauss_reduced = GAUSS_ELIM(E_gauss % 2, #E_copy, #E_id);

	for(i=1, matsize(E)[1],
		if(E_gauss_reduced[i+#E_copy] == 1,
			print(E[i,]);
		);
	);

/* Kraitchik test */
}

BUILD_FACTOR_BASE(n, b) =  {
	local(factor_base, prime_list, i);

	\\ use the sieve of eratosthenes to
	\\ find all primes under B=1229
	prime_list = SIEVE_ERATOSTHENES(b);

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

ROW_SUM(cur_row, factor_base_size) =  {
	row_sum = 0;

	for(j=1, factor_base_size,
		row_sum = row_sum + S[cur_row, j];
	);

	return(row_sum);
}

\\ for each row in the full factorization matrix, returns whether
\\ the prime at location (i, j) is a factor in corresponding column in td
PR(row, col, factor_base, td) =  {
	\\ First value 1 if Q(r) is negative
	if(col == 1,
		if( (floor(sqrt(n)) - M + row)^2 - n < 0,
			return(1),
			return(0);
		);
	);

	for(i=1, length(td[row][1]),
		if(factor_base[col] == td[row][1][i],
			\\print(td[row][1][i]);
			return(td[row][2][i]);
		);
	);

	return(0);
}

\\ Check if the row has all 0's
IS_ZERO_ROW(row, size, mat) =  {
	for(j=1, size,
		if(mat[i, j] != 0,
			break;
		);

		if(j == size,
			print(mat[i,]);
		);
	);
}

\\ Row reduce using Gaussian elimination technique described in Bressoud
GAUSS_ELIM(mat, E_length, id_length) =  {
	local(mat_local, row, col);
	mat_local = mat;
	row = 1;
	cur_row = 1;
	col = E_length;

	for(k=1, E_length,
		\\ Start at the end of the factorization matrix and find
		\\ the first row where the specified column is 1
		row = FIND_ONE(mat_local, cur_row, col);

		if(row != 0,
			cur_row = row;
			\\ add the row modulo 2 to all succeeding rows
			\\ with a 1 in the column
			for(i=row+1, matsize(mat_local)[1],
				if(mat_local[i, col] == 1,
					mat_local[i,] = ADD_ONE(mat_local[row,], mat_local[i,], #E_copy);

					\\ is the first part of the row all 0's?
					all_zero = 1;
					for(j=1, #E_copy,
						if(mat_local[i, j] == 1,
							all_zero = 0;
						);
					);

					if(all_zero == 1,
						print("h");
						return(mat_local[i,]);
					);
				);
			);

			\\ we are done with the row, so zero it out
			mat_local[row,] = vector(#mat_local, i, 0);
		);

		col--;
	);

	/*for(i=1, matsize(mat_local)[1],
		print(mat_local[i,]);
	);*/

}

FIND_ONE(mat, start_row, col) =  {
	for(i=1, matsize(mat)[1],
		if(mat[i, col] == 1,
			return(i);
		);
	);

	return(0);
}

ADD_ONE(row, row_add, l) =  {
	local(row_new);
	row_new = vector(#row);

	for(i=1, #row,
		if(row[i] == row_add[i],
			row_new[i] = 0,
			row_new[i] = 1;
		);
	);

	return(row_new);
}

\\E_reduced = lift(mattranspose(matimage(mattranspose(Mod(temp, 2)))));
