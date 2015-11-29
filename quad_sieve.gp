read("jacobi.gp");
read("tonelli.gp");
read("trial.gp");
read("fast-exp.gp");
read("gcd.gp");

QUAD_SIEVE() = {
	local(sieving_interval, B, b, S, tonelli_results, r, row_sum, trial_div_results,
		T, trial_div_ret, E, E_copy, E_gauss_reduced, x, y, x_list, in_the_list, gcd_val);

	n = 135291768612457;	\\ number to be factored (global)
	\\n = 4999486012441;

	M = 5000;	\\ used for sieving interval
	b  = 500;	\\ bound for factor base
	sieving_interval = [-M, M];

/* Choose factor base B = {-1, 2} union {p <= b: (n/p) = 1} */

	B = concat( 2, BUILD_FACTOR_BASE(n, b) );
	B = concat( -1, B );
	print("There are ", #B, " primes in the factor base.");

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

	T = 1.5;	\\ threshold a bit less than 2
	trial_div_results = List();

	\\ Attempt trial division of Q(r) for which rows
	\\ of S sum to at least a certain value
	for(i=1, 2*M,
		if(ROW_SUM(S[i,]) >= (0.5*log(n) + log(M) - T*log(b)),
			trial_div_ret = TRIAL(abs((floor(sqrt(n)) - M + i)^2 - n), 10000);

			\\ Store index along with result from trial division. Only keep
			\\ return values that factor into primes within the factor base
			if(trial_div_ret[1][#trial_div_ret[1]] <= b,
				listput(trial_div_results, [i, trial_div_ret]);
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

	for(k=1, #E_gauss_reduced,
		index_list = List();
		x_list = List();	\\ holds [base, exponent] for each distinct prime in the factorization
		listput(x_list, [2, 0]); \\ initialize list with something
		x = 1;
		y = 1;

		\\ recover index set
		for(i=1, matsize(E_id)[1],
			if(E_gauss_reduced[k][i+#E] == 1,
				listput(index_list, trial_div_results[i]);
			);
		);

		\\ populate x_list. first go through all entries of index_list
		for(i=1, #index_list,
			print(index_list[i]);
			\\ for each entry, go through all the prime bases
			for(j=1, #index_list[i][2][1],
				\\ go through x_list to check if base is already in the list
				for(z=1, #x_list,
					in_the_list = 0;

					\\ base is in the list
					if(x_list[z][1] == index_list[i][2][1][j],
						\\ update exponent
						x_list[z][2] += index_list[i][2][2][j];
						in_the_list = 1;
						break;
					);
				);

				\\ not in the list, so add it
				if(in_the_list == 0,
					listput(x_list, [index_list[i][2][1][j], index_list[i][2][2][j]]);
				);
			);
		);
		print(x_list);

		\\ go through x_list and divide all exponents by 2
		for(i=1, #x_list,
			x_list[i][2] = x_list[i][2]/2;
		);

		\\ multiply to get x
		for(i=1, #x_list,
			x = (x * x_list[i][1]^x_list[i][2])%n;
		);

		\\ multiply to get y
		for(i=1, #index_list,
			y = ((y%n) * ((index_list[i][1] + floor(sqrt(n))-M) % n));
			y = y % n;
		);

		print("x: ", x);
		print("y: ", y);

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
	factor_base = List(#prime_list);

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
	\\ First value 1 if Q(r) is negative
	if(col == 1,
		/*if( (floor(sqrt(n)) - M + row)^2 - n < 0,
			return(1),
			return(0);
		);*/
		return(1);
	);

	for(i=1, #td[row][2][1],
		if(factor_base[col] == td[row][2][1][i],
			return(td[row][2][2][i]);
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
			print("h");
			print(mat[i,]);
		);
	);
}

\\ Row reduce using Gaussian elimination technique described by Bressoud
GAUSS_ELIM(mat, E_length) =  {
	local(mat_local, mat_new, counter, row, col, all_zero, ret_list);

	mat_local = mat;
	row = 1;
	col = E_length;
	ret_list = List();	\\ holds rows where the first part is all zero

	while(col != 0,
		\\ Start at the end of the factorization matrix and look
		\\ for the first row where the specified column is 1
		row = FIND_ONE(mat_local, col);

		\\ found one
		if(row != 0,
			\\ add the row modulo 2 to all succeeding rows
			\\ with a 1 in the column
			for(i=row+1, matsize(mat_local)[1],
				if(mat_local[i, col] == 1,
					mat_local[i,] = ADD_ROWS(mat_local[row,], mat_local[i,]);

					\\ is the first part of the row all 0's?
					all_zero = 1;

					for(j=1, #E,
						if(mat_local[i, j] != 0,
							all_zero = 0;
						);
					);

					\\ all zero, add to list
					if(all_zero == 1,
						listput(ret_list, mat_local[i,]);
					);
				);
			);

			\\ we are done with the row, so remove it
			/*mat_new = matrix(matsize(mat_local)[1]-1, matsize(mat_local)[2]);
			counter = 1;
			for(i=1, matsize(mat_local)[1],
				if(i != row,
					mat_new[counter,] = mat_local[i,];
					counter++;
				);
			);
			mat_local = mat_new;*/

			mat_local[row,] = vector(#mat_local, i, 0);
		);

		col--;
	);

	return(ret_list);

}

\\ find first row in matrix with a 1 in the specified column
FIND_ONE(mat, col) =  {
	for(i=1, matsize(mat)[1],
		if(mat[i, col] == 1,
			return(i);
		);
	);

	return(0);
}

\\ add row to row_add using xor
ADD_ROWS(row, row_add) =  {
	local(row_new);

	row_new = vector(#row);

	for(o=1, #row,
		if(row[o] != row_add[o],
			row_new[o] = 1;
		);
	);

	return(row_new);
}

