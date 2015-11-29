read("jacobi.gp");
read("tonelli.gp");
read("trial.gp");
read("fast-exp.gp");
read("gcd.gp");

QUAD_SIEVE() = {
	local(sieving_interval, B, b, S, tp, i, j, r, row_sum, td_results, T, retval, x, y, temp, g);

	/* Number to be factored */
	n = 135291768612457;
	\\n = 4999486012441;

	/* Choose sieving interval [-M, M] */
	M = 5000;
	b  = 500;
	sieving_interval = [-M, M];

	/* Choose factor base B = {-1, 2} union {p <= B: (n/p) = 1} */
	B = concat( 2, BUILD_FACTOR_BASE(n, b) );
	B = concat( -1, B );
	print("There are ", #B, " primes in the factor base.");
	print(B);

/* Sieve */

	\\ Create zero matrix of size 2M x B
	S = matrix(2*M, #B);

	\\ Use tonelli's algorithm to solve
	\\ tp^2 = n mod p for all primes in
	\\ the factor base
	tp = listcreate();
	listput(tp, 0);
	listput(tp, 1);

	for(i=3, #B,
		listput(tp, TONELLI(n, B[i]));
	);
	print("tp: ", " ", tp);

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
		if(ROW_SUM(S[i,]) >= (0.5*log(n) + log(M) - T*log(b)),
			retval = TRIAL(abs((floor(sqrt(n)) - M + i)^2 - n), 10000);

			\\ only output return values that factor
			\\ into primes within the factor base
			if(retval[1][#retval[1]] <= b,
				print("num: ", abs((floor(sqrt(n)) - M + i)^2 - n));
				print("div: ", retval);
				listput(td_results, [i, retval]);
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

	if(#E_gauss_reduced == 0,
		return("none");
	);

/* Kraitchik test */

	for(k=1, #E_gauss_reduced,
		index_list = listcreate();

		\\ recover index set
		for(i=1, matsize(E_id)[1],
			if(E_gauss_reduced[k][i+#E_copy] == 1,
				listput(index_list, td_results[i]);
			);
		);

		x = 1;
		y = 1;

		ind_list = matrix(#index_list, #B, i, j, PR(i, j, B, index_list));
		ind_list = lift(mattranspose(matimage(mattranspose(Mod(ind_list, 2)))));

		i_list = List();
		for(i=1, matsize(ind_list)[1],
			for(j=1, matsize(ind_list)[2],
				if(ind_list[i, j] == 1,
					x = x * B[j];
				);
			);
		);
x = x % n;

	x = 1;
		for(i=1, #index_list,
			temp = 1;

			print(index_list[i]);
			for(j=1, #index_list[i][2][1],
				\\x = (x%n) * FASTEXP(index_list[i][2][1][j], index_list[i][2][2][j], n);
				x = x * (index_list[i][2][1][j]^index_list[i][2][2][j]);
			);

			x = (x * temp);
			y = ((y%n) * ((index_list[i][1] + floor(sqrt(n))-M) % n));
		);

		x = sqrt(x) % n;
		\\y = sqrt(y) % n;
		y = y % n;
		print("x: ", x);
		print("y: ", y);

		g = GCD(abs(x - y), n);
		if(g != 1 && g != n,
			return([g, n/g]);
		);
		print();
	);
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

ROW_SUM(cur_row) =  {
	local(row_sum);

	row_sum = 0;

	for(j=1, #cur_row,
		row_sum = row_sum + cur_row[j];
	);

	return(row_sum);
}

\\ for each row in the full factorization matrix, returns whether
\\ the prime at location (i, j) is a factor in corresponding column in td
PR(row, col, factor_base, td) =  {
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
GAUSS_ELIM(mat, E_length, id_length) =  {
	local(mat_local, mat_new, counter, row, col, all_zero, ret_list);

	mat_local = mat;
	row = 1;
	col = E_length;
	ret_list = listcreate();

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
					for(j=1, #E_copy,
						if(mat_local[i, j] != 0,
							all_zero = 0;
						);
					);

					if(all_zero == 1,
						listput(ret_list, mat_local[i,]);
					);
				);
			);

			\\ we are done with the row, so remove it
			mat_new = matrix(matsize(mat_local)[1]-1, matsize(mat_local)[2]);
			counter = 1;
			for(i=1, matsize(mat_local)[1],
				if(i != row,
					mat_new[counter,] = mat_local[i,];
					counter++;
				);
			);
			mat_local = mat_new;

			\\mat_local[row,] = vector(#mat_local, i, 0);
		);

		col--;
	);

	/*for(i=1, matsize(mat_local)[1],
		print(mat_local[i,]);
	);*/

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

\\E_reduced = lift(mattranspose(matimage(mattranspose(Mod(temp, 2)))));
