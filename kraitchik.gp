KRAITCHIK(row) =  {
	local(index_list, x_list, x, y, in_the_list);

	index_list = List();
	x_list = List();	\\ holds [base, exponent] for each distinct prime in the factorization
	listput(x_list, [2, 0]); \\ initialize list with something
	x = 1;
	y = 1;

	\\ recover index set
	for(i=1, matsize(E_id)[1],
		if(row[i+#E] == 1,
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
		x = ((x%n) * FASTEXP(x_list[i][1], x_list[i][2], n)) % n;
	);

	\\ multiply to get y
	for(i=1, #index_list,
		y = ((y%n) * ((index_list[i][1] + floor(sqrt(n))-M) % n)) % n
	);

	return([x, y]);
}
