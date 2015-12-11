\\ Row reduce using Gaussian elimination technique described by Bressoud

GAUSS_ELIM(mat, E_length) =  {
	local(mat_local, mat_new, counter, row, col, all_zero, ret_list);

	print("Begin Gaussian elimination on matrix of size: ", matsize(mat)[1], " x ", matsize(mat)[2]);
	mat_local = mat;
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
	print("Gaussian elimination complete! Found ", #ret_list, " zero rows.");

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
