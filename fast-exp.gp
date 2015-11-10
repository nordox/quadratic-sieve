FASTEXP(b, e, m) = {
	n = 1;

	while(e != 0,
		if(e % 2 == 1,
			n = (n * b) % m;
		);

		e = floor(e/2);
		b = (b * b) % m;
	);

	return(n);
}

TEST() = {
	i = 0;
	s = 100000000000;

	while(i != 10,
		val = FASTEXP(2, s-1, s);

		if(val == 1,
			i = i + 1;
			print(s, " ", val);
		);

		s = s + 1
	);
}
