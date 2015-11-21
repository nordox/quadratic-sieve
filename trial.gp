\\read("tonelli.gp");	\\ to use MAXIMAL_POWER()

\\ Inputs:
\\	n: integer to be factored
\\	m: max bound

TRIAL(n, m) =  {
	local(i, f, d, p, e);

	i = 0;
	f = n;
	d = 2;
	p = List();
	e = List();

	if(f % d == 0,
		i = i + 1;
		listput(p, d);
		listput(e, MAXIMAL_POWER(f, d)[1]);
		f = MAXIMAL_POWER(f, d)[2];
	);

	d = 3;

	while(d <= m && d^2 <= f,
		if(f % d == 0,
			i = i + 1;
			listput(p, d);
			listput(e, MAXIMAL_POWER(f, d)[1]);
			f = MAXIMAL_POWER(f, d)[2];
		);

		d = d + 2;
	);

	if(f > 1 && d^2 > f,
		i = i + 1;
		listput(p, f);
		listput(e, 1);
		f = 1;
	);

	return([p, e, f]);
}
