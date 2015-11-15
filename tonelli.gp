\\ Tonelli's algorithm for finding
\\ solutions to x^2 = a mod p for
\\ integer a and odd prime p

TONELLI(a, p) =  {
	read("jacobi.gp");
	read("fast-exp.gp");
	local(b, s, t, i, c, k, r1, r2, r);
	b = 1;
	k = 1;

	if(JACOBI(a, p) == -1,
		return(Str(a, " is a quadratic nonresidue mod ", p));
	);

	\\ find (b/p) = -1 for some b
	while(JACOBI(b, p) != -1,
		b = b + 1;
	);

	\\ Find max s such that p-1=(2^s)*t
	s = MAXIMAL_POWER(p-1, 2)[1];
	t = MAXIMAL_POWER(p-1, 2)[2];

	i = 2;
	c = (a * (b^2)) % p;
	\\c = (a % p)*(FASTEXP(b, 2, p));

	while(k <= s-1, 
		if(c^(t*(2^(s-k-1))) % p == (p - 1),
			i = i + 2^k;
			c = (c*(b^(2^k))) % p;
		);

		k = k + 1;
	);

	\\ Use fast modular exponentiation
	r1 = FASTEXP(b, i*t/2, p);
	r2 = FASTEXP(a, (t+1)/2, p);
	r = (r1 * r2) % p;

	\\print("r: " r1, " p: ", r2);
	return([r, p-r]);
}

\\ Algorithm to find maximal power e
\\ such that n = (d^e)*f
MAXIMAL_POWER(n, d) =  {
	local(e, f);

	e = 0;
	f = n;

	while(f % d == 0,
		f = f/d;
		e = e + 1;
	);

	return([e, f]);
}
