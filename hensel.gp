\\ Algorithm to perform hensel lifting on
\\ an integer a and an odd prime p

HENSEL(a, p) =  {
	local(k, sol, t1, t2, l, x1, x2);
	k = 1;

	\\ Begin by using Tonelli to find
	\\ solutions to x^2=5 mod 419
	init_solution = TONELLI(a, p);

	\\ Find t such that x = +-t mod p
	t1 = init_solution[1] % p;
	t2 = init_solution[2] % p;

	\\ Now we want to find solutions
	\\ to x^2 = 5 mod 419^(k+1) for k=1..4
	while(k < a,

		\\ Find x = +-(t + 2*t*l) mod p^k
		\\l = (((5 - t1^2)/p^k) % p)/(2*t1);
		l = FIND_L(a, t1, p, k);

		\\x1 = (t1 + l*p^k) % p^(k+1);
		x1 = FASTEXP((t1 + l*p^k), 1, p^(k+1));
		\\x2 = -(t1 + l*p^k) % p^(k+1);
		x2 = -x1;

		\\ set t to be used for the next iteration
		t1 = (t1 + l*p^k);
		t2 = -(t1 + l*p^k);

		k = k + 1;
	);

	print(x1);
	print(x2);

}

\\ Find l satisfying 2*t*l = (a-t^2)/p^k mod p
FIND_L(a, t, p, k) =  {
	local(l, z);
	z = ((a - t^2)/p^k) % p;
	l = 1;

	while(((2*t*l) % p) != z,
		l = l + 1;
	);

	return(l);
}
