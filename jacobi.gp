\\ Modification of the Legendre/Jacobi symbol program

JACOBI(a, m) = {
    a = a % m;
    t = 1;
   
    while(a != 0,
        c = 0;

        while(a % 2 == 0,
            a = a/2;
            c = 1 - c;
        );

        if(c == 1,
            if((m % 8 == 3) || (m % 8 == 5),
                t = -t;
            );
        );

        if((a % 4 == 3) && (m % 4 == 3),
            t = -t;
        );

        temp = m;
        m = a;
        a = temp;
        a = a % m;
    );

    if(m == 1,
        return(t);
    );

    return(0);
}

SIEVE_ERATOSTHENES(n) = {
    local(l, h, p, j, i, k);

    l = vector(n);
    h = vector(n);

    for(i=2, n,
        l[i] = 1;
    );

    p = 2;

    while(p^2 <= n,
        j = p + p;

        while(j <= n,
            l[j] = 0;
            j = j + p;
        );

        p = p + 1;

        while(l[p] == 0,
            p = p + 1;
        );
    );

    h[1] = 2;
    j = 2;

    for(p=3, n,
        if(l[p] == 1,
            h[j] = p;
        );

        j = j + l[p];
    );

	\\ count how many primes there are
	i = 1;
	while(h[i] != 0,
		i = i + 1;
	);

	\\ setup a matrix with only the primes (no 0s)
	retval = vector(i-1);
	for(k=1, i-1,
		retval[k] = h[k];
	);

    return(retval);
}

TEST() = {
    n = 3541905253352059459794529;
    maxN = 100;
    res = SIEVE_ERATOSTHENES(maxN);
    count = 0;
   
    z = 1;

    while(count < 100,
        while(res[z] != 0,
            if(JACOBI(n, res[z]) == 1,
                print(res[z]);
                count = count + 1;
            );
            z = z + 1;
        );
        maxN = maxN + 100;
        res = SIEVE_ERATOSTHENES(maxN);
    );
}
