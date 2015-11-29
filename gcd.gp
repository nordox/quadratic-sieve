GCD(a, b) = {
	local(x, y, temp);
	x = a;
	y = b;

	while(y != 0,
		temp = y;
		y = x % y;
		x = temp
	);

	return(abs(x))
}
