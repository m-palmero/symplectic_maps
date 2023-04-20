
void boozer_map_eqs(double *x0, double *x1, double *k)
{
    x1[0] = x0[0] - *k * x0[1] * (1 - x0[1]);
    //xn+1=xn-k*yn(1-yn)
	x1[1] = x0[1] + *k * x1[0];
    //yn+1=yn+k*xn+1
}

void boozer_map_bw(double *x0, double *x1, double *k) 
{
    x1[1] = x0[1] - *k * x0[0];
    x1[0] = x0[0] + *k * x1[1] * (1 - x1[1]);
}