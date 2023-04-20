 /** 
    * Simplified Fermi-Ulam map (sfum) (v,phi)
    * v(n+1) = |v(n) - 2eps sin(phi(n))|
    * phi(n+1) = phi(n) + 2/v(n+1)       (mod 2pi)
    */

void sfu_map_eqs(double *x0, double *x1, double *k)
{
    //x1[1] = fabs(x0[1] - (2. * (*k)) * sin(x0[0]));
    x1[1] = x0[1] - (2. * (*k)) * sin(x0[0]);
    //if (x1[1]<0.0) printf("test\n\n\n\n\n\n\n"); 
    //v(n+1) = |v(n) - 2eps sin(phi(n))|
    x1[0] = x0[0] + 2. / x1[1];
    //x1[0] = x0[0] + rand();
    //x1[0] = fmod(x1[0],2. * M_PI);
    constrainAngle(&x1[0]);
    //phi(n+1) = phi(n) + 2/v(n+1)       (mod 2pi)
}

