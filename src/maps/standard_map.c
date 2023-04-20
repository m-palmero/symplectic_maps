void standard_map_eqs(double *x0, double *x1, double *k)
{
    x1[1] = x0[1] +  *k * sin (x0[0]);
    constrainAngle(&x1[1]);	
    // if (x1[1] > 2*M_PI)
    // {
    //     x1[1] -= 2*M_PI;
    // }
    //     else if (x1[1] < 0)
    // {
    //     x1[1] += 2*M_PI;
    // }

    x1[0] = x0[0] + x1[1] + M_PI;
    constrainAngle(&x1[0]);
    //x1[0] = fmod(x1[0],2*M_PI);
    // if (x1[0] > 2*M_PI)
    // {
    //     x1[0] -= 2*M_PI;
    // }
    //     else if (x1[0] < 0)
    // {
    //     x1[0] += 2*M_PI;
    // }
}



