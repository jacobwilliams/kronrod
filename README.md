# kronrod

Fortran program for generating Gauss-Kronrod coefficients.

Uses the [MPFUN2020](https://github.com/jacobwilliams/mpfun2020-var1) arbitrary precision Fortran library.

The purpose of this is to use it for the new modernized [QUADPACK](https://github.com/jacobwilliams/quadpack) library.

### Building

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) manifest file is included, so that the application can be compiled and run with FPM. For example:

```
fpm run --profile release
```

### See also

 * Modernized [QUADPACK](https://github.com/jacobwilliams/quadpack) library.
### References

 * The main program is based on this code: https://people.sc.fsu.edu/~jburkardt/f_src/kronrod/kronrod.html, which was modified to use the MPFUN2020 library.
 * Robert Piessens, Maria Branders, "[A Note on the Optimal Addition of Abscissas to Quadrature Formulas of Gauss and Lobatto](https://www.ams.org/journals/mcom/1974-28-125/S0025-5718-1974-0343552-5/S0025-5718-1974-0343552-5.pdf)", Mathematics of Computation, Volume 28, Number 125, January 1974, pages 135-139.
 * David H. Bailey, [MPFUN2020: A new thread-safe arbitrary precision package (Full Documentation)](https://www.davidhbailey.com/dhbpapers/mpfun2020.pdf), January 9, 2022
