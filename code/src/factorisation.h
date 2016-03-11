#ifndef FACTORISATION_H
#define FACTORISATION_H

#include <gmp.h>
#include <stdbool.h>

signed long int* crible_erat(signed long int n);
mpz_t* factorisation(mpz_t n);
void step_pollard(mpz_t n, mpz_t factor, signed long int B);
void pollard(mpz_t n, signed long int B);
mpz_t* friable(int a, mpz_t n, signed long int B);
void dixon(mpz_t n, signed long int B);
void crible_quadratique(mpz_t b, signed long int B);
#endif 
