#ifndef FACTORISATION_H
#define FACTORISATION_H

#include <gmp.h>
#include <stdbool.h>

unsigned long int* crible_erat();
mpz_t* factorisation(mpz_t n);
void step_pollard(mpz_t n, mpz_t factor);
void pollard(mpz_t n);
mpz_t* friable(int a, mpz_t n, int C);
void dixon(mpz_t n);

#endif 
