#ifndef FACTORISATION_H
#define FACTORISATION_H

#include <gmp.h>

unsigned long int* crible_erat();
void factorisation(mpz_t n);
void step_pollard(mpz_t n, mpz_t factor);
void pollard(mpz_t n);

#endif 
