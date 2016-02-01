#ifndef FACTORISATION_H
#define FACTORISATION_H

#include <gmp.h>

unsigned long int* crible_erat();
void factorisation(mpz_t n);
void factor_pollard(mpz_t n);

#endif 
