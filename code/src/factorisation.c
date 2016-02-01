#include "factorisation.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

#define TAILLE 10000

unsigned long int*
crible_erat()
{
  unsigned long int* tab = calloc(TAILLE, sizeof(unsigned long int));
  unsigned long int* tab2 = malloc(TAILLE*sizeof(unsigned long int));
  int k = 0;
  
  for(int i = 2; i < TAILLE; i++)
    {
      if(tab[i] == 0)
	{
	  tab[i] = 1;
	  tab2[k] = i;
	  k++;
	  for(int j = i*2; j < TAILLE; j=j+i)
	    tab[j] = -1;
	}
    }
  printf("Il y a %d nombres premiers avant %d\n", k, TAILLE);
  free(tab);
  return tab2;
}



// fonction ultra-naïve de factorisation d'entier

void factorisation (mpz_t n)
{
  mpz_t tmp, div, q, r, root;
  //mpz_t cnt;
  unsigned int cnt2 = 0;
  unsigned long int cnt = 0;
  mpz_init(tmp);
  mpz_set(tmp, n);
  mpz_init(div);
  mpz_set_ui(div, 2);
  //mpz_init(cnt);
  mpz_init(root);
  mpz_sqrt(root,n);
  unsigned long int* prime_nbs = crible_erat();
  unsigned long int root_ui = mpz_get_ui(root);
  mpz_t* list_factor = malloc(root_ui*sizeof(mpz_t));
  for(unsigned long int i = 0; i < root_ui; i++)
    mpz_init(list_factor[i]);

  mpz_init(r);
  mpz_init(q);

  while(mpz_cmp_ui(tmp, 1) != 0)
    {
      mpz_tdiv_qr(q, r, tmp, div);
      if (mpz_cmp_ui(r, 0) == 0)
	{
	  mpz_set(tmp, q);
	  mpz_set(list_factor[cnt], div);
	  //mpz_add_ui(cnt, 1);
	  cnt++;
	}
      else
	{
	  if (cnt2 < 1229)
	    {
	      mpz_set_ui(div, prime_nbs[cnt2]);
	      cnt2++;
	    }
	  else
	    mpz_add_ui(div, div, 1);
	}
    }

  for (unsigned long int i = 0; i < cnt; i++)
    gmp_printf("%Zd ", list_factor[i]);
  gmp_printf("\n");
  
  mpz_clear(tmp);
  mpz_clear(div);
  //mpz_clear(cnt);
  for(unsigned long int i = 0; i < root_ui; i++)
    mpz_clear(list_factor[i]);
  free(list_factor);
  mpz_clear(r);
  mpz_clear(q);
  mpz_clear(root);
  free(prime_nbs);
}

/*int* 
factorisation(int n)
{
  int a = n;
  int d = 2;
  int c = 1;
  int* l = malloc(ceil(sqrt(n))*sizeof(int));
  
while(a != 1) // à modifier
    {
      printf("a=%d\n",a);
      printf("d=%d\n",d);
      printf("c=%d\n",c);
      if(a%d == 0)
	{
	  a = a/d;
	  l[c] = d;
	  c++;
	}
      else
	{
	  d++;
	}
    }
  l[0]=c;
  return l;
}
*/


// Problème de segmentation avec l'entier 49392123903


void step_pollard(mpz_t n)
{
  //Initialisation de ce dont on a besoin
  mpz_t factor;
  mpz_init(factor);
  mpz_t pgcd;
  mpz_init(pgcd);
  mpz_t cnt;
  mpz_init(cnt);
  
  //Initialisation de l'algo
  mpz_t rand;
  mpz_init(rand);
  gmp_randstate_t state;
  gmp_randinit_default(state);

  while(mpz_cmp_ui(factor, 0)==0)
    {
    mpz_urandomm(rand, state, n);
    mpz_set_ui(cnt, 2);
    for(int i =0; i < 100; i++)
      {
      mpz_gcd(pgcd,rand - 1,n);
      if(mpz_cmp_ui(pgcd,1)!=0)
	mpz_set(factor, pgcd);
      mpz_powm(rand, rand, cnt, n);
      mpz_add_ui(cnt, cnt, 1);
      }
    }
  gmp_printf("Par pollard: un facteur de %Zd est %Zd \n", n, factor);
  mpz_clear(factor);
  mpz_clear(rand);
  mpz_clear(pgcd);
  mpz_clear(cnt);
  gmp_randclear(state);
}


int
main(int argc, char *argv[])
{
  if(argc < 2)
    {
      printf("Veuillez mettre un nombre à factoriser.\n");
      return EXIT_FAILURE;
    }
  
  mpz_t n;
  mpz_init(n);
  mpz_set_ui(n, atoi(argv[1]));
  factorisation(n);
  step_pollard(n);
  mpz_clear(n);
  return EXIT_SUCCESS;
}

