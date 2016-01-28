#include "factorisation.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

#define TAILLE 10000

int*
crible_erat()
{
  int* tab = calloc(TAILLE, sizeof(int));
  int* tab2 = malloc(TAILLE*sizeof(int));
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

  free(tab);
  return tab2;
}



// fonction ultra-naïve de factorisation d'entier

void factorisation (mpz_t n)
{
  mpz_t tmp, div, q, r, root;
  //mpz_t cnt;
  unsigned long int cnt = 0;
  mpz_init(tmp);
  mpz_set(tmp, n);
  mpz_init(div);
  mpz_set_ui(div, 2);
  //mpz_init(cnt);
  mpz_init(root);
  mpz_sqrt(root,n);
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
	mpz_add_ui(div, div, 1);
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
int
main(int argc, char *argv[])
{
  mpz_t n;
  mpz_init(n);
  mpz_set_ui(n, atoi(argv[1]));
  int* tab = crible_erat();
  factorisation(n);
  mpz_clear(n);
  free(tab);
  return EXIT_SUCCESS;
}

