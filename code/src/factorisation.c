#include "factorisation.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>

#define TAILLE 10000
#define FRIABLE 10

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
  mpz_add_ui(root,root,1);
  unsigned long int* prime_nbs = crible_erat();
  unsigned long int root_ui = mpz_get_ui(root);
  mpz_t* list_factor = malloc(root_ui*sizeof(mpz_t));
  if (list_factor == NULL)
    {
      fprintf (stderr, "Malloc Problem\n");
      exit(EXIT_FAILURE);
    }
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


void step_pollard(mpz_t n, mpz_t factor)
{
  //Initialisation de ce dont on a besoin
  mpz_t pgcd;
  mpz_init(pgcd);
  mpz_set_ui(pgcd,1);
  mpz_t cnt;
  mpz_init(cnt);
  mpz_t rand_minus_one;
  mpz_init(rand_minus_one);
  int cnt2 = 0;
   mpz_set_ui(factor,0);
  //Initialisation de l'algo
  unsigned long int time_tmp;
  struct timeval t_tmp;
  gettimeofday(&t_tmp, NULL);
  time_tmp = t_tmp.tv_sec*1000000 + t_tmp.tv_usec;
  mpz_t rand;
  mpz_init(rand);  
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, time_tmp);
  mpz_set_ui(rand, 2);
  
  while(mpz_cmp_ui(factor, 0)==0)
      {
	//mpz_urandomm(rand, state, n);
      gmp_printf("Rand: %Zd\n", rand);
      mpz_set_ui(cnt, 2);
      while(cnt2 < FRIABLE && mpz_cmp_ui(pgcd, 1)==0)
	{
	  mpz_sub_ui(rand_minus_one, rand, 1);
	  mpz_gcd(pgcd,rand_minus_one,n);
	  gmp_printf("PGCD %Zd\n", pgcd);
	  if(mpz_cmp_ui(pgcd,1)!=0)
	    {
	      //gmp_printf("%Zd ", factor);
	      mpz_set(factor, pgcd);
	    }
	  mpz_powm(rand, rand, cnt, n);
	  mpz_add_ui(cnt, cnt, 1);
	  gmp_printf("2ème boucle Rand %Zd\n", rand);
	  cnt2++;
	}
    }
  gmp_printf("Par pollard: un facteur de %Zd est %Zd \n", n, factor);
  mpz_clear(rand);
  mpz_clear(pgcd);
  mpz_clear(cnt);
  mpz_clear(rand_minus_one);
  gmp_randclear(state);
  printf("FIN DE STEP POLLARD");
}

void pollard(mpz_t n)
{
  printf("\n POLLARD ALGORITHM \n");
  mpz_t tmp, root, factor, div;
  mpz_init(tmp);
  mpz_set(tmp,n);
  mpz_init(factor);
  mpz_init(root);
  mpz_init(div);
  mpz_sqrt(root,n);
  mpz_add_ui(root,root,1);
  unsigned long int root_ui = mpz_get_ui(root);
  unsigned long int cnt = 0;
  mpz_t* list_factor = malloc(root_ui*sizeof(mpz_t));
  for(unsigned long int i = 0; i < root_ui; i++)
    mpz_init(list_factor[i]);

  while(mpz_cmp_ui(tmp,1)!=0)
    {
      gmp_printf("tmp = %Zd \n factor = %Zd \n", tmp, factor);
      step_pollard(tmp, factor);
      gmp_printf("2 tmp = %Zd \n factor = %Zd \n", tmp, factor);
      mpz_cdiv_q(div, tmp, factor);
      mpz_set(list_factor[cnt], factor);
      mpz_set(tmp, div);
      gmp_printf("\n AFTER tmp = %Zd factor = %Zd div = %Zd\n",tmp, factor, div);
      cnt++;
    }
  
  for (unsigned long int i = 0; i < cnt; i++)
    gmp_printf("%Zd ", list_factor[i]);
  gmp_printf("\n");
  
  for(unsigned long int i = 0; i < root_ui; i++)
    mpz_clear(list_factor[i]);
  free(list_factor);
  mpz_clear(tmp);
  mpz_clear(root);
  mpz_clear(factor);
  }

void brand()
{
  unsigned long int time2;
  struct timeval t1;
  gettimeofday(&t1, NULL);
  time2 = t1.tv_sec*1000000 + t1.tv_usec;
  mpz_t n;
  mpz_init(n);
  mpz_set_ui(n, 124213);
  mpz_t rand;
  mpz_init(rand);  
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, time2);
  mpz_urandomm(rand, state, n);
  gmp_printf("%Zd ", rand);
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
  //brand();
  factorisation(n);
  pollard(n);
  mpz_clear(n);
  return EXIT_SUCCESS;
}

