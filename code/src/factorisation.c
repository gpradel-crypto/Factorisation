#include "factorisation.h"

#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

// Renvoie la liste des nombres premiers inférieurs à n (avec -1 en premier)
signed long int *
crible_erat (signed long int n)
{
//Initialisation des tableaux que l'on va utiliser
  signed long int *tab = calloc ((n + 2), sizeof (signed long int));
  if (tab == NULL)
    {
      printf ("Dans crible_erat : calloc de tab impossible.\n");
      exit (EXIT_FAILURE);
    }

  signed long int *tab2 = calloc (((n + 3) / 2), sizeof (signed long int));
  if (tab2 == NULL)
    {
      printf ("Dans crible_erat : calloc de tab2 impossible.\n");
      exit (EXIT_FAILURE);
    }

//Algorithme
  int k = 1;
  tab2[0] = -1;

  for (int i = 2; i <= n; i++)
    {
      if (tab[i] == 0)
	{
	  tab[i] = 1;
	  tab2[k] = i;
	  k++;
	  for (int j = i * 2; j <= n; j = j + i)
	    tab[j] = -1;
	}
    }
  free (tab);
  return tab2;
}


//Calcule un facteur (pas forcément premier) de n 
void
step_pollard (mpz_t n, mpz_t factor, signed long int B)
{
  //Initialisation de ce dont on a besoin
  mpz_t pgcd;
  mpz_init (pgcd);
  mpz_set_ui (pgcd, 1);
  mpz_t rest;
  mpz_init (rest);
  mpz_t cnt, tmp_prime;
  mpz_init (cnt);
  mpz_init (tmp_prime);
  mpz_t rand_minus_one;
  mpz_init (rand_minus_one);
  int cnt2 = 0;
  mpz_set_ui (factor, 0);

  //Initialisation de l'algorithme
  unsigned long int time_tmp;
  struct timeval t_tmp;
  gettimeofday (&t_tmp, NULL);
  time_tmp = t_tmp.tv_sec * 1000000 + t_tmp.tv_usec;
  mpz_t rand;
  mpz_init (rand);
  gmp_randstate_t state;
  gmp_randinit_default (state);
  gmp_randseed_ui (state, time_tmp);

  mpz_set_ui (rand, 2);
  mpz_cdiv_r_ui (rest, n, 2);
  if (mpz_cmp_ui (rest, 0) == 0)
    mpz_set_ui (rand, 3);

  //Algorithme
  while (cnt2 < B && mpz_cmp_ui (factor, 0) == 0)
    {
      if (mpz_cmp_ui (rand, 0) == 0 || mpz_cmp_ui (rand, 1) == 0)
	mpz_urandomm (rand, state, n);
      mpz_set_ui (cnt, 2);
      while (cnt2 < B && mpz_cmp_ui (pgcd, 1) == 0)
	{
	  if (mpz_cmp_ui (rand, 1) == 0 || mpz_cmp_ui (rand, 0) == 0)
	    mpz_urandomm (rand, state, n);
	  mpz_sub_ui (rand_minus_one, rand, 1);
	  mpz_gcd (pgcd, rand_minus_one, n);
	  if (mpz_cmp_ui (pgcd, 1) != 0)
	    mpz_set (factor, pgcd);
	  mpz_powm (rand, rand, cnt, n);
	  mpz_nextprime (cnt, tmp_prime);
	  mpz_set (tmp_prime, cnt);
	  cnt2++;
	}
    }

  //Aucun facteur trouvé. Sûrement premier ou échec de l'algorithme. On renvoie le même nombre.
  if (mpz_cmp_ui (factor, 0) == 0)
    mpz_set (factor, n);

  //On libère la mémoire
  mpz_clear (rand);
  mpz_clear (pgcd);
  mpz_clear (rest);
  mpz_clear (cnt);
  mpz_clear (tmp_prime);
  mpz_clear (rand_minus_one);
  gmp_randclear (state);
}

//Factorise n petit à petit en utilisant step_pollard
void
pollard (mpz_t n, signed long int B)
{
  //Initialisation
  mpz_t tmp, root, factor, div, check;
  mpz_init (tmp);
  mpz_set (tmp, n);
  mpz_init (factor);
  mpz_init (root);
  mpz_init (check);
  mpz_init (div);
  mpz_sqrt (root, n);
  mpz_add_ui (root, root, 1);
  unsigned long int cnt = 0;
  int cnt2 = 0;
  mpz_t *list_factor = malloc (1000 * sizeof (mpz_t));
  if (list_factor == NULL)
    {
      printf ("Dans pollard : malloc de list_factor impossible.\n");
      exit (EXIT_FAILURE);
    }
  for (unsigned long int i = 0; i < 1000; i++)
    mpz_init (list_factor[i]);

  //Utilisation de step_pollard pour factoriser
  while (mpz_cmp_ui (tmp, 1) != 0)
    {
      step_pollard (tmp, factor, B);
      while (cnt2 < B && mpz_probab_prime_p (factor, 47) == 0)
	{
	  mpz_set (check, factor);
	  step_pollard (check, factor, B);
	  cnt2++;
	}
      mpz_cdiv_q (div, tmp, factor);
      mpz_set (list_factor[cnt], factor);
      mpz_set (tmp, div);
      cnt++;

    }

  //Ici nous avons la liste de tous les facteurs. Maintenant on les regroupe
  int nb_divisors = 1;
  unsigned long int divisor = mpz_get_ui (list_factor[0]);
  for (unsigned int i = 0; i < cnt; i++)
    if (mpz_cmp_ui (list_factor[i], divisor) != 0)
      {
	nb_divisors++;
	divisor = mpz_get_ui (list_factor[i]);
      }
  divisor = mpz_get_ui (list_factor[0]);
  unsigned int cnt_tmp = 0;
  int cnt_tmp2 = 0;
  int *coef = calloc (nb_divisors, sizeof (int));
  if (coef == NULL)
    {
      printf ("Dans pollard : calloc de coef impossible.\n");
      exit (EXIT_FAILURE);
    }
  while (cnt_tmp < cnt)
    {
      if (mpz_cmp_ui (list_factor[cnt_tmp], divisor) == 0)
	coef[cnt_tmp2]++;
      else
	{
	  cnt_tmp2++;
	  coef[cnt_tmp2]++;
	}
      cnt_tmp++;
      divisor = mpz_get_ui (list_factor[cnt_tmp - 1]);
    }

  //Affichage des facteurs avec leur multiplicité
  cnt_tmp2 = 0;
  for (int i = 0; i < nb_divisors; i++)
    {
      cnt_tmp2 += coef[i];
      gmp_printf ("[%Zd, %d]; ", list_factor[cnt_tmp2 - 1], coef[i]);
    }
  gmp_printf ("\n");

  //On libère la mémoire
  for (unsigned long int i = 0; i < 1000; i++)
    mpz_clear (list_factor[i]);
  free (list_factor);
  mpz_clear (tmp);
  mpz_clear (root);
  mpz_clear (factor);
  mpz_clear (check);
  mpz_clear (div);
  free (coef);
}

// Si a = 0, imprime les facteurs premiers de n, imprime et renvoie leurs multiplicités (factorisation naïve)
// Si a = 1, renvoie les multiplicités des facteurs de n inférieurs à B
mpz_t *
friable (int a, mpz_t n, signed long int B)
{
  //Initialisation
  int nb_primes = 0;
  int k = 1;
  signed long int m = mpz_get_si (n);
  signed long int *prime_nbs = crible_erat (B);
  while (prime_nbs[nb_primes] <= B && prime_nbs[nb_primes] != 0)
    nb_primes++;

  mpz_t *smooth_list = malloc ((nb_primes + 1) * sizeof (mpz_t));
  if (smooth_list == NULL)
    {
      printf ("Dans friable : malloc de smooth_list impossible.\n");
      exit (EXIT_FAILURE);
    }
  for (int l = 0; l < nb_primes + 1; l++)
    mpz_init (smooth_list[l]);

  if (mpz_cmp_ui (n, 0) == 0)
    return smooth_list;

  if (mpz_cmp_ui (n, 0) < 0)
    {
      mpz_add_ui (smooth_list[k], smooth_list[k], 1);
      m = abs (m);
    }

  //Division petit à petit de n en utilisant l'amélioration via le crible d'Eratosthène
  while (k < nb_primes)
    {
      if (m % prime_nbs[k] == 0)
	{
	  m = m / prime_nbs[k];
	  mpz_add_ui (smooth_list[k + 1], smooth_list[k + 1], 1);
	}
      else
	k++;
    }
  //Si le nombre est B-friable
  if (m == 1)
    {
      mpz_set_ui (smooth_list[0], 1);
      if (a == 0)
	{
	  //printf ("%ld-FRIABLE\n", B);
	  for (int l = 1; l < nb_primes + 1; l++)
	    if (mpz_cmp_ui (smooth_list[l], 0) != 0)
	      gmp_printf ("[%d, %Zd] ", prime_nbs[l - 1], smooth_list[l]);
	  printf ("\n");
	}
    }
  //S'il n'est pas B-friable, on fait vraiment naïvement
  else
    {
      mpz_set_ui (smooth_list[0], 0);
      if (a == 0)
	{
	  mpz_t **list_factor = malloc (1000 * sizeof (mpz_t *));
	  if (list_factor == NULL)
	    {
	      printf ("Dans friable : malloc de list_factor impossible.\n");
	      exit (EXIT_FAILURE);
	    }
	  for (int l = 0; l < 1000; l++)
	    {
	      list_factor[l] = malloc (2 * sizeof (mpz_t));
	      if (list_factor[l] == NULL)
		{
		  printf
		    ("Dans friable : malloc de list_factor[l] impossible.\n");
		  exit (EXIT_FAILURE);
		}
	      for (unsigned long int i = 0; i < 2; i++)
		mpz_init (list_factor[l][i]);
	    }
	  unsigned long int cnt = prime_nbs[nb_primes - 1] + 1;
	  int cnt2 = 0;
	  while (m != 1)
	    {
	      if (m % cnt == 0)
		{
		  m = m / cnt;
		  mpz_set_ui (list_factor[cnt2][0], cnt);
		  mpz_add_ui (list_factor[cnt2][1], list_factor[cnt2][1], 1);
		  if (m % cnt != 0)
		    cnt2++;
		}
	      else
		cnt++;
	    }

	  //On libère la mémoire
	  for (int l = 1; l < nb_primes + 1; l++)
	    if (mpz_cmp_ui (smooth_list[l], 0) != 0 && a == 0)
	      gmp_printf ("[%d, %Zd]\n", prime_nbs[l - 1], smooth_list[l]);
	  for (int l = 0; l < cnt2; l++)
	    if (a == 0)
	      gmp_printf ("[%Zd, %Zd]\n", list_factor[l][0],
			  list_factor[l][1]);
	  for (int l = 0; l < 1000; l++)
	    {
	      for (int i = 0; i < 2; i++)
		mpz_clear (list_factor[l][i]);
	      free (list_factor[l]);
	    }
	  free (list_factor);
	}
    }
  free (prime_nbs);
  return smooth_list;
}


// Applique le crible de Dixon à n (imprime seulement la matrice et les relations, faire avec Sage ou PARI/GP ensuite)
void
dixon (mpz_t n, signed long int B)
{
  //Initialisation
  int nb_primes = 0;
  int cnt = 0;
  signed long int *prime_nbs = crible_erat (B);
  mpz_t rand;
  mpz_init (rand);
  mpz_t sq_rand;
  mpz_init (sq_rand);
  unsigned long int time_tmp;
  struct timeval t_tmp;
  gmp_randstate_t state;
  gmp_randinit_default (state);
  while (prime_nbs[nb_primes] <= B && prime_nbs[nb_primes] != 0)
    nb_primes++;

  unsigned long int **tab = malloc (nb_primes * sizeof (unsigned long int *));
  if (tab == NULL)
    {
      printf ("Dans dixon : malloc de tab impossible.\n");
      exit (EXIT_FAILURE);
    }

  for (int j = 0; j < nb_primes; j++)
    {
      tab[j] = malloc ((nb_primes + 10) * sizeof (unsigned long int));
      if (tab[j] == NULL)
	{
	  printf ("Dans dixon : malloc de tab[j] impossible.\n");
	  exit (EXIT_FAILURE);
	}
    }

  //Algorithme
  while (cnt < nb_primes + 10)
    {
      mpz_set_ui (sq_rand, 0);
      while (mpz_cmp_ui (sq_rand, 0) == 0 || mpz_cmp_ui (sq_rand, 1) == 0)
	{
	  gettimeofday (&t_tmp, NULL);
	  time_tmp = t_tmp.tv_sec * 1000000 + t_tmp.tv_usec;
	  gmp_randseed_ui (state, time_tmp);
	  mpz_urandomm (rand, state, n);
	  mpz_powm_ui (sq_rand, rand, 2, n);
	}

      mpz_t *smooth_list = friable (1, sq_rand, B);
      if (mpz_cmp_ui (smooth_list[0], 1) == 0)
	{
	  for (int j = 0; j < nb_primes; j++)
	    tab[j][cnt] = mpz_get_ui (smooth_list[j + 1]);
	  cnt++;
	  gmp_printf ("X_%d = %Zd\n", cnt, rand);
	}
      for (int j = 0; j < nb_primes + 1; j++)
	mpz_clear (smooth_list[j]);
      free (smooth_list);
    }

  //Impression de la matrice
  for (int k = 0; k < nb_primes; k++)
    {
      printf ("[");
      for (int l = 0; l < nb_primes + 9; l++)
	printf ("%lu, ", tab[k][l] % 2);
      printf ("%lu", tab[k][nb_primes + 9] % 2);
      printf ("],");
      printf ("\n");
    }

  //On libère la mémoire
  for (int k = 0; k < nb_primes; k++)
    free (tab[k]);
  free (tab);
  mpz_clear (rand);
  mpz_clear (sq_rand);
  gmp_randclear (state);
}

//Crible quadratique appliqué à n (imprime seulement la matrice, faire avec Sage ou PARI/GP ensuite)
void
crible_quadratique (mpz_t n, signed long int B)
{
  //Initialisation
  mpz_t root, m, a, tmp, tmp2, tmp3, two;
  mpz_init (root);
  mpz_sqrt (root, n);
  //mpz_add_ui (root, root, 1);
  mpz_init (m);
  mpz_pow_ui (m, root, 2);
  gmp_printf ("m est egal a %Zd\n", root);
  mpz_init (a);
  mpz_set_si (a, -150);
  mpz_init (tmp);
  mpz_init (tmp2);
  mpz_init (tmp3);
  mpz_init (two);
  mpz_set_ui (two, 2);
  mpz_t *list_tmp = friable (1, two, B);
  unsigned long int convert;
  signed long int *prime_nbs = crible_erat (B);
  int nb_prime = 0;
  while (prime_nbs[nb_prime] <= B && prime_nbs[nb_prime] != 0)
    nb_prime++;
  unsigned long int **matrix = calloc (301, sizeof (unsigned long int *));
  if (matrix == NULL)
    {
      printf ("Dans crible_quadratique : malloc de matrix impossible.\n");
      exit (EXIT_FAILURE);
    }
  for (int i = 0; i < 301; i++)
    {
      matrix[i] = calloc (nb_prime, sizeof (unsigned long int));
      if (matrix[i] == NULL)
	{
	  printf
	    ("Dans crible_quadratique : malloc de matrix[i] impossible.\n");
	  exit (EXIT_FAILURE);
	}
    }
  int a_tmp = 0;
  signed long int *list_a = calloc (301, sizeof (signed long int));
  if (list_a == NULL)
    {
      printf ("Dans crible_quadratique : malloc de list_a impossible.\n");
      exit (EXIT_FAILURE);
    }


  //Début de l'algorithme
  bool ligne = false;
  for (int i = 0; i < 301; i++)
    {
      //Calculs initiaux. A la fin, tmp corresponds au calcul de tmp = m^2 - n + a^2 + 2am mod n (voir rapport)
      mpz_pow_ui (tmp2, a, 2);
      mpz_mul (tmp3, a, root);
      mpz_mul (tmp3, tmp3, two);
      mpz_sub (tmp, m, n);
      mpz_add (tmp, tmp, tmp2);
      mpz_add (tmp, tmp, tmp3);
      if (mpz_cmp_ui (tmp, 0) < 0)
	{
	  mpz_mod (tmp, tmp, n);
	  mpz_sub (tmp, tmp, n);
	}
      else
	mpz_mod (tmp, tmp, n);

      //On regarde si tmp est friable
      list_tmp = friable (1, tmp, B);

      //Si tmp est friable on rentre les données qu'il faut dans la matrice
      if (mpz_cmp_ui (list_tmp[0], 1) == 0)
	{
	  for (int j = 1; j < nb_prime + 1; j++)
	    {
	      convert = mpz_get_ui (list_tmp[j]);
	      if (convert % 2 == 1)
		ligne = true;
	      matrix[i][j - 1] = convert % 2;
	    }
	  if (ligne == true)
	    {
	      list_a[a_tmp] = i - 100;
	      a_tmp++;
	      //On imprime les relations
	      gmp_printf ("%Zd =", a);
	      list_tmp = friable (0, tmp, B);
	    }
	  ligne = false;
	}
      mpz_add_ui (a, a, 1);
    }


  //On imprime seulement les lignes pour lesquels tmp est friable et lorsque au moins une valuation n'est pas pair.
  //Le modèle d'impression est utile pour Sage
  for (int i = 0; i < a_tmp; i++)
    {
      printf ("[");
      for (int j = 0; j < nb_prime - 1; j++)
	printf ("%lu,", matrix[list_a[i] + 100][j]);
      printf ("%lu", matrix[list_a[i] + 100][nb_prime - 1]);
      printf ("],");
      printf ("\n");
    }
  // On libère toute la mémoire utilisée
  for (int j = 0; j < nb_prime; j++)
    mpz_clear (list_tmp[j]);
  free (list_tmp);
  for (int i = 0; i < 301; i++)
    free (matrix[i]);
  free (matrix);
  free (prime_nbs);
  free (list_a);
  mpz_clear (root);
  mpz_clear (m);
  mpz_clear (a);
  mpz_clear (tmp);
  mpz_clear (tmp2);
  mpz_clear (tmp3);
  mpz_clear (two);
}

int
main (void)
{
  //Initialisation
  mpz_t n, tmp;
  mpz_init (n);
  mpz_init (tmp);
  signed long int algo, B, B_tmp;

  //Récupération des données
  gmp_printf ("Entrez un nombre Ã  factoriser: ", n);
  gmp_scanf ("%Zd", n);
  printf
    ("Quel algorithme de factorisation voulez vous utilisez?\n[1] Factorisation naive\n[2] P - 1 de Pollard\n[3] Crible de Dixon\n[4] Crible quadratique\nRentrez le nombre correspondant Ã  l'algorithme choisi.\n");
  gmp_scanf ("%Zd", tmp);
  algo = mpz_get_si (tmp);

  //Calcul de la borne B pour qu'elle soit optimale (deux méthodes différentes)
  if (algo == 4)
    {
      B_tmp = mpz_get_si (n);
      double log_tmp = log (log (B_tmp));
      log_tmp *= log(B_tmp);
      B = ceil (exp ( sqrt (log_tmp)/2 ));
      printf ("B = %ld\n", B);
    }
  else
    {
      B_tmp = mpz_get_si (n);
      B_tmp = abs (B_tmp);
      B = ceil (exp (sqrt (log (B_tmp))));
      printf ("B = %ld\n", B);
    }
  //Calcule du nombres de premiers avant B
  int nb_primes = 0;
  signed long int *prime_nbs = crible_erat (B);
  while (prime_nbs[nb_primes] <= B && prime_nbs[nb_primes] != 0)
      nb_primes++;
  free (prime_nbs);

  //En fonction de l'algorithme demandé, on fait la factorisation
  switch (algo)
    {
    case 1:
      printf ("Par la factorisation naive, nous obtenons:\n");
      mpz_t *list = friable (0, n, B);
      for (int k = 0; k < nb_primes + 1; k++)
	mpz_clear (list[k]);
      free (list);
      break;

    case 2:
      printf ("Par l'algorithme p-1 de Pollard, nous obtenons:\n");
      pollard (n, B);
      break;

    case 3:
      printf ("Par le crible de dixon, nous obtenons:\n");
      dixon (n, B);
      break;

    case 4:
      printf ("Par le crible quadratique, nous obtenons:\n");
      crible_quadratique (n, B);
      break;

    default:
      printf
	("Les algorithmes sont numerotes de 1 a 4. Choisissez l'un de ces nombres s'il vous plait.\n");
      exit (EXIT_FAILURE);
      break;
    }

  mpz_clear (n);
  mpz_clear (tmp);
  return EXIT_SUCCESS;
}
