/*
 */

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <error.h>
#include <math.h>

#include "init.h"

/* Algoritmo de Verlet (`Velocity Verlet'):

   r(t + dt) = r(t) + v(t) * dt + (dt^2/2m) * f(t)
   v(t + dt) = v(t) + dt/2m * [ f(t) + f(t + dt) ]

*/


#define N ((int)250)

int main (int argc, char *argv[]) {

  int status, n;
  unsigned long s;
  FILE *frand;
  gsl_rng *rng;
  system_t sys;
  
  frand = fopen("/dev/urandom", "r");

  if (frand == NULL)
    s = (unsigned long)time(NULL);
  else {
    fread(&s, sizeof(unsigned long), 1, frand);
    fclose(frand);
  }

  /* Pruebas. */

  rng = gsl_rng_alloc(gsl_rng_ranlxs2);

  if (rng != NULL) {

    gsl_rng_set(rng, s);

    initialise(&sys, N, 40, 10, 0, rng, INIT_SEITZ);
    printf("Escrito: %lu\n", save_state("prueba", sys, N, 1));
    free(sys);

    sys = NULL;
    printf("Leído: %lu\n", load_state("prueba.bin", &sys, 0));

    for (n = 0; n < N; ++n)
      printf("%3i %.5f %.5f %.5f %.5f %.5f %.5f "
	     "%.5f %.5f %.5f %.5f %.5f %.5f\n",
	     n, sys[n].x, sys[n].y, sys[n].z,
	     sys[n].r, sys[n].theta, sys[n].phi,
	     sys[n].px, sys[n].py, sys[n].pz,
	     sys[n].p, sys[n].p_theta, sys[n].p_phi);

    free(sys);
    gsl_rng_free(rng);
  }

  return 0;
}

