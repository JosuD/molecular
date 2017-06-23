/*

 */

#include "init.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h>

#include "gsl_math.h"
#include "gsl_complex.h"
#include "gsl_complex_math.h"

#include "table.h"

#define INIT_MAX_TRIES ((int)100)

void
initialise(system_t * restrict sys, unsigned int N,
	   double rho, double kT, double thresh, int flags) {
  
  int tries, n, k;
  unsigned long seed;
  double sigma;
  FILE *frand;
  register particle_t *swarm;
  double d_sq, x, y, z, dl, L;
  gsl_rng *rng;
  gsl_complex c;

  /* Radio de Seitz:
   * 2D: rs^2 = L^2/(PI * N).
   * 3D: rs^3 = 3 * L^3/(4 * PI * N).
   */

  const double rs = cbrt(3/(4 * M_PI * rho));
 
  sys->swarm = (particle_t *)calloc(N, sizeof(particle_t));
  rng = gsl_rng_alloc(gsl_rng_ranlxs2);

  swarm = sys->swarm;

  if (swarm != NULL && rng != NULL) {

    frand = fopen("/dev/urandom", "r");
      
    if (frand != NULL) {
      fread(&seed, sizeof(unsigned long), 1, frand);
      fclose(frand);
    }
    else
      seed = (unsigned long)time(NULL);

    gsl_rng_set(rng, seed);

    // La estimación es aproximadamente ~ L/cbrt(N) = 1/cbrt(rho).
    // dl es el diferencial de longitud

    /* Propiedades del sistema. */

    sys->N = N;
    sys->L = L;
    sys->kT = kT;
    sys->u = 0;

    n = 0;
    L = cbrt(N/rho);
    dl = 1/(cbrt(rho) + 2/L);

    sigma = sqrt(PARTICLE_MASS * kT);

    for (x = dl; x <= L - dl; x += dl) {

      for (y = dl; y <= L - dl ; y += dl) {

	for (z = dl; z <= L - dl; z += dl) {

	  // Ubica las partículas.
	  swarm[n].x = x;
	  swarm[n].y = y;
	  swarm[n].z = z;
	  
	  /* Momentos. */
	  swarm[n].px = gsl_ran_gaussian(rng, sigma);
	  swarm[n].py = gsl_ran_gaussian(rng, sigma);
	  swarm[n].pz = gsl_ran_gaussian(rng, sigma);

	  /* Energía. */
	  
	  /* Partícula. */
	  swarm[n].r = sqrt(x * x + y * y + z * z);
	  swarm[n].p = swarm[n].px * swarm[n].px +
	    swarm[n].py * swarm[n].py +
	    swarm[n].pz * swarm[n].pz;

	  swarm[n].K = swarm[n].p/(2 * PARTICLE_MASS);
	  swarm[n].p = sqrt(swarm[n].p);

	  /* Sistema. */
	  sys->u += swarm[n].K;

	  for (k = 0; k < n; ++k) { /* Todas las anteriores. */

	    d_sq = sqrt((swarm[k].x - x) * (swarm[k].x - x) +
		    (swarm[k].y - y) * (swarm[k].y - y) +
		    (swarm[k].z - z) * (swarm[k].z - z));

	    sys->u += potencial(d_sq);
	  }	  

	  n++;
	  if (n > 511)
	    goto _a;

	}
      }
    }

  _a:    
    printf("Delta l:%i/%i.\n", n, N);

    gsl_rng_free(rng);
  }
}

ssize_t
save_state (const char *path, system_t * restrict state, int table) {

  int fbin, n;
  ssize_t cwrote, wrote;
  size_t total;
  int8_t *buff;
  FILE *ftable;
  char *tpath;

  if (table != 0) {  /* Si se desea, guarda el estado actual como tabla. */

    wrote = strlen(path) + 5; /* Longitud de la extensión + NULL. */
    tpath = (char *)malloc(wrote);
    snprintf(tpath, wrote, "%s.txt", path);

    ftable = fopen(tpath, "w");
    if (ftable == NULL)
      ftable = stdout;

    fprintf(ftable, "# particula K x y z r theta phi px py pz p p_theta p_phi\n");

    for (n = 0; n < state->N; ++n)
      fprintf(ftable, "%3i %.5f %.5f %.5f %.5f %.5f %.5f "
	      "%.5f %.5f %.5f %.5f %.5f %.5f %.5f\n",
	      n, state->swarm[n].K,
	      state->swarm[n].x, state->swarm[n].y, state->swarm[n].z,
	      state->swarm[n].r, state->swarm[n].theta, state->swarm[n].phi,
	      state->swarm[n].px, state->swarm[n].py, state->swarm[n].pz,
	      state->swarm[n].p, state->swarm[n].p_theta, state->swarm[n].p_phi);
      
    fflush(ftable);
    if (ftable != stdout)
      fclose(ftable);

    snprintf(tpath, wrote, "%s.bin", path);
  }
  else
    tpath = path;

  /* Almacena el bloque actual que representa el estado del sistema. */

  wrote = -1;
  fbin = open(tpath, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);

  if (fbin > 0) {

    /* Formato en bloque: guarda las propiedades y luego el enjambre.
     * Nota: guarda también la dirección del puntero `*swarm'. */

    write(fbin, state, sizeof(system_t)); /* Propiedades. */
    fsync(fbin);

    wrote = 0;
    buff = (int8_t *)state->swarm;
    total = state->N * sizeof(particle_t);

    do {

      cwrote = write(fbin, buff + wrote , total - wrote);

      if (cwrote < 0) { /* En caso de error, deja de escribir y cierra. */
	close(fbin);
	break;
      }
      
      wrote += cwrote;
    } while (wrote < total);

    /* Cierra (si no se hizo antes). */
    if (cwrote >= 0)
      close(fbin);
  }

  if (tpath != path)
    free(tpath);
  
  return (wrote + sizeof(system_t));
}

ssize_t
load_state(const char *path, system_t * restrict state) {

  int fdes;
  off_t flen;
  size_t sw_size;
  int8_t *buff;
  ssize_t sread, cread;

  /* Lee el estado almacenado en `path' directamente sobre `state'.
   * Es importante preservar el orden de la estructura `system', por
   * condiciones de `padding' y alineamiento (al parecer). */

  sread = 0;
  fdes = open(path, O_RDONLY);
  
  if (fdes > 0) {

    flen = lseek(fdes, 0, SEEK_END);

    if (flen > sizeof(system_t)) {

      lseek(fdes, 0, SEEK_SET);
      buff = (int8_t *)state;
      sread = 0;

      do {

	cread = read(fdes, buff + sread, sizeof(system_t) - sread);

	if (cread < 0) {
	  close(fdes);
	  fdes = -1;
	  break;
	}

	sread += cread;
      } while (sread < sizeof(system_t));

      if (fdes > 0) {

	/* Lee el enjambre, tomando el tamaño indicado en las propiedades. */

	state->swarm = (particle_t *)calloc(state->N, sizeof(particle_t));
	sw_size = state->N * sizeof(particle_t);

	if (state->swarm != NULL) {

	  sread = 0;
	  buff = (int8_t *)state->swarm;

	  do {
	    cread = read(fdes, buff + sread, sw_size - sread);

	    if (cread < 0) {
	      close(fdes);
	      break;
	    }

	    sread += cread;
	  } while (sread < sw_size);
	}

	sread += sizeof(system_t);

	/* Cierra (si no se hizo antes). */
	if (cread >= 0)
	  close(fdes);
      }
      else
	sread = -sread;
    }
  }

  return sread;
}

void
sys_free(system_t * restrict system) {

  free(system->swarm);
  memset(system, 0, sizeof(system_t));
}
