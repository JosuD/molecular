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
	   double L, double kT, double thresh,
	   const gsl_rng *restrict rng, int flags) {
  
  int tries, n, k;
  unsigned long seed;
  double sigma;
  FILE *frand;
  register particle_t *swarm;
  double d_sq, x, y, z;
  gsl_rng *p_rng;
  gsl_complex c;
  
  /* IMPORTANTE: POR AHORA, EL RADIO DE SEITZ LO CALCULA PARA 2D,
   * PERO EL RESTO LO HACE PARA 3D. */

  /* Radio de Seitz:
   * 
   * 2D: rs^2 = L^2/(PI * N).
   * 3D: rs^3 = 3 * L^3/(4 * PI * N).
   *
   */

  const double rs_sq = L * L/(M_PI * N);
 
  sys->swarm = (particle_t *)calloc(N, sizeof(particle_t));
  swarm = sys->swarm;

  if (swarm != NULL) {

    /* Propiedades. */

    sys->N = N;
    sys->kT = kT;
    sys->u = 0;

    if (thresh == 0)
      thresh = rs_sq; /* Usa el radio de Seitz como thresh. */
    
    if ((flags & INIT_P_GEN) != 0) {
      
      /* Usa otro generador para las velocidades (a partir
	 del existente). Inicializa con una nueva semilla. */

      p_rng = gsl_rng_clone(rng);
      frand = fopen("/dev/urandom", "r");

      if (frand != NULL) {
	fread(&seed, sizeof(unsigned long), 1, frand);
	gsl_rng_set(p_rng, seed); /* Nueva semilla. */
	fclose(frand);
      }
      else
	gsl_rng_set(p_rng, (unsigned long)time(NULL));
    }
    else
      p_rng = rng;

    /* .:: Distribución espacial y de velocidades ::. */

    sigma = sqrt(PARTICLE_MASS * kT);

    if (thresh < 0) { /* Para agilizar... */

      for (n = 0; n < N; n++) {

	/* Espacial. */

	swarm[n].x = L * gsl_rng_uniform(rng) - L/2;
	swarm[n].y = L * gsl_rng_uniform(rng) - L/2;
	swarm[n].z = L * gsl_rng_uniform(rng) - L/2;
	
	/* De momentos. */

	swarm[n].px = gsl_ran_gaussian(p_rng, sigma);
	swarm[n].py = gsl_ran_gaussian(p_rng, sigma);
	swarm[n].pz = gsl_ran_gaussian(p_rng, sigma);

	/* Para la energía, es necesario hallar |p| y |r|, por lo que
	 * se calculan aún si INIT_POLAR no está en `flags'. */

	swarm[n].r = sqrt(swarm[n].x * swarm[n].x +
			  swarm[n].y * swarm[n].y +
			  swarm[n].z * swarm[n].z);

	/* En realidad, se calcula |p|^2 para la energía y 
	 * luego se toma la raíz. */
	swarm[n].p = swarm[n].px * swarm[n].px +
	  swarm[n].py * swarm[n].py +
	  swarm[n].pz * swarm[n].pz;
	
	/* .:: Energía ::. */

	/* Partícula. */
	swarm[n].K = swarm[n].p/(2 * PARTICLE_MASS);
	swarm[n].p = sqrt(swarm[n].p);

	/* Sistema. */
	sys->u += swarm[n].K;

	for (k = 0; k < n; ++k) { /* Potencial (con todos los anteriores). */

	  /* En este caso, no es `sq' :-D */
	  d_sq = sqrt((swarm[k].x - swarm[n].x) * (swarm[k].x - swarm[n].x) +
		      (swarm[k].y - swarm[n].y) * (swarm[k].y - swarm[n].y) +
		      (swarm[k].z - swarm[n].z) * (swarm[k].z - swarm[n].z));

	  sys->u += potencial(d_sq);
	}	  
      }
    }
    else { /* Hay que aplicar un umbral. */

      /* El procedimiento consiste en:
       *
       * 1. Generar un punto nuevo (x, y, z) ___bajo la suposición___ de que
       *    inicialmente, las posiciones de las partículas son estadísticamente
       *    independientes, por lo que se puede construir la posición a través de
       *    una muestra de tres números consecutivos del generador `gsl_rng_uniform'.
       *
       * 2. Comparar la distancia de este nuevo `centro' de partícula a las generadas
       *    anteriormente, con el thresh. Si la distancia es menor a este último 
       *    (esto es, si ambas partículas se encuentran en el mismo volumen), se
       *    descarta la terna generada `(x, y, z)' y se re-intenta con una nueva.
       *    El proceso se repite una cantidad de veces `INIT_MAX_TRIES', para evitar
       *    generar bucles infinitos. Si la densidad es muy alta, puede que se alcance
       *    rápido este límite. En cualquier caso, si se supera el límite de intentos,
       *    se guarda la posición generada aunque no cumpla con la condición de separación
       *    impuesta por el thresh. En todos los casos, el thresh es esférico, centrado
       *    en cada una de las ternas anteriormente generadas.
       *
       * 3. Se guarda la terna y se calculan las coordenadas polares (si así se desea).
       */

      tries = 0; /* Intentos. */

      swarm[0].x = L * gsl_rng_uniform(rng) - L/2;
      swarm[0].y = L * gsl_rng_uniform(rng) - L/2;
      swarm[0].z = L * gsl_rng_uniform(rng) - L/2;

      for (n = 1; n < N;) {

	x = L * gsl_rng_uniform(rng) - L/2;
	y = L * gsl_rng_uniform(rng) - L/2;
	z = L * gsl_rng_uniform(rng) - L/2;

	for (k = 0; k < n; k++) {
	
	  /* Distancia cuadrática entre ambos centros. */

	  d_sq = (swarm[k].x - x) * (swarm[k].x - x) +
	    (swarm[k].y - y) * (swarm[k].y - y);
	    //	    (swarm[k].z - z) * (swarm[k].z - z);
	
	  if (d_sq < thresh) /* muy cerca... */
	    break;
	}
	
	if (k < n && tries < INIT_MAX_TRIES) {
	  /* Intentos...tampoco esto debe convertirse en un bucle infinito. */
	  tries++; 
	  continue;
	}

	/* Almacena las ternas. */

	/* Posiciones. */
	swarm[n].x = x;
	swarm[n].y = y;
	swarm[n].z = z;

	/* Momentos. */
	swarm[n].px = gsl_ran_gaussian(p_rng, sigma);
	swarm[n].py = gsl_ran_gaussian(p_rng, sigma);
	swarm[n].pz = gsl_ran_gaussian(p_rng, sigma);

	/* .:: Energía ::. */

	/* Partícula. */
	swarm[n].r = sqrt(x * x + y * y + z * z);
	swarm[n].p = swarm[n].px * swarm[n].px +
	  swarm[n].py * swarm[n].py +
	  swarm[n].pz * swarm[n].pz;

	swarm[n].K = swarm[n].p/(2 * PARTICLE_MASS);
	swarm[n].p = sqrt(swarm[n].p);

	/* Sistema. */
	sys->u += swarm[n].K;
	
	/* Calculo la energía con todos los anteriores.*/

	for (k = 0; k < n; ++k) { /* Todas las anteriores. */

	  d_sq = sqrt((swarm[k].x - x) * (swarm[k].x - x) +
		      (swarm[k].y - y) * (swarm[k].y - y) +
		      (swarm[k].z - z) * (swarm[k].z - z));

	  sys->u += potencial(d_sq);
	}	  

	n++;
	tries = 0;
      }
    }
    
    /* ¿Hay que calcular polares? */

    if ((flags & INIT_POLAR) != 0) {

      for (n = 0; n < N; ++n) {

	/* Tal vez no haga falta usar la GSL para calcular los ángulos;
	 * en el sentido de que, quizá, la misma STD provea la misma
	 * precisión. De todas formas, ya que usamos la GSL para los rands... */

	/* Ángulos cenitales. */
	swarm[n].theta = GSL_REAL(gsl_complex_arccos_real(swarm[n].z/swarm[n].r));
	swarm[n].p_theta = GSL_REAL(gsl_complex_arccos_real(swarm[n].pz/swarm[n].p));

	/* Ángulos azimutales. */
	GSL_SET_COMPLEX(&c, swarm[n].x, swarm[n].y);
	swarm[n].phi = gsl_complex_arg(c);

	GSL_SET_COMPLEX(&c, swarm[n].px, swarm[n].py);
	swarm[n].p_phi = gsl_complex_arg(c);
      }
    }
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
