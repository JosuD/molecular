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

#define INIT_MAX_TRIES ((int)100)

void
initialise(system_t * restrict sys_p, int N,
	   double L, double kT, double thresh,
	   const gsl_rng *restrict rng, int flags) {
  
  int tries, n, k;
  unsigned long seed;
  double sigma;
  FILE *frand;
  system_t sys;
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
  
  *sys_p = (particle_t *)calloc(N, sizeof(particle_t));
  sys = *sys_p;

  if (sys != NULL) {

    if ((flags & INIT_SEITZ) != 0)
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

    if (thresh <= 0) { /* Para agilizar... */

      for (n = 0; n < N; n++) {

	/* Espacial. */

	sys[n].x = L * gsl_rng_uniform(rng) - L/2;
	sys[n].y = L * gsl_rng_uniform(rng) - L/2;
	sys[n].z = L * gsl_rng_uniform(rng) - L/2;
	
	/* De momentos. */

	sys[n].px = gsl_ran_gaussian(p_rng, sigma);
	sys[n].py = gsl_ran_gaussian(p_rng, sigma);
	sys[n].pz = gsl_ran_gaussian(p_rng, sigma);
      }
    }
    else { /* Hay que aplicar un umbral. */

      /* El procedimiento consyste en:
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

      sys[0].x = L * gsl_rng_uniform(rng) - L/2;
      sys[0].y = L * gsl_rng_uniform(rng) - L/2;
      sys[0].z = L * gsl_rng_uniform(rng) - L/2;

      for (n = 1; n < N;) {

	x = L * gsl_rng_uniform(rng) - L/2;
	y = L * gsl_rng_uniform(rng) - L/2;
	z = L * gsl_rng_uniform(rng) - L/2;

	for (k = 0; k < n; k++) {
	
	  /* Distancia cuadrática entre ambos centros. */

	  d_sq = (sys[k].x - x) * (sys[k].x - x) +
	    (sys[k].y - y) * (sys[k].y - y);
	    //	    (sys[k].z - z) * (sys[k].z - z);
	
	  if (d_sq < thresh) /* muy cerca... */
	    break;
	}
	
	if (k < n && tries < INIT_MAX_TRIES) {
	  /* Intentos...tampoco esto debe convertirse en un bucle infinito. */
	  tries++; 
	  continue;
	}

	/* Guarda la terna generada. */

	/* Posiciones. */
	sys[n].x = x;
	sys[n].y = y;
	sys[n].z = z;

	/* Momentos. */
	sys[n].px = gsl_ran_gaussian(p_rng, sigma);
	sys[n].py = gsl_ran_gaussian(p_rng, sigma);
	sys[n].pz = gsl_ran_gaussian(p_rng, sigma);

	n++;
	tries = 0;
      }
    }

    /* ¿Hay que calcular polares? 
       Nota: Ver si es mejor poner esto dentro de cada bucle. (¿Muchos ifs?). */

    if ((flags & INIT_POLAR) != 0) {

      for (n = 0; n < N; ++n) {

	/* Tal vez no haga falta usar la GSL para calcular los ángulos;
	 * en el sentido de que, quizá, la misma STD provea la misma
	 * precisión. De todas formas, ya que usamos la GSL para los rands... */

	/* Posiciones. */

	sys[n].r = sqrt(sys[n].x * sys[n].x +
			sys[n].y * sys[n].y +
			sys[n].z * sys[n].z);

	/* Ángulo cenital. */
	sys[n].theta = GSL_REAL(gsl_complex_arccos_real(sys[n].z/sys[n].r));

	/* Ángulo azimutal. */
	GSL_SET_COMPLEX(&c, sys[n].x, sys[n].y);
	sys[n].phi = gsl_complex_arg(c);

	/* Momentos. */
	sys[n].p = sqrt(sys[n].px * sys[n].px +
			sys[n].py * sys[n].py +
			sys[n].pz * sys[n].pz);

	sys[n].p_theta = GSL_REAL(gsl_complex_arccos_real(sys[n].pz/sys[n].p));

	GSL_SET_COMPLEX(&c, sys[n].px, sys[n].py);
	sys[n].p_phi = gsl_complex_arg(c);
      }
    }
  }
}

ssize_t
save_state (const char *path, system_t state, size_t size_n, int table) {

  int fbin, n;
  ssize_t cwrote, wrote;
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

    for (n = 0; n < size_n; ++n)
      fprintf(ftable, "%3i %.5f %.5f %.5f %.5f %.5f %.5f "
	      "%.5f %.5f %.5f %.5f %.5f %.5f\n",
	      n, state[n].x, state[n].y, state[n].z,
	      state[n].r, state[n].theta, state[n].phi,
	      state[n].px, state[n].py, state[n].pz,
	      state[n].p, state[n].p_theta, state[n].p_phi);
      
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

    buff = (int8_t *)state;
    size_n *= sizeof(particle_t);
    wrote = 0;

    do {

      cwrote = write(fbin, buff + wrote , size_n - wrote);

      if (cwrote < 0) { /* En caso de error, deja de escribir y cierra. */
	close(fbin);
	break;
      }
      
      wrote += cwrote;
    } while (wrote < size_n);

    /* Cierra (si no se hizo antes). */
    if (cwrote >= 0)
      close(fbin);
  }

  if (tpath != path)
    free(tpath);
  
  return wrote;
}

ssize_t
load_state(const char *path, system_t * restrict state, size_t size_n) {

  int fdes;
  off_t flen;
  int8_t *buff;
  ssize_t sread, cread;

  /* Lee el estado almacenado en `path' directamente sobre `state'. */

  sread = -1;
  fdes = open(path, O_RDONLY);
  
  if (fdes > 0) {

    /* Busca el tamaño del archivo usando `lseek'; si hay un error,
     * calcula la cantidad de bytes a leer usando el parámetro size_n.
     *
     * Importante:
     *
     *   1. Si `state' es un bloque ya alojado debe ser suficientemente
     *      grande como para alojar la cantidad de bytes en el archivo o
     *      (si no se pudiera hallar esta cantidad con `lseek') calculados
     *      con `size_n'.
     *
     *   2. Si existe alguna diferencia en el estado almacenado y el estado
     *      especificado mediante `size_n', o si se leen menos bytes desde
     *      el archivo, es posible que la estructura en state[0] esté
     *      __a medio escribir__. Tener cuidado con esto.
     */

    flen = lseek(fdes, 0, SEEK_END);

    if (flen < 0)
      size_n *= sizeof(particle_t);
    else
      size_n = flen;

    lseek(fdes, 0, SEEK_SET);

    if (*state == NULL) /* Hay que alojar un bloque nuevo. */
      *state = (system_t)malloc(size_n);

    if (*state != NULL) {

      /* Lee bloques de bytes del archivo y va copiándolos a memoria. */
      buff = (int8_t *)(*state);
      sread = 0;

      do {

	cread = read(fdes, buff + sread, size_n - sread);

	/* En caso de error, deja de leer inmediatamente y cierra
	 * el archivo. La cantidad de bytes devueltos (por la rutina)
	 * indica el total de bytes leídos. */

	if (cread < 0) {
	  close(fdes);
	  break;
	}
	
	sread += cread;

      } while (cread > 0 && sread < size_n);

      /* Cierra el archivo (si no se hizo antes). */
      if (cread >= 0)
	close(fdes);
    }
  }

  return sread;
}
