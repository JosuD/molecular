/*
 */

#ifndef _INIT_H_
#define _INIT_H_

#include <time.h>
#include <unistd.h>

#include "gsl_rng.h"
#include "gsl_randist.h"

#define PARTICLE_MASS ((double)1)

struct particle {

  double K; /* Energía cinética. */

  double x, y, z,
    r, theta, phi,
    px, py, pz,
    p, p_theta, p_phi;
};

typedef struct particle particle_t;

struct system {
  unsigned int N;
  double kT, u;
  particle_t *swarm;
};

typedef struct system system_t;

/* Opciones. */
#define INIT_P_GEN  ((int)1)
#define INIT_POLAR  ((int)2)

/* Rutinas. */

/* initialise -- inicia el sistema pasado en system_t:
 *
 * Argumentos:
 *   system_t *  : puntero a una estructura de tipo system_t.
 *   unsigned int: cantidad de partículas.
 *   double      : lado del cubo de la caja contenedora (en unidades de sigma).
 *   double      : kT -- Temperatura.
 *   double      : umbral. 
 *               : Si es < 0: no se usa un umbral.
 *               : si es = 0: se usa el Radio de Seitz.
 *               : si es > 0: se usa el valor especificado.
 *   gsl_rng *   : estado del generador de números aleatorios.
 *   int         : opciones, definidas más arriba con el prefijo INIT_
 *                 INIT_P_GEN: usa otro generador para los momentos.
 *                 INIT_POLAR: calcula los ángulos en coordenadas polares.
 */

void initialise(system_t * restrict, unsigned int, double, double,
		double, const gsl_rng * restrict, int);

/* save_state -- guarda el estado del sistema en system_t en un archivo binario.
 *
 * Argumentos:
 *   const char * : ruta al archivo donde guardar los datos (tomado como 'tallo').
 *   system_t   * : puntero al estado del sistema que se busca guardar.  
 *   int          : booleano. Si 1, guarda una tabla con los datos de cada partícula.
 *                  Agrega en este caso la extensión '.txt' al nombre de archivo y
 *                  '.bin' al nombre binario.
 *
 * Devuelve: número de bytes escritos en el archivo binario.
 */

ssize_t save_state(const char *, system_t * restrict, int);

/* load_state -- carga el estado del sistema en el archivo binario especificado.
 *
 * Argumentos:
 *   const char * : ruta del archivo.
 *   system_t   * : puntero a la estructura donde almacenar los bytes leídos.
 *
 * Devuelve: número de bytes leídos.
 */

ssize_t load_state(const char *, system_t * restrict);

/* sys_free -- libera recursos usados por un sistema.
 *
 * Argumentos:
 *   system_t * : puntero a la estructura que contiene recursos a liberar.
 */

void sys_free(system_t * restrict);

#endif /* _INIT_H_ */
