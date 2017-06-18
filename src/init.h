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

  double x, y, z,
    r, theta, phi,
    px, py, pz,
    p, p_theta, p_phi;
  
  double K; /* Energía cinética. */
  
};

typedef struct particle particle_t;

struct system {
  unsigned int N;
  double kT, u;
  particle_t *swarm;
};

typedef struct system system_t;

/* Opciones. */
#define INIT_SEITZ  ((int)1)
#define INIT_P_GEN  ((int)2)
#define INIT_POLAR  ((int)4)

/* Rutinas. */
void initialise(system_t * restrict, unsigned int, double, double,
		double, const gsl_rng * restrict, int);

ssize_t save_state(const char *, system_t * restrict, int);
ssize_t load_state(const char *, system_t * restrict, int);

#endif /* _INIT_H_ */
