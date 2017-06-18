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
};

typedef struct particle particle_t;
typedef particle_t * system_t;

/* Opciones. */
#define INIT_SEITZ  ((int)1)
#define INIT_P_GEN  ((int)2)
#define INIT_POLAR  ((int)4)

/* Rutinas. */
void initialise(system_t * restrict, int, double, double,
		double, const gsl_rng * restrict, int);

ssize_t save_state(const char *, system_t, size_t, int);
ssize_t load_state(const char *, system_t * restrict, size_t);

#endif /* _INIT_H_ */
