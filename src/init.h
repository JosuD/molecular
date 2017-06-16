/*
 */

#ifndef _INICIAL_H_
#define _INICIAL_H_

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

enum errors = {

  ERROR_RNG_NOT_INITIATED

};

typedef enum errors error_t;

struct particula {

  int id;

  double x, y, z,
    r, theta, phi,
    px, py, pz,
    p, p_theta, p_phi;

};

typedef struct particula particula_t;

#endif /* _INICIAL_H_ */
