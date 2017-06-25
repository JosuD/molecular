#include "init.h"

void verlet(system_t * sys, unsigned int N, double L, float **tforce, double dt);
void fuerza_vec(particle_t *swarmc, double *force_x, double *force_y, double *force_z, int N, int n, float **tforce, double L);
void vecf(particle_t *swarmc,double *force_x, double *force_y, double *force_z, int n, int m, int a, int b, int c, double L, float **tforce);

