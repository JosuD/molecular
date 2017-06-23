/*
 */

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <error.h>
#include <math.h>

#include "init.h"
#include "table.h"
#include "verlet.h"

/* Algoritmo de Verlet (`Velocity Verlet'):

   r(t + dt) = r(t) + v(t) * dt + (dt^2/2m) * f(t)
   v(t + dt) = v(t) + dt/2m * [ f(t) + f(t + dt) ]
*/

#define N ((int)250)

int main (int argc, char *argv[]) {

  int i;
  system_t sys;

  initialise(&sys, 512, 0.8442, 1, 0, 0);

  for (i = 0; i < 512; i++)
    printf("%3i. r = %.5f\n", i, sys.swarm[i].r);
  
  /* float **tpot, **tforce; */
  /* load_table(&tpot, &tforce); */
  /* verlet(&sys, N, 40, tforce, 0.1);  */
  sys_free(&sys);

  return 0;
}

