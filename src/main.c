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


int main (int argc, char *argv[]) {

  int i, N;
  system_t sys;

  /* N = strtol(argv[1], NULL, 10); */

  N = 255;

  /* Argumentos: 
   * Puntero a system_t.
   * Número de partículas.
   * Densidad.
   * kT.
   * Distribución:
   *    SP_GRID: Grilla no necesariamente uniforme.
   *    SP_RAND: Aleatoria (en principio, uniforme)
                 en espacio reducido por 'Umbral/2'.
   * Umbral: usada sólo si SP_RAND. Define el umbral.
   *         si <= 0, se usa Seitz.
   */

  /* Devuelve el número de partículas agregadas. */
  i = initialise(&sys, N, 0.12, 1, SP_GRID, 0.8);
  printf("Agregadas: %3i/%3i.\n", i, N);


  sys_free(&sys);

  /* FILE *tabla; */
  /* tabla = fopen("tabla.dat", "w"); */

  /* for (i = 0; i < N; ++i) { */
  /*   fprintf(tabla, "%.4f %.4f %.4f\n", sys.swarm[i].x, */
  /* 	    sys.swarm[i].y, sys.swarm[i].z); */
  /* } */

  /* fflush(tabla); */
  /* fclose(tabla); */

  return 0;
}

