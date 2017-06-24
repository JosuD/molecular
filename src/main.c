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

int main (int argc, char *argv[]) {

  int i, j, N;
  system_t sys;

  /* N = strtol(argv[1], NULL, 10); */

  N = 512;

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
  i = initialise(&sys, N, 0.8442, 1, SP_GRID, 0);

  /* Descomentar si se quiere guardar una tabla. */

  /* FILE *tabla; */
  /* tabla = fopen("tabla.dat", "w"); */

  /* for (i = 0; i < N; ++i) { */
  /*   fprintf(tabla, "%.4f %.4f %.4f\n", sys.swarm[i].x, */
  /* 	    sys.swarm[i].y, sys.swarm[i].z); */

  /*     printf("%.3f\n", sys.swarm[i].px * sys.swarm[i].px + */
  /* 	     sys.swarm[i].py * sys.swarm[i].py + */
  /* 	     sys.swarm[i].pz * sys.swarm[i].pz); */

  /* } */

  /* fflush(tabla); */
  /* fclose(tabla); */

  printf("Agregadas: %3i/%3i.\n", i, N);

  float **tforce, **tpot;
  load_table(&tpot, &tforce);

  for (i = 0; i < 100; i++) {
    verlet(&sys, N, sys.L, tforce, 1e-4);

    for (j = 0; j < N; j++) /* Imprime |p| para cada partícula */
      printf("%.3f\n", sys.swarm[j].px * sys.swarm[j].px +
  	     sys.swarm[j].py * sys.swarm[j].py +
  	     sys.swarm[j].pz * sys.swarm[j].pz);

    getchar(); // Pausa.

  }

  sys_free(&sys);

  return 0;
}

