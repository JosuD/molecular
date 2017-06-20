/*
 * Este es un script de prueba para las funciones de table.c/h
 */
#include <stdlib.h>
#include "table.h"
#include <stdio.h>

int main()
{
    float **tpot, **tforce;
    int lineas;
    lineas = load_table(&tpot, &tforce);
    printf("%f %f\n", tpot[0][0], tpot[1][0]);
    double potencial = apppot(tpot, 1.07);
    double fuerza = appforce(tforce, 1.07);

    printf("potencial %f, fuerza %f\n", potencial, fuerza);
}

