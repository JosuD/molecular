/*
 * Este es un script de prueba para las funciones de table.c/h
 */
#include <stdlib.h>
#include "table.h"
#include <stdio.h>

int main()
{
    float **tpot, **tforce;
    gen_tabla();
    int lineas;
    lineas = load_table(&tpot, &tforce);
    double potencial = apppot(tpot, 1.07);
    double fuerza = appforce(tforce, 1.07);

    printf("potencial %f, fuerza %f\n", potencial, fuerza);
}

