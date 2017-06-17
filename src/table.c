#include <stdio.h>
#include <math.h>
#include "table.h"
#include <stdlib.h>

#define RMAX 10
#define SIGMA 2.5
#define DELTA 0.01


int main(){
    double r;
    int i;
    FILE * tabla;
    tabla = fopen("tabla.dat", "w+");
    for(i=0;i<RMAX/DELTA; i++){
        r = i*DELTA;
        fprintf(tabla, "%f %f\n", r, force(r));
    }
    fclose(tabla);
    return 0;
}

double force(double r){
    double fuerza;
    fuerza = 4*(-12/pow(r,13)+6/pow(r,7));
    return fuerza;
}
