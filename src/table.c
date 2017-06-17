#include <stdio.h>
#include <math.h>
#include "table.h"
#include <stdlib.h>

#define RMAX 10
#define SIGMA 2.5
#define DELTA 0.01
#define ROUT 3


int main(){
    double r;
    int i;
    FILE * tabla;
    tabla = fopen("tabla.dat", "w+");
    for(i=0;i<RMAX/DELTA; i++){
        r = DELTA + i*DELTA;
        fprintf(tabla, "%f %f\n", r, force(r));
    }
    fclose(tabla);
    return 0;
}

double force(double r){
    double fuerza;
    double st;
    if(r<SIGMA){
        fuerza = -4*(-12/pow(r,13)+6/pow(r,7));
    } else if (r>SIGMA && r<ROUT){
        fuerza = -4*(-12/pow(r,13)+6/pow(r,7));
        fuerza *= pow((pow(ROUT,2)-pow(r,2)),2);
        fuerza *= (pow(r,2)+2*pow(r,2)-3*pow(SIGMA,2));
        fuerza *= 1/(pow((pow(ROUT,2)-pow(SIGMA,2)),3));
        st = -4*(1/pow(r,12)-1/pow(r,6));
        st *= (pow(ROUT,2)-pow(r,2))/pow(pow(ROUT,2)-pow(SIGMA,2),3);
        st *= 4*r*(2*pow(ROUT,2)+pow(r,2)-3*pow(SIGMA,2));
        fuerza = fuerza + st;
    } else {
        fuerza = 0;
    }
    return fuerza;
}
