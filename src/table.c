#include <stdio.h>
#include <math.h>
#include "table.h"
#include <stdlib.h>

#define RMAX 10
#define SIGMA 2.5
#define DELTA 0.01
#define ROUT 3


int gen_tabla(){
    double r;
    int i;
    FILE * tabla_fuerza;
    FILE * tabla_potencial;
    tabla_fuerza = fopen("tabla_fuerza.dat", "w+");
    tabla_potencial = fopen("tabla_potencial.dat", "w+");
    for(i=0;i<RMAX/DELTA; i++){
        r = DELTA + i*DELTA;
        fprintf(tabla_fuerza, "%f %f\n", r, force(r));
        fprintf(tabla_potencial, "%f %f\n", r, potencial(r));
    }
    fclose(tabla_fuerza);
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

double potencial(double r){
    double v;
    if(r<SIGMA){
        v =-4*(1/pow(r,12)+1/pow(r,6));
    } else if (r>SIGMA && r<ROUT){
        v = 4*(1/pow(r,12)+1/pow(r,6));
        v *= pow((pow(ROUT,2)-pow(r,2)),2);
        v *= (pow(r,2)+2*pow(r,2)-3*pow(SIGMA,2));
        v *= 1/(pow((pow(ROUT,2)-pow(SIGMA,2)),3));
    } else {
        v = 0;
    }
    return v;
}

int load_table(float **tpot, float **tforce){
    FILE * tabla_potencial = fopen("tabla_potencial.dat","r");
    FILE * tabla_fuerza = fopen("tabla_fuerza.dat","r");

    int i = 0;
    int lineas = 0;
    float temp1, temp2;
    while(EOF != fscanf(tabla_potencial, " %f %f",  &temp1, &temp2))
    {
        printf("prmero %f segundo %f\n", temp1, temp2);
        lineas++;
    }
    tpot = malloc(sizeof(float)*lineas);
    if (tpot)
    {
        for (i = 0; i < lineas; i++)
        {
            tpot[i] = malloc(sizeof(float)*2);
        }
    }
    rewind(tabla_potencial);
    i = 0;
    while( EOF != fscanf(tabla_potencial, " %f %f", &tpot[i][0], &tpot[i][1])){
        i++;
    }

    tpot = malloc(sizeof(float)*lineas);
    if (tforce)
    {
        for (i = 0; i < lineas; i++)
        {
            tforce[i] = malloc(sizeof(float)*2);
        }
    }
    i = 0;
    while( EOF != fscanf(tabla_fuerza, " %f %f", &tforce[i][0], &tforce[i][1])){i++;}
    

    return lineas;
}

//double apppot(double *tpot, double r){
//    return potencial;
//}
//
//double appforce(double *tforce, double r){
//    return fuerza;
//}
