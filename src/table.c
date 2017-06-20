#include <stdio.h>
#include <math.h>
#include "table.h"
#include <stdlib.h>

#define RMAX 3
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
    fclose(tabla_potencial);
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

int load_table(float ***tpot, float ***tforce){
    FILE * tabla_potencial = fopen("tabla_potencial.dat","r");
    FILE * tabla_fuerza = fopen("tabla_fuerza.dat","r");

    int i = 0;
    int lineas = 0;
    float temp1, temp2;
    while(EOF != fscanf(tabla_potencial, " %f %f",  &temp1, &temp2))
    {
        lineas++;
    }

    *tpot = malloc(sizeof(float *)*lineas);
    for (i = 0; i < lineas; i++)
    {
        (*tpot)[i] = malloc(sizeof(float)*2);
    }

    rewind(tabla_potencial);
    i = 0;
    while( EOF != fscanf(tabla_potencial, " %f %f", &(*tpot)[i][0], &(*tpot)[i][1])){
        i++;
    }

    *tforce = malloc(sizeof(float *)*lineas);
    for (i = 0; i < lineas; i++)
    {
        (*tforce)[i] = malloc(sizeof(float)*2);
    }

    i = 0;
    while( EOF != fscanf(tabla_fuerza, " %f %f", &(*tforce)[i][0], &(*tforce)[i][1])){i++;}
    

    return lineas;
}

double apppot(float **tpot, double r){

    if(r>ROUT) return (double) 0.0;
    float h = tpot[1][0] - tpot[0][0];
    int candidato = (int)floor(r/h);
    float potencial;
    if(tpot[candidato][0] == (float)r){
        return (double) tpot[candidato][1];
    }else{
        float pendiente = (tpot[candidato+1][1]-tpot[candidato][1])/(tpot[candidato+1][0]-tpot[candidato][0]);
        potencial = tpot[candidato][1] + pendiente * (r - tpot[candidato][0]);

    }

    return (double) potencial;
}

double appforce(float **tforce, double r){

    if(r>ROUT) return (double) 0.0;
    float h = tforce[1][0] - tforce[0][0];
    int candidato = (int)floor(r/h);
    float fuerza;
    if(tforce[candidato][0] == (float)r){
        return (double) tforce[candidato][1];
    }else{
        float pendiente = (tforce[candidato+1][1]-tforce[candidato][1])/(tforce[candidato+1][0]-tforce[candidato][0]);
        fuerza = tforce[candidato][1] + pendiente * (r - tforce[candidato][0]);
    }

    return (double) fuerza;
}
