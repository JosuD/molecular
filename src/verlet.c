#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <error.h>
#include <math.h>

#include "verlet.h"
#include "init.h"
#include "table.h"

#define ROUT 2.95

void verlet(system_t * sys, unsigned int N, double L, float **tforce, double dt){
    int n;
    particle_t *swarmc;
    particle_t *swarm;
    swarm = sys->swarm;
    swarmc = (particle_t *)calloc(N, sizeof(particle_t));
    double force_x, force_y, force_z;
    double force_xdt, force_ydt, force_zdt;

    for(n = 0; n < N; n++){
        swarmc[n].x = swarm[n].x;
        swarmc[n].y = swarm[n].y;
        swarmc[n].z = swarm[n].z;
        swarmc[n].px = swarm[n].px;
        swarmc[n].py = swarm[n].py;
        swarmc[n].pz = swarm[n].pz;
        swarmc[n].r = swarm[n].r;
        swarmc[n].theta = swarm[n].theta;
        swarmc[n].phi = swarm[n].phi;
        swarmc[n].p = swarm[n].p;
        swarmc[n].p_theta = swarm[n].p_theta;
        swarmc[n].p_phi = swarm[n].p_phi;
    }
    for(n = 0; n < N; n++){
        fuerza_vec(swarmc, &force_x, &force_y, &force_z, N, n, tforce);

        swarm[n].x = swarm[n].x + swarm[n].px*dt+force_x*dt*dt/2;
        swarm[n].y = swarm[n].y + swarm[n].py*dt+force_y*dt*dt/2;
        swarm[n].z = swarm[n].z + swarm[n].pz*dt+force_z*dt*dt/2;

        fuerza_vec(swarm, &force_xdt, &force_ydt, &force_zdt, N, n, tforce);

        swarm[n].px = swarm[n].px + (force_x+force_xdt)*dt/2;
        swarm[n].py = swarm[n].py + (force_y+force_ydt)*dt/2;
        swarm[n].pz = swarm[n].pz + (force_z+force_zdt)*dt/2;

        swarm[n].r = sqrt(swarm[n].x * swarm[n].x +
                  swarm[n].y * swarm[n].y +
                  swarm[n].z * swarm[n].z);

        swarm[n].p = swarm[n].px * swarm[n].px +
            swarm[n].py * swarm[n].py +
          swarm[n].pz * swarm[n].pz;
        
        swarm[n].K = swarm[n].p/(2 * PARTICLE_MASS);
        swarm[n].p = sqrt(swarm[n].p);


    }
}


void fuerza_vec(particle_t *swarmc,
        double *force_x, double *force_y, double *force_z,
        int N, int n, float **tforce)
{
    double force = 0;
    *force_x = 0, *force_y = 0, *force_z = 0;
    double distancia;
    int m;
    
    double dir_x, dir_y, dir_z;

    for( m = 0; m < N; m++){
        if ( m == n) continue;
        distancia = sqrt(pow(swarmc[m].x-swarmc[n].x,2)+pow(swarmc[m].x-swarmc[n].x,2)+pow(swarmc[m].x-swarmc[n].x,2));
        if(distancia>ROUT) continue;
        dir_x = (swarmc[m].x-swarmc[n].x)/distancia;
        dir_y = (swarmc[m].y-swarmc[n].y)/distancia;
        dir_z = (swarmc[m].y-swarmc[n].z)/distancia;
       force = appforce(tforce, distancia);
       *force_x += dir_x * force;
       *force_y += dir_y * force;
       *force_z += dir_z * force;
    }


}
