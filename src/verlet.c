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

// swarmc será el estado en un tiempo t
// primero copio el estado actual
    for(n = 0; n < N; n++){
        swarmc[n].x = swarm[n].x;
        swarmc[n].y = swarm[n].y;
        swarmc[n].z = swarm[n].z;
        swarmc[n].px = swarm[n].px;
        swarmc[n].py = swarm[n].py;
        swarmc[n].pz = swarm[n].pz;
        swarmc[n].p = 0;
    }
    for(n = 0; n < N; n++){
        //para el primer paso necesito la fuerza en tiempo t
        //fuerza_vec devuelve las componente sde la fuerza
        fuerza_vec(swarmc, &force_x, &force_y, &force_z, N, n, tforce, L);

        //con la fuerza calculo la nueva posicion y lo almaceno en swarm
        //ahi pongo la posicion en tiempo t+dt
        swarm[n].x = swarm[n].x + swarm[n].px*dt+force_x*dt*dt/2;
        swarm[n].y = swarm[n].y + swarm[n].py*dt+force_y*dt*dt/2;
        swarm[n].z = swarm[n].z + swarm[n].pz*dt+force_z*dt*dt/2;
//        printf("fuerza %f\n", force_x);
//condiciones periodicas de contorno
        swarm[n].x = fabs(swarm[n].x) - L*floor(fabs(swarm[n].x)/L);
        swarm[n].y = fabs(swarm[n].y) - L*floor(fabs(swarm[n].y)/L);
        swarm[n].z = fabs(swarm[n].z) - L*floor(fabs(swarm[n].z)/L);

        //necesito la fuerza en t+dt para actualizar las velocidades
        fuerza_vec(swarm, &force_xdt, &force_ydt, &force_zdt, N, n, tforce, L );

        swarm[n].px = swarm[n].px + (force_x+force_xdt)*dt/2;
//        printf("force_x %f\n", force_x);
        swarm[n].py = swarm[n].py + (force_y+force_ydt)*dt/2;
        swarm[n].pz = swarm[n].pz + (force_z+force_zdt)*dt/2;

 //       swarm[n].r = sqrt(swarm[n].x * swarm[n].x +
//                  swarm[n].y * swarm[n].y +
//                  swarm[n].z * swarm[n].z);

        swarm[n].p = swarm[n].px * swarm[n].px +
            swarm[n].py * swarm[n].py +
          swarm[n].pz * swarm[n].pz;
        
        swarm[n].K = swarm[n].p/(2 * PARTICLE_MASS);
        swarm[n].p = sqrt(swarm[n].p);


    }
}


void fuerza_vec(particle_t *swarmc,
        double *force_x, double *force_y, double *force_z,
        int N, int n, float **tforce, double L)
{
    *force_x = 0, *force_y = 0, *force_z = 0;
    int m;
    int a, b, c;
    
    /*
     * para cada particula necesito buscar todos sus vecinos y calcular la
     * fuerza a ellos en sus tres componentes, eso se guarda en fuerza_i
     */

    for( m = 0; m < N; m++){
        if ( m == n) continue;
        for(a=-1;a<2;a++){
            for(b=-1;b<2;b++){
                for(c=-1;c<2;c++){
                    vecf(swarmc, force_x, force_y, force_z, n, m, a, b ,c, L, tforce);
                }
            }
        }

    }


}

//calcula la fuerza entre la particula m y n
void vecf(particle_t *swarmc,double *force_x, double *force_y, double *force_z, int n, int m,
        int a, int b, int c, double L, float **tforce){
        double distancia;
        double dir_x, dir_y, dir_z;
        double force = 0;
        distancia = sqrt(pow(swarmc[m].x+a*L-swarmc[n].x,2)+pow(swarmc[m].y+b*L-swarmc[n].y,2)+pow(swarmc[m].z+c*L-swarmc[n].z,2));
    if(distancia<ROUT) {
        dir_x = -(swarmc[m].x+a*L-swarmc[n].x)/distancia;
        dir_y = -(swarmc[m].y+b*L-swarmc[n].y)/distancia;
        dir_z = -(swarmc[m].z+c*L-swarmc[n].z)/distancia;
       force = appforce(tforce, distancia);
//       printf("fuerza %f, distancia %f\n", force, distancia);
//       printf("dentro a %d b %d c %d\n", a, b,c);
       *force_x += dir_x * force;
       *force_y += dir_y * force;
       *force_z += dir_z * force;
    }

}
