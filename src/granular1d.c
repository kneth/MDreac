/*
 * granular1d - simulation of granular media in one dimension
 *
 * (C) Copyright 2012 by Kenneth Geisshirt <http://kenneth.geisshirt.dk/>
 * Released under GNU General Public License v2
 *
 * References:
 * [1] Computer Simulation of Liquids. M.P. Allen and D.J. Tildesley. Claredon
 *     Press, 1987.
 * [2] Non-Equilibrium Molecular Dynamics Simulations: Applications to
 *     Oscillating Chemical Reactions and Inelastic Colliding Particles. K.
 *     Geisshirt. Roskilde University, 1997.
 * [3] P. Eshuis et al. Physical Review E 80, 011302 (2009).
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
 
// particles
double *x;          // positions
double *v;          // velocities
 
// system
size_t N;           // number of particles
double rho;         // density
double T_init;      // initial temperature
double T_wall;      // temperature of left wall
double R;           // radius of particles
double eps;         // elasticity
size_t N_coll;      // number of collisions
 
// derived parameters
double L;           // length of system

// statistics
size_t N_pp;        // particle/particle collisions
size_t N_pl;        // collisions with left wall
size_t N_pr;        // collisions with right wall 
double time;        // master time

int dcompare(const void *p1, const void *p2) {
    double a = *(double *)p1;
    double b = *(double *)p2;
    if (a > b)
        return 1;
    else if (a < b)
        return -1;
    else
        return 0;
}

void Initialize(void) {
    double dL, sum = 0.0, sum_sq = 0.0;
    size_t i;
    
    x = (double *)calloc(N, sizeof(double));
    v = (double *)calloc(N, sizeof(double));
    
    L = (double)N/rho;
    dL = (L-4.0*R)/(double)N;
    for(i=0; i<N; ++i) {
        x[i]    = (double)(i+1)*dL;
        v[i]    = 2.0*drand48()-1.0;
        sum    += v[i];
    }

    for(i=0; i<(N-1); ++i) {
        if (x[i+1]-x[i] <= 2.0*R) {
            printf("Overlap\n");
            exit(-1);
        }
    }
        
    // adjust velocities
    sum /= (double)N;
#pragma omp parallel for
    for(i=0; i<N; ++i) {
        v[i] -= sum;
        sum_sq += v[i]*v[i];
    }
    
    // scale to temperature
    double sc = sqrt(T_init*(double)(N-1)/sum_sq);
#pragma omp parallel for
    for(i=0; i<N; ++i) {
        v[i] *= sc;
    }
    
    time = 0.0;
}
 
void Bump(void) {
    size_t i, which;
    double dt, dx, dv, tau, tmp;
    enum { P_NONE, P_P, P_LEFT, P_RIGHT } what;
    
    which = N+1;       // out of range indicating undefined
    dt    = HUGE_VAL;  // just some large value
    what  = P_NONE;
    
    // particle-particle collisions
    for(i=0; i<(N-2); ++i) {
        dv  = v[i+1]-v[i];
        dx  = x[i]-x[i+1]-2.0*R;
        tau = dx/dv;
        if (tau >= 0.0 && dt >= tau) {
            dt    = tau;
            which = i;
            what  = P_P;
        }
    }
    
    // left wall
    tau = (R-x[0])/v[0];
    if (tau >= 0.0 && dt >= tau) {
        dt    = tau;
        which = 0;
        what  = P_LEFT;
    }
    
    // right wall
    tau = (L-R-x[N-1])/v[N-1];
    if (tau >= 0.0 && dt >= tau) {
        dt    = tau;
        which = N-1;
        what  = P_RIGHT;
    }
    
    // move particles
#pragma omp parallel for
    for(i=0; i<N; ++i) {
        x[i] += v[i]*dt;
        if (x[i] > L) {
            printf("Moving right: %ld %ld %ld\n", i, which, what);
            exit(-1);
        }
        if (x[i] < 0.0) {
            printf("Moving left: %ld %ld %ld\n", i, which, what);
            exit(-1);
        }
    }
    
    // collision
    switch (what) {
    case P_P:
        tmp        = v[which]; // TODO: inelastic collision
        v[which]   = v[which+1];
        v[which+1] = tmp;
        N_pp++;
        break;
    case P_RIGHT:
        v[N-1] = -v[N-1]; // mirror
        N_pr++;
        break;
    case P_LEFT:
        v[0] = -v[0]; // TODO: thermostat
        N_pl++; 
        break;
    case P_NONE:
        printf("No collision found!\n");
        exit(-1);
    }
    
    time += dt;
}

void Usage(char *prgname) {
    printf("Usage: %s\n", prgname);
    printf("  -h     this text\n");
    printf("  -T     temperature at left wall\n");
    printf("  -t     initial temperature\n");
    printf("  -N     number of particles\n");
    printf("  -n     number of collision to simulate\n");
    printf("  -d     (number) density\n");
    printf("  -R     radius of particles\n");
    printf("  -e     elasticity\n");
    exit(-1);
}

void ParseCmdLine(int argc, char *argv[]) {
    int          c;
    char        *progname;
    extern char *optarg;

    progname = strdup(argv[0]);
    while ((c=getopt(argc, argv, "he:T:N:n:d:R:t:")) != EOF) {
        switch (c) {
        case 'h':
            Usage(progname);
            break;
        case 'e':
            eps = atof(optarg);
            break;
        case 'T':
            T_wall = atof(optarg);
            break;
        case 't':
            T_init = atof(optarg);
            break;
        case 'N':
            N = strtoul(optarg, NULL, 0);
            break;
        case 'n':
            N_coll = strtoul(optarg, NULL, 0);
            break;
        case 'd':
            rho = atof(optarg);
            break;
        case 'R':
            R = atof(optarg);
            break;
        }
    }
}

void PrintReport(size_t n) {
    size_t i;
    double x_mean = 0.0;
    double v_sq   = 0.0;
    double v_mean = 0.0;

#pragma omp parallel for
    for(i=0; i<N; ++i) {
        x_mean += x[i];
        v_mean += v[i];
        v_sq   += v[i]*v[i];
    }
    x_mean /= (double)N;
    
    printf("%ld %e %e %e %e %ld %ld %ld\n", n, time, x_mean, v_mean, v_sq, N_pp, N_pl, N_pr);
}

void PrintConfig(size_t n) {
    size_t i;
    for(i=0; i<N; ++i) {
        printf("%ld %e %e\n", i, x[i], v[i]);
    }
}

int main(int argc, char *argv[]) {
    size_t i;
    
    ParseCmdLine(argc, argv);
    Initialize();
    PrintConfig(0);
    N_pp = N_pl = N_pr = 0;
    for(i=0; i<N_coll; ++i) {
        Bump();
        PrintReport(i);
    }
//    PrintConfig(i);
    exit(0); 
}
