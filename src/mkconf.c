/*******************************************************************************
 *
 * mkconf - generate configuration for md2
 *
 * (C) Copyright 2008 by Kenneth Geisshirt <http://kenneth.geisshirt.dk/>
 * Released under GNU General Public License v2 or later.
 * 
 * References:
 * [1] Computer Simulation of Liquids. M.P. Allen and D.J. Tildesley. Claredon
 *     Press, 1987.
 * [2] Non-Equilibrium Molecular Dynamics Simulations: Applications to
 *     Oscillating Chemical Reactions and Inelastic Colliding Particles. K.
 *     Geisshirt. Roskilde University, 1997.
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

int main(int argc, char *argv[]) {
    extern char *optarg;
    char c;
    FILE *outfile;
    
    
    double rho = 0.8;
    int    n   = 1024;

    double rx, ry;
    double vx, vy;
    int m, i, j;
    double L, dL;

    while ((c=getopt(argc, argv, "n:d:")) != EOF) {
        switch (c) {
        case 'n':
            n = atoi(optarg);
            break;
        case 'd':
            rho = atof(optarg);
            break;
        default:
            fprintf(stderr, "Unknown option\n");
            break;
        }
    }
    L = sqrt(rho/(double)n);
    m = (int)sqrt((double)n);
    dL = L/m;
    outfile = fopen("test.conf", "w");
    for(i=0; i<m; i++) {
        rx = -0.5*L+dL*(double)i;
        for(j=0; j<m; j++) {
            ry = -0.5*L+dL*(double)j;
            vx = 2.0*drand48()-1.0; /* fairly low quality RNG */
            vy = 2.0*drand48()-1.0;
            fprintf(outfile, "%e %e %e %e\n", rx, ry, vx, vy);
        }
    }
    fclose(outfile);
    return 1;
}
