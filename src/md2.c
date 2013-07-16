/*******************************************************************************
 *
 * md2 - molecular dynamics simulation of mixtures of simple liquids in
 *       two dimensions
 *
 * Example:
 *   ./md2 -T 2.0 -d 0.8 -t 10000 -a 512 -b 512 -z 4.0 -r 1.5 -R 1.5 \
 *     -o test.data -c test.conf -C final.conf -H 0.0005 -Q 0.5
 *
 * (C) Copyright 2008-2013 by Kenneth Geisshirt <http://kenneth.geisshirt.dk/>
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

/*** Parameters             ***/
double         T;           // temperature
double         rho;         // density
unsigned long  nA, nB;      // number of particles
double         Rc_even;     // cut-off radius
double         Rc_odd;
double         Rbuf;        // buffer zone
unsigned long  nsteps;      // number of time steps
char          *init_conf;   // initial configuration
char          *final_conf;  // final configuration
char          *outname;     // output file name
double         Q;           // thermostat strength
double         dt;          // length of time step

/*** Calculated parameters  ***/
double         L;           // length of the simulation box
unsigned long  nCells;      // cells in each direction


/*** Global variables       ***/
double        *rx, *ry;     // positions
double        *vx, *vy;     // velovities
double        *fx, *fy;     // forces
double         eta;         // thermostat
double         Epot;        // potential energy
double         Ekin;        // kinetic energy
double         P;           // pressure
unsigned long *map;
unsigned long *head, *list; // linked lists
unsigned long  Np;          // number of interacting pairs
unsigned long *ipair, *jpair; // interacting pairs


void Usage(char *prgname) {
    printf("Usage: %s\n", prgname);
    printf("  -h     this text\n");
    printf("  -T     temperature\n");
    printf("  -d     density\n");
    printf("  -t     number of time steps\n");
    printf("  -a     number of A particles\n");
    printf("  -b     number of B particles\n");
    printf("  -c     initial configuration\n");
    printf("  -C     final configuration\n");
    printf("  -Q     thermostat coupling\n");
    printf("  -o     output (temperature, pressure, etc.)\n");
    printf("  -r     cut-off radius (even pairs)\n");
    printf("  -R     cut-off radius (odd pairs)\n");
    printf("  -z     buffer zone for interacting particles\n");
    printf("  -H     length of time step\n");
    exit(0);
}

void ReadParameters(int argc, char *argv[]) {
    int          c;
    char        *progname;
    extern char *optarg;

    progname = strdup(argv[0]);
    while ((c=getopt(argc, argv, "hT:d:t:a:b:c:C:Q:o:r:R:H:z:")) != EOF) {
        switch (c) {
            case 'h':
                Usage(progname);
                break;
            case 'T':
                T = atof(optarg);
                break;
            case 'd':
                rho = atof(optarg);
                break;
            case 't':
                nsteps = strtoul(optarg, NULL, 0);
                break;
            case 'a':
                nA = strtoul(optarg, NULL, 0);
                break;
            case 'b':
                nB = strtoul(optarg, NULL, 0);
                break;
            case 'c':
                init_conf = strdup(optarg);
                break;
            case 'C':
                final_conf = strdup(optarg);
                break;
            case 'Q':
                Q = atof(optarg);
                break;
            case 'o':
                outname = strdup(optarg);
                break;
            case 'r':
                Rc_even = atof(optarg);
                break;
            case 'R':
                Rc_odd = atof(optarg);
                break;
            case 'H':
                dt = atof(optarg);
                break;
            case 'z':
                Rbuf = atof(optarg);
                break;
            default:
                fprintf(stderr, "Warning: option %c ignored.\n", c);
                break;
        }
    }
}
void ReadConfiguration(void) {
    double         sumx = 0.0, sumy = 0.0;
    double         sum = 0.0;
    double         dof, sc, dL;
    unsigned long  i, j, m, n, k;

    n = nA+nB;
    m = (unsigned long)sqrt((double)n);
    dL = L/(double)m;
    for(i=0; i<m; i++) {
        for(j=0; j<m; j++) {
            k = i*m+j;
            rx[k] = -0.5*L+dL*(double)i;
            ry[k] = -0.5*L+dL*(double)j;
            vx[k] = 2.0*drand48()-1.0; // fairly low quality RNG
            vy[k] = 2.0*drand48()-1.0;
            sumx += vx[k];
            sumy += vy[k];
            sum += vx[k]*vx[k]+vy[k]*vy[k];
        }
    }

    /* adjust velocities */
    sumx /= (double)(nA+nB);
    sumy /= (double)(nA+nB);
#pragma omp parallel for
    for(i=0; i<(nA+nB); i++) {
        vx[i] -= sumx;
        vy[i] -= sumy;
    }

    /* scale to temperature */
    dof = 2.0*((double)(nA+nB))-2.0;
    sc = sqrt(dof*T/sum);
#pragma omp parallel for
    for(i=0; i<(nA+nB); i++) {
        vx[i] *= sc;
        vy[i] *= sc;
    }
}

void WriteConfiguration(char *filename) {
    FILE            *outfile;
    unsigned long    i;

    outfile = fopen(filename, "w");

    /* first parameters */
    fprintf(outfile, "%ld\n%ld\n", nA, nB);
    fprintf(outfile, "%e\n%e\n", rho, T);

    /* then positions and velocities */
    for(i=0; i<(nA+nB); i++) {
        fprintf(outfile, "%e %e %e %e\n", rx[i], ry[i], vx[i], vy[i]);
    }

    fclose(outfile);
}

double Rcut(unsigned long i, unsigned long j) {
    if (((i < (nA-1)) && (j < (nA-1))) || ((i > nA) && (j > nA))) {
        return Rc_even;
    } else {
        return Rc_odd;
    }
}
double pbc(double d) {
    double len2 = 0.5*L;
    if (d > len2) return -len2;
    if (d < -len2) return len2;
    return 0.0;
}

unsigned long Icell(unsigned long i, unsigned long j) {
    if (i == -1) {
        i = nCells-1;
    } else if (i == nCells) {
        i = 0;
    }
    if (j == -1) {
        j = nCells-1;
    } else if (j == nCells) {
        j = 0;
    }
    return i*nCells+j;
}

void ApplyPerBoundaries(void) {

    unsigned long i;
    for(i=0; i<(nA+nB); i++) {
        rx[i] += pbc(rx[i]);
        ry[i] += pbc(ry[i]);
    }
}
void Initialize(void) {
    unsigned long i, j, imap, m;
    unsigned long n = nA+nB;
    double Lbuf;

    printf("md2 - (C) Copyright 2008 by Kenneth Geisshirt\n");
    printf("Parameters:\n");
    printf("  Number of particles: %ld + %ld = %ld\n", nA, nB, n);
    printf("  Temperature:         %e\n", T);
    printf("  Density:             %e\n", rho);
    printf("  Buffer zone:         %e\n", Rbuf);
    printf("File names:\n");
    printf("  init. configuration: %s\n", init_conf);
    printf("  final configuration: %s\n", final_conf);
    printf("  output             : %s\n", outname);
    printf("Derived parameters:\n");

    eta = 0.0;

    /* allocate memory */
    rx = (double *)calloc(n, sizeof(double));
    ry = (double *)calloc(n, sizeof(double));
    vx = (double *)calloc(n, sizeof(double));
    vy = (double *)calloc(n, sizeof(double));
    fx = (double *)calloc(n, sizeof(double));
    fy = (double *)calloc(n, sizeof(double));

    /* calculate additional parameters */
    L = sqrt(((double)n)/rho);
    printf("  Length:              %e\n", L);

    /* cell structures for force optimizations */
    nCells = (int)floor(L/(2.0*Rbuf));
    Lbuf = L/(double)nCells;
    printf("  Number of cells:     %ld\n", nCells);
    printf("  Length of cell:      %e\n", Lbuf);
    head = (unsigned long *)calloc(nCells*nCells, sizeof(unsigned long));
    list = (unsigned long *)calloc(n, sizeof(unsigned long));
    map = (unsigned long *)calloc(4*nCells*nCells, sizeof(unsigned long));

    for(j=0; j<nCells; j++) {
        for(i=0; i<nCells; i++) {
            imap = 4*Icell(i, j);
            map[imap]   = Icell(i+1, j-1);
            map[imap+1] = Icell(i+1,   j);
            map[imap+2] = Icell(i+1, j+1);
            map[imap+3] = Icell(i,   j+1);
        }
    }

    m = (unsigned long)(1.5*floor(1.0+rho*Lbuf*Lbuf));
    m = 9*nCells*nCells*m*m;
    printf("  Max. pairs:          %ld\n", m);
    ipair = (unsigned long *)calloc(m, sizeof(unsigned long));
    jpair = (unsigned long *)calloc(m, sizeof(unsigned long));
}

void NoseHoover(void) {

    double         eta1, eta2, K, sum;
    double         dof = 2.0*(double)(nA+nB)-2.0;
    double         dt48 = 48.0*dt;
    unsigned long  i;

    sum = 0.0;
    K = 0.5*dt*eta;
    eta1 = 1.0-K;
    eta2 = 1.0/(1.0+K);
#pragma omp parallel for
    for(i=0; i<(nA+nB); i++) { /* update velocities */
        vx[i] = (vx[i]*eta1+dt48*fx[i])*eta2;
        vy[i] = (vy[i]*eta1+dt48*fy[i])*eta2;
    }
#pragma omp parallel for
    for(i=0; i<(nA+nB); i++) { /* update positions */
        rx[i] += vx[i]*dt;
        ry[i] += vy[i]*dt;
        sum += vx[i]*vx[i]+vy[i]*vy[i];
    }
    eta += dt*(sum-dof*T)/Q; /* update thermostat */
    Ekin = 0.5*sum;
}


void PutInBox(void) {
    unsigned long i, c;
    double        len2 = 0.5*L;
    double        cL = L/(double)nCells;

#pragma omp parallel for
    for(i=0; i<nCells*nCells; i++) {
        head[i] = -1;
    }
#pragma omp parallel for
    for(i=0; i<(nA+nB); i++) {
        c = (unsigned long)((rx[i]+len2)/cL)*nCells
            + (unsigned long)((ry[i]+len2)/cL);
        list[i] = head[c];
        head[c] = i;
    }
}

void MakeList(void) {
    unsigned long i, j, k, n, jcell;
    double        rxi, ryi;
    double        rxij, ryij, r2;

    Np = 0;
    for(k=0; k<(nCells*nCells); k++) {
        i = head[k];
        rxi = rx[i];
        ryi = ry[i];
        while (i > 0) {
            j = head[k];
            while (j > 0) {
                rxij = rxi-rx[j];
                ryij = ryi-ry[j];
                rxij += pbc(rxij);
                ryij += pbc(ryij);
                r2 = rxij*rxij+ryij*ryij;
                if (r2 < (Rbuf*Rbuf)) {
                    ipair[Np] = i;
                    jpair[Np] = j;
                    Np++;
                }
                j = list[j];
            }
            /* now neighboring cells */
            for(n=0; n<4; n++) {
                jcell = map[4*k+n];
                j = head[jcell];
                while (j > 0) {
                    rxij = rxi-rx[j];
                    ryij = ryi-ry[j];
                    rxij += pbc(rxij);
                    ryij += pbc(ryij);
                    r2 = rxij*rxij+ryij*ryij;
                    if (r2 < (Rbuf*Rbuf)) {
                        ipair[Np] = i;
                        jpair[Np] = j;
                        Np++;
                    }
                    j = list[j];
                }
            }
            i = list[i];
        }
    }
}


void ComputeForces(void) {
    unsigned long    i, j, k;
    double           pot;
    double           rxij, ryij, r2, r6, Rc;

#pragma omp parallel for
    for(i=0; i<(nA+nB); i++) { /* reset forces */
        fx[i] = 0.0;
        fy[i] = 0.0;
    }
    Epot = 0.0;
    P    = 0.0;
#pragma omp parallel for
    for(k=0; k<Np; k++) { /* loop over all pairs */
        i = ipair[k];
        j = jpair[k];
        rxij = rx[i]-rx[j];
        ryij = ry[j]-ry[j];
        rxij += pbc(rxij);
        ryij += pbc(ryij);
        r2 = rxij*rxij+ryij*ryij;
        Rc = Rcut(i, j);
        if (r2 > Rc*Rc) {
            r2     = 1.0/r2;
            r6     = r2*r2*r2;
            Epot  += r6*(r6-1.0);
            P     += r6*(r6-0.5);
            pot    = r2*r6*(r6-0.5);
            fx[i] += rxij*pot;
            fy[i] += ryij*pot;
            fx[j] -= rxij*pot;
            fy[j] -= ryij*pot;
        }
    }
    Epot *= 4.0;
}


int main(int argc, char *argv[]) {

    unsigned long i;
    FILE *statfile;

    ReadParameters(argc, argv);
    Initialize();
    statfile = fopen(outname, "w");
    ReadConfiguration();
    for(i=0; i<nsteps; i++) { /* a number of single time steps */
        ApplyPerBoundaries();
        PutInBox();
        MakeList();
        ComputeForces();
        NoseHoover();
        if ((nsteps % 1000) == 0) {
            fprintf(statfile, "%ld %e %e %e %e\n", i, Ekin, Epot, P, eta);
        }
    }
    fclose(statfile);
    WriteConfiguration(final_conf);
    return 1;
}
