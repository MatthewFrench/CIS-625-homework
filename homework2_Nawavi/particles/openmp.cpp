#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"
#include <pthread.h>



#if !defined density
#define density 0.0005
#endif

#if !defined mass
#define mass 0.01
#endif

#if !defined cutoff
#define cutoff 0.01
#endif

#if !defined min_r
#define min_r (cutoff/100)
#endif

#if !defined dt
#define dt      0.0005
#endif

//
// particle data structure
//
struct particleNode
{
    particle_t *particle;
    int gridX, gridY;
    int index;
    struct particleNode *next, *prev;
};
typedef struct particleNode ParticleNode;

//My updated force function that updates both particles at the same time
void static inline apply_forceBoth( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
    if (r2 != 0)
    {
        if (r2/(cutoff*cutoff) < *dmin * (*dmin))
            *dmin = sqrt(r2)/cutoff;
        (*davg) += sqrt(r2)/cutoff;
        (*navg) ++;
    }

    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;

    double coefdx = coef * dx;
    double coefdy = coef * dy;

    particle.ax += coefdx;
    particle.ay += coefdy;

    neighbor.ax -= coefdx;
    neighbor.ay -= coefdy;
}
//This is to debug our grid to make sure particle nodes are where they should be
void printGrid(ParticleNode ***grid, int gridSize) {
    printf("***PRINT GRID %dx%d***\n", gridSize, gridSize);
    ParticleNode* particleNode;
    for (int x = 0; x < gridSize; x++) {
        printf("[");
        for (int y = 0; y < gridSize; y++) {
            particleNode = grid[x][y];
            printf("%d, %d=", x, y);
            bool correct = true;
            int nodeCount = 0;
            while (particleNode != NULL) {
                nodeCount++;
                if (!(particleNode->gridX == x && particleNode->gridY == y)) {
                    printf("(invalid location: %d, %d)\n", particleNode->gridX, particleNode->gridY);
                    correct = false;
                }
                particleNode = particleNode->next;
            }
            if (correct) {
                printf("correct(%d nodes), ", nodeCount);
            } else  {
                particleNode = grid[x][y];
                while (particleNode != NULL) {
                    printf("|Node Details: %d <-- %d, %d --> %d|\n", particleNode->prev, particleNode->gridX, particleNode->gridY, particleNode->next);
                    particleNode = particleNode->next;
                }
                printf("incorrect(%d nodes), ", nodeCount);
            }
        }
        printf("]\n");
    }
}

inline void addParticleNodeToGrid(int particleIndex, int x, int y, ParticleNode* addNode, ParticleNode ***grid) {
    //Get the first node at the current grid spot
    ParticleNode* node = grid[x][y];
    if (node == NULL) { //Grid is empty so just set the node in the grid
        addNode->gridX = x;
        addNode->gridY = y;
        grid[x][y] = addNode;
    } else { //Add the particle node to the grid but place it ordered in the linked list based on index to make memory lookup slightly faster
        addNode->gridX = x;
        addNode->gridY = y;
        ParticleNode* firstNode = node;
        ParticleNode* pastNode = NULL;
        while (node->index < particleIndex) {
            pastNode = node;
            node = node->next;
            if (node == NULL) break;
        }
        //Node should equal either null or a node we need to place our particle node in front of
        if (node == firstNode) {
            //Just push back the first node
            addNode->next = node;
            node->prev = addNode;
            grid[x][y] = addNode;
        } else if (node != NULL) {
            //Put the particle node between two nodes
            ParticleNode* leftNode = node->prev;
            ParticleNode* rightNode = node;
            addNode->prev = leftNode;
            leftNode->next = addNode;
            addNode->next = rightNode;
            rightNode->prev = addNode;
        } else {
            //Put the particle node at the end
            ParticleNode* leftNode = pastNode;
            addNode->prev = leftNode;
            leftNode->next = addNode;
        }
    }
}

inline void removeParticleNodeFromGrid(int x, int y, ParticleNode* node, ParticleNode ***grid) {
    //Get the neighbor nodes in the same grid as this particle node
    ParticleNode* nextNode = node->next;
    ParticleNode* prevNode = node->prev;
    //We need to remove the Node from it's neighbors and from the grid
    if (nextNode != NULL && prevNode != NULL) { //Particle is in the middle of a linked list
        //Splice the particle node out by setting our neighbors to be next to each other
        nextNode->prev = prevNode;
        prevNode->next = nextNode;
    } else if (nextNode != NULL && prevNode == NULL) { //Particle is first in the linked list and is connected to another node
        //Set the second node to be the first on the grid and remove the second node's reference to this particle node
        nextNode->prev = NULL;
        grid[x][y] = nextNode;
    } else if (nextNode == NULL && prevNode != NULL) { //Particle is the last in the linked list
        //Remove the neighbor's connection to this particle node
        prevNode->next = NULL;
    } else if (nextNode == NULL && prevNode == NULL) { //Grid is empty except for this particle node
        //Set the grid at this location to NULL cause the particle node is being removed and it is the only node there
        grid[x][y] =  NULL;
    }
    //Particle node was removed so remove references to neighbors
    node->next = NULL;
    node->prev = NULL;
}

struct searchArea {
    int startX, startY, endX, endY;
};
typedef struct searchArea SearchArea;

inline SearchArea getSearchArea(int x, int y, int searchArea, int gridSize) {
    SearchArea search;
    search.startX = x - searchArea;
    if (search.startX < 0) search.startX = 0;
    search.startY = y - searchArea;
    if (search.startY < 0) search.startY = 0;
    search.endX = x + searchArea;
    if (search.endX >= gridSize) search.endX = gridSize-1;
    search.endY = y + searchArea;
    if (search.endY >= gridSize) search.endY = gridSize-1;

    return search;
}
/*
 void swap(void *i, void *j) {
 void t = *i;
 *i = *j;
 *j = t;
 }*/

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg,nabsavg=0,numthreads;
    double dmin, absmin=1.0,davg,absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    particle_t *particlesCopy = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    for (int i = 0; i < n; i++) {
        particlesCopy[i] = particles[i];
    }




    //Create the field grid
    double blockSize = cutoff; //Use to be cutoff
    int searchArea = 2;
    double fieldSize = sqrt(density * n);
    int gridSize = ceil(fieldSize / blockSize);
    //Create particle grid
    ParticleNode ***grid = (ParticleNode ***)malloc(gridSize * sizeof(ParticleNode **));
    for (int i=0; i<gridSize; i++)
        grid[i] = (ParticleNode **)malloc(gridSize * sizeof(ParticleNode*));

    ParticleNode ***gridCopy = (ParticleNode ***)malloc(gridSize * sizeof(ParticleNode **));
    for (int i=0; i<gridSize; i++)
        gridCopy[i] = (ParticleNode **)malloc(gridSize * sizeof(ParticleNode*));
    //Create particle node array
    ParticleNode **particleNodes = (ParticleNode**) malloc( n * sizeof(ParticleNode*) );
    ParticleNode **particleNodesCopy = (ParticleNode**) malloc( n * sizeof(ParticleNode*) );

    //Clear the particle grid
    for (int x = 0; x < gridSize; x++) {
        for (int y = 0; y < gridSize; y++) {
            grid[x][y] = NULL;
            gridCopy[x][y] = NULL;
        }
    }

    //Sort particles in linked lists in the grid
    for (int i = 0; i < n; i++) {
        particle_t *particle = &particles[i];
        //Get location in grid
        int locationX = floor(particle->x / blockSize);
        int locationY = floor(particle->y / blockSize);
        //Get node at that grid location
        ParticleNode* node = grid[locationX][locationY];
        //Create new node for particle we are adding
        ParticleNode* newNode = (ParticleNode*) malloc(sizeof(ParticleNode));
        newNode->particle = particle;
        newNode->next = NULL;
        newNode->prev = NULL;
        newNode->index = i;

        addParticleNodeToGrid(i, locationX, locationY, newNode, grid);
        //Set the new particle node in the particle nodes array
        particleNodes[i] = newNode;
    }
    ///Sort copies
    for (int i = 0; i < n; i++) {
        particle_t *particle = &particlesCopy[i];
        //Get location in grid
        int locationX = floor(particle->x / blockSize);
        int locationY = floor(particle->y / blockSize);
        //Get node at that grid location
        ParticleNode* node = gridCopy[locationX][locationY];
        //Create new node for particle we are adding
        ParticleNode* newNode = (ParticleNode*) malloc(sizeof(ParticleNode));
        newNode->particle = particle;
        newNode->next = NULL;
        newNode->prev = NULL;
        newNode->index = i;

        addParticleNodeToGrid(i, locationX, locationY, newNode, gridCopy);
        //Set the new particle node in the particle nodes array
        particleNodesCopy[i] = newNode;
    }

    //Simulate

    //Create the locks for the grid

    /*
     omp_lock_t **gridLocks = (omp_lock_t **)malloc(gridSize * sizeof(omp_lock_t *));
     for (int i=0; i<gridSize; i++) {
     gridLocks[i] = (omp_lock_t *)malloc(gridSize * sizeof(omp_lock_t));
     for (int j=0; j<gridSize; j++) {
     omp_init_lock(&gridLocks[i][j]);
     }
     }*/
    pthread_spinlock_t** gridLocks;
    //Create grid locks
    gridLocks = (pthread_spinlock_t **)malloc(gridSize * sizeof(pthread_spinlock_t *));
    for (int i=0; i<gridSize; i++) {
        gridLocks[i] = (pthread_spinlock_t *)malloc(gridSize * sizeof(pthread_spinlock_t));
        for (int j = 0; j < gridSize; j++) {
            pthread_spin_init(&gridLocks[i][j], PTHREAD_PROCESS_PRIVATE);
        }
    }
#pragma omp parallel
    numthreads = omp_get_num_threads();

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    #pragma omp parallel
    {
        for( int step = 0; step < 1000; step++ )
        {
            //{
            navg = 0;
            davg = 0.0;
            dmin = 1.0;

                //Loop through all particles
#pragma omp for reduction (+:navg) reduction(+:davg)
                for( int i = 0; i < n; i++ )
                {
                    //Get the particle Node and particle to simulate
                    ParticleNode* particleNodeI = particleNodes[i];
                    particle_t *particleI = particleNodeI->particle;

                    //Create a grid search area for particle collision
                    int locationX = particleNodeI->gridX;
                    int locationY = particleNodeI->gridY;

                    SearchArea search = getSearchArea(locationX, locationY, searchArea, gridSize);

                    for (int x = search.startX; x <= search.endX; x++) {
                        for (int y = search.startY; y <= search.endY; y++) {
                            //pthread_spin_lock(&(gridLocks[x][y]));
                            ParticleNode* startingNode = grid[x][y];
                            while (startingNode != NULL) { //Go down the entire linked list on the grid location and calculate particle forces
                                //if (startingNode->index > particleNodeI->index) {
                                //  apply_forceBoth( *particleI, *(startingNode->particle),&dmin,&davg,&navg);
                                //}
                                apply_force( *particleI, *(startingNode->particle),&dmin,&davg,&navg);
                                //Set the node to the next node so we continue down the chain of nodes
                                startingNode = startingNode->next;
                            }
                            //pthread_spin_unlock(&(gridLocks[x][y]));
                        }
                    }
                }
#pragma omp for
            for( int i = 0; i < n; i++ ) {
                
                ParticleNode* particleNodeI = particleNodes[i];
                particle_t *particleI = particleNodeI->particle;
                
                //Move a copy of the particle
                ParticleNode* particleNodeICopy = particleNodesCopy[i];
                particle_t *particleICopy = particleNodeICopy->particle;
                particleICopy->ax = particleI->ax;
                particleICopy->ay = particleI->ay;
                move( *particleICopy );
                particleICopy->ax = particleICopy->ay = 0;
                //Place the copy of the particle node into the right spot in the copy of the grid array
                int oldLocationX = particleNodeICopy->gridX;
                int oldLocationY = particleNodeICopy->gridY;
                int newLocationX = floor(particleICopy->x / blockSize);
                int newLocationY = floor(particleICopy->y / blockSize);
                //Update the grid cause the particle moved
                if (oldLocationX != newLocationX || oldLocationY != newLocationY) {
                    pthread_spin_lock(&(gridLocks[oldLocationX][oldLocationY]));
                    removeParticleNodeFromGrid(oldLocationX, oldLocationY, particleNodeICopy, gridCopy);
                    pthread_spin_unlock(&(gridLocks[oldLocationX][oldLocationY]));
                    
                    pthread_spin_lock(&(gridLocks[newLocationX][newLocationY]));
                    addParticleNodeToGrid(i, newLocationX, newLocationY, particleNodeICopy, gridCopy);
                    pthread_spin_unlock(&(gridLocks[newLocationX][newLocationY]));
                }
            }
                #pragma omp for
                for( int i = 0; i < n; i++ ) {
                  particles[i] = particlesCopy[i];
                }

#pragma omp master
{
            //Swap
            ParticleNode **particleNodesTemp = particleNodesCopy;
            particleNodesCopy = particleNodes;
            particleNodes = particleNodesTemp;
            //Swap
            particle_t *particlesTemp = particlesCopy;
            particlesCopy = particles;
            particles = particlesTemp;
            //Swap
            ParticleNode ***gridTemp = gridCopy;
            gridCopy = grid;
            grid = gridTemp;
          }
            /*
             #pragma omp for
             for( int i = 0; i < n; i++ )
             {
             //Get the particle Node and particle to simulate
             ParticleNode* particleNodeI = particleNodes[i];
             particle_t *particleI = particleNodeI->particle;

             //Move the particle
             move( particles[i] );
             particleI->ax = particleI->ay = 0;

             //Find the old and new location in the grid of the particle
             int oldLocationX = particleNodeI->gridX;
             int oldLocationY = particleNodeI->gridY;
             int newLocationX = floor(particleI->x / blockSize);
             int newLocationY = floor(particleI->y / blockSize);
             //Update the grid cause the particle moved
             if (oldLocationX != newLocationX || oldLocationY != newLocationY) {
             pthread_spin_lock(&(gridLocks[oldLocationX][oldLocationY]));
             removeParticleNodeFromGrid(oldLocationX, oldLocationY, particleNodeI, grid);
             pthread_spin_unlock(&(gridLocks[oldLocationX][oldLocationY]));

             pthread_spin_lock(&(gridLocks[newLocationX][newLocationY]));
             addParticleNodeToGrid(i, newLocationX, newLocationY, particleNodeI, grid);
             pthread_spin_unlock(&(gridLocks[newLocationX][newLocationY]));
             }
             }
             }*/


            if( find_option( argc, argv, "-no" ) == -1 )
            {
                //
                //  compute statistical data
                //
                //Run this section of code on the master thread
                #pragma omp master
                if (navg) {
                    absavg += davg/navg;
                    nabsavg++;
                }
                //Can only be ran on a single thread at a time
                #pragma omp critical
                if (dmin < absmin) absmin = dmin;

                //
                //  save if necessary
                //
                //Run this section of code on the master thread
                #pragma omp master
                if( fsave && (step%SAVEFREQ) == 0 )
                    save( fsave, n, particles );
                //}
            }
        }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
        if (nabsavg) absavg /= nabsavg;
        //
        //  -the minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
        if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    /*
     for (int i=0; i<gridSize; i++) {
     for (int j=0; j<gridSize; j++) {
     omp_destroy_lock(&gridLocks[i][j]);
     }
     free(gridLocks[i]);
     }
     free(gridLocks);
     */
    for (int i=0; i<gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            pthread_spin_destroy(&gridLocks[i][j]);
        }
        //free(gridLocks[i]);
    }
    free(gridLocks);

    for (int i=0; i<gridSize; i++) {
        free(grid[i]);
    }
    free(grid);
    for (int i=0; i<gridSize; i++) {
        free(gridCopy[i]);
    }
    free(gridCopy);

    for (int i = 0; i < n; i++) {
        free(particleNodes[i]);
    }
    free(particleNodes);
    for (int i = 0; i < n; i++) {
        free(particleNodesCopy[i]);
    }
    free(particleNodesCopy);

    if( fsum )
        fclose( fsum );

    free( particles );
    free( particlesCopy );
    if( fsave )
        fclose( fsave );

    return 0;
}

