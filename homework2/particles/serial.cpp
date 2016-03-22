#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

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


// changed
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

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
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
    set_size( n );
    init_particles( n, particles );

    //Create the field grid
    double blockSize = cutoff; //Use to be cutoff
    int searchArea = 4;//33
    double fieldSize = sqrt(density * n);
    int gridSize = ceil(fieldSize / blockSize);
    //Create particle grid
    ParticleNode ***grid = (ParticleNode ***)malloc(gridSize * sizeof(ParticleNode **));
    for (int i=0; i<gridSize; i++)
        grid[i] = (ParticleNode **)malloc(gridSize * sizeof(ParticleNode*));
    //Create particle node array
    ParticleNode **particleNodes = (ParticleNode**) malloc( n * sizeof(ParticleNode*) );

    //Clear the particle grid
    for (int x = 0; x < gridSize; x++) {
        for (int y = 0; y < gridSize; y++) {
            grid[x][y] = NULL;
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
        if (node == NULL) { //Add the new particle node
            newNode->gridX = locationX;
            newNode->gridY = locationY;
            grid[locationX][locationY] = newNode;
        } else { //Add the new node to the grid but place it ordered in the linked list based on index to make memory lookup slightly faster
            newNode->gridX = locationX;
            newNode->gridY = locationY;
            ParticleNode* firstNode = node;
            ParticleNode* pastNode = NULL;
            while (node->index < i) {
                pastNode = node;
                node = node->next;
                if (node == NULL) break;
            }
            //Node should equal either null or a node we need to place our new node in front of
            if (node == firstNode) {
                //Just push back the first node
                newNode->next = node;
                node->prev = newNode;
                grid[locationX][locationY] = newNode;
            } else if (node != NULL) {
                //Put the new node between two nodes
                ParticleNode* leftNode = node->prev;
                ParticleNode* rightNode = node;
                newNode->prev = leftNode;
                leftNode->next = newNode;
                newNode->next = rightNode;
                rightNode->prev = newNode;
            } else {
                //Put the new node at the end
                ParticleNode* leftNode = pastNode;
                newNode->prev = leftNode;
                leftNode->next = newNode;
            }
        }
        //Set the new particle node in the particle nodes array
        particleNodes[i] = newNode;
    }

    //Simulate
    /* How simulation works
     Past: Loop through all the particles and loop through all particles against that particle to calculate forces, then loop through particles and move the particle.
        for particles
            clear forces particle
        for particles
            for particles
                applyForce particle
        for particles
            move particle
     Now: Loop through all the particles, move the particle, clear forces and loop through nearby grid blocks and calculate forces of particles in those blocks.
        for particles
            for search grid x
                for search grid y
                    for particle linked list
                        applyForce particle
            move
            update particle node in grid
            clear forces

     */
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        //Loop through all particles
        for( int i = 0; i < n; i++ )
        {
            //Get the particle Node and particle to simulate
            ParticleNode* particleNodeI = particleNodes[i];
            particle_t *particleI = particleNodeI->particle;

            //Create a grid search area for particle collision
            int oldLocationX = particleNodeI->gridX;
            int oldLocationY = particleNodeI->gridY;
            //Only search in a 2 block radius around the block the particle node is currently in
            int searchStartX = oldLocationX - searchArea;
            if (searchStartX < 0) searchStartX = 0;
            int searchStartY = oldLocationY - searchArea;
            if (searchStartY < 0) searchStartY = 0;
            int searchEndX = oldLocationX + searchArea;
            if (searchEndX >= gridSize) searchEndX = gridSize-1;
            int searchEndY = oldLocationY + searchArea;
            if (searchEndY >= gridSize) searchEndY = gridSize-1;

            for (int x = searchStartX; x <= searchEndX; x++) {
                for (int y = searchStartY; y <= searchEndY; y++) {
                    ParticleNode* startingNode = grid[x][y];
                    while (startingNode != NULL) { //Go down the entire linked list on the grid location and calculate particle forces
                        //if (startingNode->index > particleNodeI->index) {
                            //apply_forceBoth( *particleI, *(startingNode->particle),&dmin,&davg,&navg);
                        apply_force( particles[i], *(startingNode->particle),&dmin,&davg,&navg);
                        //}
                        //Set the node to the next node so we continue down the chain of nodes
                        startingNode = startingNode->next;
                    }
                }
            }
        }
        for( int i = 0; i < n; i++ )
        {
            //Get the particle Node and particle to simulate
            ParticleNode* particleNodeI = particleNodes[i];
            particle_t *particleI = particleNodeI->particle;

            //Move the particle
            move( particles[i] );
        }
        for( int i = 0; i < n; i++ )
        {
            //Get the particle Node and particle to simulate
            ParticleNode* particleNodeI = particleNodes[i];
            particle_t *particleI = particleNodeI->particle;

            //Find the old and new location in the grid of the particle
            int oldLocationX = particleNodeI->gridX;
            int oldLocationY = particleNodeI->gridY;
            int newLocationX = floor(particleI->x / blockSize);
            int newLocationY = floor(particleI->y / blockSize);
            //Update the grid cause the particle moved
            if (oldLocationX != newLocationX || oldLocationY != newLocationY) {
                //Get the neighbor nodes in the same grid as this particle node
                ParticleNode* nextNode = particleNodeI->next;
                ParticleNode* prevNode = particleNodeI->prev;
                //We need to remove the Node from it's neighbors and from the grid
                if (nextNode != NULL && prevNode != NULL) { //Particle is in the middle of a linked list
                    //Splice the particle node out by setting our neighbors to be next to each other
                    nextNode->prev = prevNode;
                    prevNode->next = nextNode;
                } else if (nextNode != NULL && prevNode == NULL) { //Particle is first in the linked list and is connected to another node
                    //Set the second node to be the first on the grid and remove the second node's reference to this particle node
                    nextNode->prev = NULL;
                    grid[oldLocationX][oldLocationY] = nextNode;
                } else if (nextNode == NULL && prevNode != NULL) { //Particle is the last in the linked list
                    //Remove the neighbor's connection to this particle node
                    prevNode->next = NULL;
                } else if (nextNode == NULL && prevNode == NULL) { //Grid is empty except for this particle node
                    //Set the grid at this location to NULL cause the particle node is being removed and it is the only node there
                    grid[oldLocationX][oldLocationY] =  NULL;
                }
                //Particle node was removed so remove references to neighbors
                particleNodeI->next = NULL;
                particleNodeI->prev = NULL;

                //Put the node back in the next grid
                ParticleNode* node = grid[newLocationX][newLocationY];
                if (node == NULL) { //Grid is empty so just set the node in the grid
                    particleNodeI->gridX = newLocationX;
                    particleNodeI->gridY = newLocationY;
                    grid[newLocationX][newLocationY] = particleNodeI;
                } else { //Add the particle node to the grid but place it ordered in the linked list based on index to make memory lookup slightly faster
                    particleNodeI->gridX = newLocationX;
                    particleNodeI->gridY = newLocationY;
                    ParticleNode* firstNode = node;
                    ParticleNode* pastNode = NULL;
                    while (node->index < i) {
                        pastNode = node;
                        node = node->next;
                        if (node == NULL) break;
                    }
                    //Node should equal either null or a node we need to place our particle node in front of
                    if (node == firstNode) {
                        //Just push back the first node
                        particleNodeI->next = node;
                        node->prev = particleNodeI;
                        grid[newLocationX][newLocationY] = particleNodeI;
                    } else if (node != NULL) {
                        //Put the particle node between two nodes
                        ParticleNode* leftNode = node->prev;
                        ParticleNode* rightNode = node;
                        particleNodeI->prev = leftNode;
                        leftNode->next = particleNodeI;
                        particleNodeI->next = rightNode;
                        rightNode->prev = particleNodeI;
                    } else {
                        //Put the particle node at the end
                        ParticleNode* leftNode = pastNode;
                        particleNodeI->prev = leftNode;
                        leftNode->next = particleNodeI;
                    }
                }
            }

            //Reset the particle's forces because we just moved the particle and we don't need the forces anymore
            //particleI->ax = particleI->ay = 0;
        }
        for( int i = 0; i < n; i++ )
        {
            //Get the particle Node and particle to simulate
            ParticleNode* particleNodeI = particleNodes[i];
            particle_t *particleI = particleNodeI->particle;
            //Reset the particle's forces because we just moved the particle and we don't need the forces anymore
            particleI->ax = particleI->ay = 0;
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {
            //
            // Computing statistical data
            //
            if (navg) {
                absavg +=  davg/navg;
                nabsavg++;
            }
            if (dmin < absmin) absmin = dmin;

            //
            //  save if necessary
            //
            if( fsave && (step%SAVEFREQ) == 0 )
                save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

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
        fprintf(fsum,"%d %g\n",n,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );
    for (int i=0; i<gridSize; i++) {
        free(grid[i]);
    }
    free(grid);
    for (int i = 0; i < n; i++) {
        free(particleNodes[i]);
    }
    free(particleNodes);
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
