#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include "common.h"
#include <sched.h>


static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}



int P = 0;
int bar = 0; // Counter of threads, faced barrier.
volatile int passed = 0; // Number of barriers, passed by all threads.

void barrier_wait()
{
    int passed_old = passed; // Should be evaluated before incrementing *bar*!
    
    if(__sync_fetch_and_add(&bar,1) == (P - 1))
    {
        // The last thread, faced barrier.
        bar = 0;
        // *bar* should be reseted strictly before updating of barriers counter.
        __sync_synchronize();
        passed++; // Mark barrier as passed.
    }
    else
    {
        // Not the last thread. Wait others.
        unsigned long long cycles = rdtsc(); //1
        
        //10^9 per second or 1000000000
        //1000000 per millisecond
        //
        while(passed == passed_old) {
            //unsigned long long cycles2 = rdtsc() - cycles;           //2
            //if (cycles >= 1000000) {
                sched_yield();
                //printf("Time is %d\n", (unsigned)cycles);
            //}
        }; /*sched_yield();*/
        
        
        
        
        
        // Need to synchronize cache with other threads, passed barrier.
        __sync_synchronize();
        
        //cycles = rdtsc() - cycles;           //2
        //printf("Time is %d\n", (unsigned)cycles);
    }
}










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
    particle_t particleCopy;
    int gridX, gridY;
    int index;
    struct particleNode *next, *prev;
    
    //List nearbyParticles;
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

ParticleNode ***grid;
ParticleNode ***gridCopy;
double blockSize;
int searchArea;
double fieldSize;
int gridSize;
ParticleNode **particleNodes;
ParticleNode **particleNodesCopy;

volatile bool* threadBools;






//
//  global variables
//
int n, n_threads,no_output=0;
particle_t *particles;
particle_t *particlesCopy;
FILE *fsave,*fsum;
pthread_barrier_t barrier;
bool *myThreadBooleanLock, *myThreadBooleanLocks2;
pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
pthread_spinlock_t** gridLocks;
pthread_spinlock_t* spinLockArray1;
double gabsmin=1.0,gabsavg=0.0;

//
//  check that pthreads routine call was successful
//
#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}


//
//  This is where the action happens
//
void *thread_routine( void *pthread_id )
{
    int navg,nabsavg=0;
    double dmin,absmin=1.0,davg,absavg=0.0;
    int thread_id = *(int*)pthread_id;
    
    int particles_per_thread = (n + n_threads - 1) / n_threads;
    int first = min(  thread_id    * particles_per_thread, n );
    int last  = min( (thread_id+1) * particles_per_thread, n );
    
    
    //
    //  simulate a number of time steps
    //
    for( int step = 0; step < NSTEPS; step++ )
    {
        //for (int i = 0; i < numthreads; i++) {
        threadBools[thread_id] = false;
        pthread_barrier_wait( &barrier );
        //}
        dmin = 1.0;
        navg = 0;
        davg = 0.0;
        
        //Loop through all particles
        for( int i = first; i < last; i++ )
        {
            //Get the particle Node and particle to simulate
            ParticleNode* particleNodeI = particleNodes[i];
            particle_t *particleI = particleNodeI->particle;
            
            int locationX = particleNodeI->gridX;
            int locationY = particleNodeI->gridY;
            
            SearchArea search = getSearchArea(locationX, locationY, searchArea, gridSize);
            
            //for each particle in nearbyParticles
            for (int x = search.startX; x <= search.endX; x++) {
                for (int y = search.startY; y <= search.endY; y++) {
                    ParticleNode* startingNode = grid[x][y];
                    while (startingNode != NULL) { //Go down the entire linked list on the grid location and calculate particle forces
                        //if (startingNode->index > particleNodeI->index) {
                            apply_force( *particleI, *(startingNode->particle),&dmin,&davg,&navg);
                        //}
                        //Set the node to the next node so we continue down the chain of nodes
                        startingNode = startingNode->next;
                    }
                }
            }
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
        
        //Lets make a staggered thread descent into moving the particles
        //if (thread_id != 0) {
        //    while (threadBools[thread_id - 1] == false) {
        //    sched_yield();
                //asm ("");
        //    }
        //}
        pthread_barrier_wait( &barrier );
        for( int i = first; i < last; i++ ) {
            particles[i] = particlesCopy[i];
        }
        
/*
        for (int i = first; i < last; i++) {
            //Get the particle Node and particle to simulate
            ParticleNode* particleNodeI = particleNodes[i];
            particle_t *particleI = particleNodeI->particle;
            
            move( *particleI );
            //Reset the particle's forces because we just moved the particle and we don't need the forces anymore
            particleI->ax = particleI->ay = 0;
            
        }*/
        //threadBools[thread_id] = true;
        
        //pthread_barrier_wait( &barrier );
        /*
        for (int i = first; i < last; i++) {
            //Get the particle Node and particle to simulate
            ParticleNode* particleNodeI = particleNodes[i];
            particle_t *particleI = particleNodeI->particle;
            
            move( *particleI );
            //Reset the particle's forces because we just moved the particle and we don't need the forces anymore
            particleI->ax = particleI->ay = 0;
            
            //Create a grid search area for particle collision
            int oldLocationX = particleNodeI->gridX;
            int oldLocationY = particleNodeI->gridY;
            
            //Find the old and new location in the grid of the particle
            int newLocationX = floor(particleI->x / blockSize);
            int newLocationY = floor(particleI->y / blockSize);
            
            
            
            //Update the grid cause the particle moved
            if (oldLocationX != newLocationX || oldLocationY != newLocationY) {
                //Lock the two sections of grid we are modifying
                pthread_spin_lock(&(gridLocks[oldLocationX][oldLocationY]));
                
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
                pthread_spin_unlock(&(gridLocks[oldLocationX][oldLocationY]));
                
                pthread_spin_lock(&(gridLocks[newLocationX][newLocationY]));
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
                pthread_spin_unlock(&(gridLocks[newLocationX][newLocationY]));
            }
        }
        */
        
        pthread_barrier_wait( &barrier );
        if (thread_id == 0) {
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
        pthread_barrier_wait( &barrier );
        
        //barrier_wait();
        
        if( no_output == 0 )
        {
            //
            // Computing statistical data
            //
            if (navg) {
                absavg +=  davg/navg;
                nabsavg++;
            }
            if (dmin < absmin) absmin = dmin;
        }
        
        //
        //  save if necessary
        //
        if (no_output == 0)
            if( thread_id == 0 && fsave && (step%SAVEFREQ) == 0 )
                save( fsave, n, particles );
        
    }
    
    if (no_output == 0 )
    {
        absavg /= nabsavg;
        //printf("Thread %d has absmin = %lf and absavg = %lf\n",thread_id,absmin,absavg);
        pthread_mutex_lock(&mutex);
        gabsavg += absavg;
        if (absmin < gabsmin) gabsmin = absmin;
        pthread_mutex_unlock(&mutex);
    }
    
    return NULL;
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    //
    //  process command line
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-p <int> to set the number of threads\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    n = read_int( argc, argv, "-n", 1000 );
    n_threads = read_int( argc, argv, "-p", 2 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    fsave = savename ? fopen( savename, "w" ) : NULL;
    fsum = sumname ? fopen ( sumname, "a" ) : NULL;
    
    if( find_option( argc, argv, "-no" ) != -1 )
        no_output = 1;
    
    //
    //  allocate resources
    //
    particles = (particle_t*) malloc( n * sizeof(particle_t) );
    particlesCopy = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    for (int i = 0; i < n; i++) {
        particlesCopy[i] = particles[i];
    }
    
    
    
    
    
    
    
    //Create the field grid
    blockSize = cutoff; //Use to be cutoff
    searchArea = 2;
    fieldSize = sqrt(density * n);
    gridSize = ceil(fieldSize / blockSize);
    //Create particle grid
    grid = (ParticleNode ***)malloc(gridSize * sizeof(ParticleNode **));
    for (int i=0; i<gridSize; i++)
        grid[i] = (ParticleNode **)malloc(gridSize * sizeof(ParticleNode*));
    gridCopy = (ParticleNode ***)malloc(gridSize * sizeof(ParticleNode **));
    for (int i=0; i<gridSize; i++)
        gridCopy[i] = (ParticleNode **)malloc(gridSize * sizeof(ParticleNode*));
    
    //Create grid locks
    gridLocks = (pthread_spinlock_t **)malloc(gridSize * sizeof(pthread_spinlock_t *));
    for (int i=0; i<gridSize; i++) {
        gridLocks[i] = (pthread_spinlock_t *)malloc(gridSize * sizeof(pthread_spinlock_t));
        for (int j = 0; j < gridSize; j++) {
            pthread_spin_init(&gridLocks[i][j], PTHREAD_PROCESS_PRIVATE);
        }
    }
    spinLockArray1 = (pthread_spinlock_t *)malloc(n_threads * sizeof(pthread_spinlock_t));
    for (int j = 0; j < n_threads; j++) {
        pthread_spin_init(&spinLockArray1[j], PTHREAD_PROCESS_PRIVATE);
    }
    
    
    //Create particle node array
    particleNodes = (ParticleNode**) malloc( n * sizeof(ParticleNode*) );
    particleNodesCopy = (ParticleNode**) malloc( n * sizeof(ParticleNode*) );
    
    P = n_threads;
    
    myThreadBooleanLock = (bool *) malloc( n_threads * sizeof(bool) );
    for (int i = 0; i < n_threads; i++) {
        myThreadBooleanLock[i] = false;
    }
    myThreadBooleanLocks2 = (bool *) malloc( n_threads * sizeof(bool) );
    for (int i = 0; i < n_threads; i++) {
        myThreadBooleanLocks2[i] = true;
    }
    
    //Clear the particle grid
    for (int x = 0; x < gridSize; x++) {
        for (int y = 0; y < gridSize; y++) {
            grid[x][y] = NULL;
        }
    }
    
    threadBools = (bool*) malloc(n_threads * sizeof(bool) );;
    
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
    
    
    
    
    
    
    
    
    
    
    pthread_attr_t attr;
    P( pthread_attr_init( &attr ) );
    P( pthread_barrier_init( &barrier, NULL, n_threads ) );
    
    
    int *thread_ids = (int *) malloc( n_threads * sizeof( int ) );
    for( int i = 0; i < n_threads; i++ )
        thread_ids[i] = i;
    
    pthread_t *threads = (pthread_t *) malloc( n_threads * sizeof( pthread_t ) );
    
    //
    //  do the parallel work
    //
    double simulation_time = read_timer( );
    for( int i = 1; i < n_threads; i++ )
        P( pthread_create( &threads[i], &attr, thread_routine, &thread_ids[i] ) );
    
    thread_routine( &thread_ids[0] );
    
    for( int i = 1; i < n_threads; i++ )
        P( pthread_join( threads[i], NULL ) );
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);
    
    if( find_option( argc, argv, "-no" ) == -1 )
    {
        gabsavg /= (n_threads*1.0);
        //
        //  -the minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf( ", absmin = %lf, absavg = %lf", gabsmin, gabsavg);
        if (gabsmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting ");
        if (gabsavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting ");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,n_threads,simulation_time);
    
    //
    //  release resources
    //
    P( pthread_barrier_destroy( &barrier ) );
    free(myThreadBooleanLock);
    free(myThreadBooleanLocks2);
    P( pthread_attr_destroy( &attr ) );
    free( thread_ids );
    free( threads );
    for (int i=0; i<gridSize; i++) {
        free(grid[i]);
    }
    free(grid);
    for (int i=0; i<gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            pthread_spin_destroy(&gridLocks[i][j]);
        }
        //    free(gridLocks[i]);
    }
    for (int j = 0; j < n_threads; j++) {
        pthread_spin_destroy(&spinLockArray1[j]);
    }
    free(gridLocks);
    for (int i = 0; i < n; i++) {
        free(particleNodes[i]);
    }
    free(particleNodes);
    free( particles );
    if( fsave )
        fclose( fsave );
    if( fsum )
        fclose ( fsum );
    
    return 0;
}

