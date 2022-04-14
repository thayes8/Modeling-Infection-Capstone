/* 
  HODGE-C -- A C implementation of Martin Gerhard & Heike Schuster's hodge-podge machine.

  This is a version drastically cut down by John Burkardt.
  The only other file it needs is "hodge.map", which is used to pick the colors.

  Parameters for the model:

    n:    dimenstion of the grid of hospital beds
    k:    number of days of contagion
    tau:  transmission rate for the infected
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <iostream>
#include <fstream>

#include <trng/mt19937_64.hpp>
#include <trng/mt19937.hpp>
#include <trng/lcg64_shift.hpp>
#include <trng/normal_dist.hpp>
#include <trng/uniform_dist.hpp>
#include <trng/exponential_dist.hpp>

/********************************************
 * Need at least this many rows and columns *
 ********************************************/
const int MINIMUM_DIMENSION = 1;
const int MINIMUM_K = 0;
const int MINIMUM_TAU = 0;

void getArguments(int argc, char *argv[], int *n, int *k, float *tau);
int assert_minimum_value(char which_value[16], int actual_value, int expected_value);
void pluralize_value_if_needed(int value); 
void printGrid(int **our_pop, int n);
void initialize_grid(int **grid, int n);
void spread_infection(int **pop, int **npop, int n, int k, float tau);
bool infect(int **pop, int i, int j, float tau, int nRows, int nCols);

using namespace std;
int main(int argc, char **argv) {

  trng::mt19937_64 RNengine1;
  trng::uniform_dist<> uni(0, 1);

  float rand = uni(RNengine1);
  float rand2 = uni(RNengine1);

  printf("rand = %f, rand2 = %f\n", rand, rand2);

  
  // default value updated by command line argument
  int n = 5;
  int k = 5;
  float tau = .1;

  // for checking if n, k, tau are sensible
  int return_value;

  // the 2D grids of integers
  int **pop, **npop;

  // loop variables
  int row;

  // set random engine
  
  // get command line arguments
  getArguments(argc, argv, &n, &k, &tau);

  // make sure dimension is not <= 0
  return_value = assert_minimum_value("n", n, MINIMUM_DIMENSION);
  // make sure number of days of contagion is not < 1
  return_value += assert_minimum_value("k", k, MINIMUM_K);
  // make sure tau is >= 0
  return_value += assert_minimum_value("tau", tau, MINIMUM_TAU);

  if (return_value != 0) {
    exit(-1);
  }
  
  pop = (int**)malloc(n * n * sizeof(int));
  npop = (int**)malloc(n * n * sizeof(int));
  for (row = 0; row < n; row ++) {
    pop[row] = (int*)malloc(n * sizeof(int));
    npop[row] = (int*)malloc(n * sizeof(int));
  }

  initialize_grid(pop, n);
  initialize_grid(npop, n);
  spread_infection(pop, npop, n, k, tau);

  free(pop);
  free(npop);
}

void spread_infection(int **pop, int **npop, int n, int k, float tau) {
  int t, i, j, new_value;
  
  int ninfected = 1;

  pop[1][2] = 1; // set first patient to infected (probably change to random nums later?

  t = 0;

  int a = 0;
  printGrid(pop, n);
  while (a++ < 3) {
    t = t + 1;

    for (i = 0; i < n; i ++) {

      for (j = 0; j < n; j++) {
        
        new_value = pop[i][j];
        if (new_value > 0) {
          new_value = new_value + 1;

          if (new_value > k) {
            new_value = -1;
            ninfected--;
          }

        }

        else {
          if (new_value == 0) {
            new_value = infect(pop, i, j, tau, n, n);
            ninfected;
          }
        }

        npop[i][j] = new_value;

      }


    }

    for (i = 0; i < n; i ++) {

      for (j = 0; j < n; j++) {

        pop[i][j] = npop[i][j];

      }
    }
    printGrid(pop, n);
  } 

}

/**
Looks at all of a current cell's neighbors and gets a random num for each
infected neighbor and generates a random num from 0 to 1, then compares that
num to the infection rate to determine if it gets infected.
Variable Definitions:
i = current x pos
j = current y pos
nRows = numRows
nCols = numColumns
tau = Transmission Rate
**/
bool infect(int **pop, int i, int j, float tau, int nRows, int nCols){
    trng::mt19937_64 RNengine1;
    trng::uniform_dist<> uni(0, 1);

    float rand = uni(RNengine1);
    printf("rand = %f", rand);

    //Tracks whether current cell has been infected
    int t = 0;

    //if not the leftmost wall
    if(i > 0) {
        //if left neighbor is sick
        if(pop[i-1][j] > 0){
            t = (rand < tau);
        }
    }
    //if i is not the rightmost wall
    if(i < nRows - 1) {
        //if left neighbor is sick
        if(pop[i+1][j] > 0){
            t = t + (rand < tau);
        }
    }

    if(j > 0) {
        //if left neighbor is sick
        if(pop[i][j-1] > 0){
            t = t + (rand < tau);
        }
    }

    if(j < nCols - 1) {
      
        if(pop[i][j+1] > 0){
            t = t + (rand < tau);
        }
    }

    bool p = 0;
    
    if(t > 0){
        p = 1;
    }
    
    return p;

}

void initialize_grid(int **grid, int n) {
  int row;
  int column;

  for (row = 0; row < n; row ++) {
    for (column = 0; column < n; column++) {
      grid[row][column] = 0;
    }
  }
}


void getArguments(int argc, char *argv[], int *n, int *k, float *tau) {
  char *nvalue, *kvalue, *tvalue;
  int c;    // result from getopt calls

  while ((c = getopt(argc, argv, "n:k:t:")) != -1) {
    switch (c) {
      case 'n':
        nvalue = optarg;
        *n = atoi(nvalue);
        break;
      case 'k':
        kvalue = optarg;
        *k = atoi(kvalue);
        break;
      case 't':
        tvalue = optarg;
        *tau = atof(tvalue);
        break;
      case '?':
      default:
        fprintf(stderr, "Usage %s [-n dimension of grid] [-k number of days in contagion] [-t transmissablitiy rate]\n", argv[0]);
        
        exit(-1);  
    }
  }
}

/*******************************************************************************
 * Make sure a value is >= another value, print error and return -1 if it isn't
 ******************************************************************************/
int assert_minimum_value(char which_value[16], int actual_value,
        int expected_value)
{
    int retval;

    if(actual_value < expected_value)
    {
        fprintf(stderr, "ERROR: %d %s", actual_value, which_value);
        pluralize_value_if_needed(actual_value);
        fprintf(stderr, "; need at least %d %s", expected_value, which_value);
        pluralize_value_if_needed(expected_value);
        fprintf(stderr, "\n");
        retval = -1;
    }
    else
        retval = 0;

    return retval;
}

/*****************************************************
 * Add an "s" to the end of a value's name if needed *
 *****************************************************/
void pluralize_value_if_needed(int value)
{
    if(value != 1)
        fprintf(stderr, "s");

    return;
}

void printGrid(int **pop, int n) {
  int current_row, current_column;
  
  for (current_row = 0; current_row < n; current_row++) {

    for (current_column = 0; current_column < n; current_column ++) {
      
      printf("%d", pop[current_row][current_column]);

    }

    printf("\n");
  }
  printf("\n");
}