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


int main(int argc, char **argv) {
  
  // default value updated by command line argument
  int n = 5;
  int k = 1;
  float tau = .1;

  // for checking if n, k, tau are sensible
  int return_value;

  // the 2D grids of integers
  int **pop, **npop;

  // loop variables
  int i, j, new, t, ninfected;
  
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

  pop[1][2] = 1; // set first patient to infected (probably change to random nums later?)
  ninfected = 1;

  t = 0;

  while (ninfected > 0) {
    t = t + 1;
    
    for (i = 0; i < n; i ++) {

      for (j = 0; j < n; j++) {
        new = pop[i][j];

        if (new > 0) {
          new = new + 1;

          if (new > k) {
            new = -1;
          }

        }

        else {
          if (new == 0) {
            new = infect(pop, i, j, n, n, tau);
            if (new == 1) {
              ninfected++;
            }
          }
        }

        npop[i][j] = new;
      }
    }

    pop = npop;
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
}