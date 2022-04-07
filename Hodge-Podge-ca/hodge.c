/* 
  HODGE-C -- A C implementation of Martin Gerhard & Heike Schuster's hodge-podge machine.

  This is a version drastically cut down by John Burkardt.
  The only other file it needs is "hodge.map", which is used to pick the colors.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#define MAXCELLS	128

#define S_SANE		0
#define S_INFECTED	1
#define S_ILL		2

#define bool		int
#define TRUE		1
#define FALSE		0

typedef enum 
{
  SANE = S_SANE, INFECTED = S_INFECTED, ILL = S_ILL
} cell_s;

typedef unsigned char cell_v;

typedef struct 
{
  cell_s state;
  cell_v value;
} cell_t, CELL;

typedef cell_t grid_t[MAXCELLS][MAXCELLS];

typedef unsigned char color_t[3];

typedef struct {
	unsigned char red;
	unsigned char green;
	unsigned char blue;
} rgb_t, RGB;

typedef rgb_t cmap_t[256];


/* 
  flags 
*/
bool moore = FALSE;      /* von-Neumann default */
bool torus = FALSE;      /* bounded by default */
/* 
  NCELLS is the number of cells on each side of the grid.
*/
int ncells = 128;      
/*
  Storage for the cellular automaton.
*/
grid_t cell;
grid_t ncell;
/*
  A temporary way to get the compiler to quit complaining.
*/
color_t rgb;
/*
  The color map.
*/
cmap_t cmap;

/*
  Function prototypes.
*/
int main ( int argc, char *argv[] );
cell_s get_state ( int *val, int n );
void grid2ppm ( int t, int n, grid_t grid );
void hodge ( int g, int k1, int k2, int n );
void initialize ( int seed, int n );
void load_color_map ( char *s );
int map ( int i );
color_t *rgb_color ( int val );
int scount ( int i, int j, cell_s state );
int vsum ( int i, int j );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for the Hodge-Podge program.

  Modified:

    29 December 2007

  Parameters:

    char *CFILENAME, the name of the color map file.

    int G, the infection rate.  Dewdney says this should be between 1 and 20.

    int K1, the infected divisor.

    int K2, the ill divisor.

    int N, the number of states, from SANE to ILL.
*/
{
  char *cfilename = "hodge.map";
  int g;
  int k1;
  int k2;
  int n;

  printf ( "\n" );
  printf ( "HODGE\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  A cut-down version of Gerhardt and Schuster's \n" );
  printf ( "  hodgepodge machine.\n" );

  n = 100;
  g = 10;
  k1 = 1;
  k2 = 1;

  printf ( "\n" );
  printf ( "  The parameters for this run include:\n" );
  printf ( "\n" );
  printf ( "    the number of states N =  %d\n", n );
  printf ( "    the infection rate G =    %d\n", g );
  printf ( "    the infected divisor K1 = %d\n", k1 );
  printf ( "    the ill divisor K2 =      %d\n", k2 );

  printf ( "\n" );
  printf ( "  The color map will be read from the file %s\n", cfilename );

  load_color_map ( cfilename );

  hodge ( g, k1, k2, n );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HODGE\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

cell_s get_state ( int *val, int n )

/******************************************************************************/
/*
  Purpose:

    GET_STATE computes the state from a cell's value.

  Modified:

    29 December 2007

  Parameters:
*/
{
  if ( *val <= 0 )
  {
    *val = 0;
    return SANE;
  } 
  else if ( n <= *val )
  {
    *val = n;
    return ILL;
  } 
  else
  {
    return INFECTED;
  }

}
/******************************************************************************/

void grid2ppm ( int t, int n, grid_t grid )

/******************************************************************************/
/*
  Purpose:

    GRID2PPM writes out the grid converted to a binary color portable pixmap.

  Modified:

    29 December 2007

  Parameters:
*/
{
  static char *afilename = "hodge_%0003d.ppm";
  static int fcount = 0;
  FILE *fp;
  int i;
  int j;
  char s[255];
  unsigned char val;
/* 
  Generate the frame's file name 
*/
  sprintf ( s, afilename, fcount++ );

  fp = fopen ( s, "w" );

  if ( fp == NULL )
  {
    printf ( "\n" );
    printf ( "GRID2PPM - Fatal error!\n" );
    printf ( "  Can't open file %s\n", s );
  }
/* 
  Write the raw ppm(5) header 
*/
  fprintf ( fp, "P6\n" );
  fprintf ( fp, "# File:    %s\n", s );
  fprintf ( fp, "# Creator: Hodge (C) 1993 by Joerg Heitkoetter\n" );
  fprintf ( fp, "%d %d\n", ncells, ncells );
  fprintf ( fp, "%d\n", n );
/* 
  Turn the bitimage into RGB colors 
*/
  for ( i = 0; i < ncells; i++ )
  {
    for ( j = 0; j < ncells; j++ )
    {
      val = cell[i][j].value;

      rgb[0] = cmap[val].red;
      rgb[1] = cmap[val].green;
      rgb[2] = cmap[val].blue;

      fwrite ( rgb_color ( cell[i][j].value ), sizeof (color_t), 1, fp );
    }
  }
/* 
  Clean up 
*/
  fclose (fp);

  printf ( "\n" );
  printf ( "GRID2PPM wrote (%dx%d) file `%s' (%d colors)\n",
    ncells, ncells, s, n + 1);

}
/******************************************************************************/

void hodge ( int g, int k1, int k2, int n )

/******************************************************************************/
/*
  Purpose:

    HODGE carries out the main loop.

  Modified:

    29 December 2007

  Parameters:

    int G, the infection rate.  Dewdney says this should be between 1 and 20.

    int K1, the infected divisor.

    int K2, the ill divisor.

    int N, the number of states, from SANE to ILL.

  Local parameters:

    int NFRAMES, the number of PPM image files to write.

    int RANDOM_SEED, a seed for the random number generator.
*/
{
  int A;
  int B;
  int framestotake;
  int i;
  int j;
  int nframes = 10;
  int nill;
  int ninfected;
  int nsane;
  int nticks = 1000;
  int nval;  
  int random_seed = 12345678; 
  int S; 
  int skipbetweentakes;
  int startoffilm;
  int t = 0;

  framestotake = nframes;
  startoffilm = 100;
  skipbetweentakes = 100;

/* 
  NVAL has the range: [0..8*n + n] 
*/
  initialize ( random_seed, n );

  while ( t++ < nticks )
  {
/*
  Count the number of cells of each of the three states.
*/
    nsane = 0;
    ninfected = 0;
    nill = 0;
    for ( i = 0; i < ncells; i++)
    {
      for ( j = 0; j < ncells; j++) 
      {
        switch ( cell[i][j].state )
        {
          case SANE:
            nsane++;
            break;

          case INFECTED:
            ninfected++;
            break;

          case ILL:
            nill++;
            break;
        }
      }
    }
    printf ( "%d\t%d\t%d\t%d\t%lf\n",
      t, nsane, ninfected, nill, 
     ( double ) ninfected / ( double ) ( ncells * ncells ) );
/*
  Update the state and value of each cell.
*/
    for ( i = 0; i < ncells; i++)
    {
      for ( j = 0; j < ncells; j++) 
      {
        A = scount ( i, j, ILL );
        B = scount ( i, j, INFECTED );
        S = vsum ( i, j ) + cell[i][j].value;

        switch ( cell[i][j].state )
        {
          case SANE:
            nval = ( cell_v ) ( A / k1 + B / k2 );
            break;

          case INFECTED:
            B++;
            nval = ( cell_v ) ( S / B + g );
            break;

          case ILL:
            nval = 0;
            break;
        }
        ncell[i][j].state = get_state ( &nval, n );
        ncell[i][j].value = (cell_v) nval;
      }
    }
/* 
  Write out a PPM image file. 
*/ 
    if 
    ( 
      startoffilm <= t && 
      0 < framestotake &&
      ( t - startoffilm ) % skipbetweentakes == 0 
    )
    {
      grid2ppm ( t, n, cell );
      --framestotake;
    }
/* 
  Copy new cells to old 
*/
    memcpy ( cell, ncell, sizeof ( grid_t ) );
  }
}
/******************************************************************************/

void initialize ( int seed, int n )

/******************************************************************************/
/*
  Purpose:

    INITIALIZE initializes the cells.

  Modified:

    29 December 2007

  Parameters:
*/
{
  int i;
  int j;

  srand ( seed );

  for ( i = 0; i < ncells; i++ )
  {
   for ( j = 0; j < ncells; j++ )
   {
     cell[i][j].state = ( cell_s ) ( random () % 3 );

      if ( cell[i][j].state == SANE)
      {
        cell[i][j].value = 0;
      }
      else if ( cell[i][j].state == INFECTED)
      {
        cell[i][j].value = ( cell_v ) ( rand () % ( n - 1 ) + 1 );
      }
      else
      {
        cell[i][j].value = n;
      }

      ncell[i][j].state = SANE;
      ncell[i][j].value = 0;
    }
  }
  return;
}
/******************************************************************************/

void load_color_map ( char *cfilename )

/******************************************************************************/
/*
  Purpose:

    LOAD_COLOR_MAP loads a color map.

  Modified:

    29 December 2007

  Parameters:

    Input, char *CFILENAME, the name of the color map file.
*/
{
# define MAXLINELEN 1024

  int b;
  FILE *fp;
  int g;
  int i;
  char line[MAXLINELEN];
  int r;

  fp = fopen ( cfilename, "r");

  if ( fp == NULL )
  {
    printf ( "\n" );
    printf ( "LOAD_COLOR_MAP - Fatal error!\n" );
    printf ( "  Can't load the color map file \"%s\".\n", cfilename );
    exit ( 1 );
  }

  for ( i = 0; i < 256; i++ )
  {
    if ( fgets ( line, MAXLINELEN, fp ) == NULL )
    {
      printf ( "\n" );
      printf ( "LOAD_COLOR_MAP - Fatal error!n" );
      printf ( "  Not enough colors in color map file \"%s\"\n", cfilename );
    }

    if ( *line == '#' )
    {
      continue;
    }

    sscanf ( line, " %d %d %d", &r, &g, &b );

    cmap[i].red = r & 0xff;
    cmap[i].green = g & 0xff;
    cmap[i].blue = b & 0xff;
  }

  printf ( "\n" );
  printf ( "LOAD_COLOR_MAP:\n" );
  printf ( "  Read the color map \"%s\"\n", cfilename );

# undef MAXLINELEN
}
/******************************************************************************/

int map ( int i )

/******************************************************************************/
/*
  Purpose:

    MAP maps an index to the grid structure.

  Modified:

    29 December 2007

  Parameters:
*/
{
  if ( torus )
  {
    return ( i % ncells );
  }
  else
  {
    if ( i < 0 )
    {
      i = 0;
    }
    if ( ncells - 1 < i )
    {
      i = ncells - 1;
    }
    return ( i );
  } 

}
/******************************************************************************/

color_t *rgb_color ( int val )

/******************************************************************************/
/*
  Purpose:

    RGB_COLOR turns a value into an RGB string.

  Modified:

    29 December 2007

  Parameters:
*/
{
  rgb[0] = cmap[val].red;
  rgb[1] = cmap[val].green;
  rgb[2] = cmap[val].blue;

  return ( ( color_t *) rgb );
}
/******************************************************************************/

int scount ( int i, int j, cell_s state )

/******************************************************************************/
/*
  Purpose:

    SCOUNT counts the number of cells in a specified state.

  Modified:

    29 December 2007

  Parameters:
*/
{
  int count = 0;

/* 
  Always count von-Neumann cells 
*/
  if ( cell[map (i - 1)][j].state == state )
  {
    ++count;
  }
  if ( cell[i][map (j - 1)].state == state )
  {
    ++count;
  }
  if ( cell[i][map (j + 1)].state == state )
  {
    ++count;
  }
  if ( cell[map (i + 1)][j].state == state )
  {
    ++count;
  }
/* 
  Sometimes count Moore cells 
*/
  if ( moore )
  {
    if ( cell[map (i - 1)][map (j - 1)].state == state )
    {
      ++count;
    }
    if ( cell[map (i + 1)][map (j - 1)].state == state )
    {
      ++count;
    }
    if ( cell[map (i - 1)][map (j + 1)].state == state )
    {
      ++count;
    }
    if ( cell[map (i + 1)][map (j + 1)].state == state )
    {
      ++count;
    }
  }
  return ( count );
}
/******************************************************************************/

int vsum ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    VSUM sums the neighbor cell values.

  Modified:

    29 December 2007

  Parameters:
*/
{
  int sum = 0;
/* 
  always sum von-Neumann cells 
*/
  sum = cell[map (i - 1)][j].value
      + cell[i][map (j - 1)].value
      + cell[i][map (j + 1)].value
      + cell[map (i + 1)][j].value;

/* 
  Sometimes sum Moore cells 
*/
  if ( moore )
  {
    sum += cell[map (i - 1)][map (j - 1)].value
        + cell[map (i + 1)][map (j - 1)].value
        + cell[map (i - 1)][map (j + 1)].value
        + cell[map (i + 1)][map (j + 1)].value;
  }

  return sum;
}
