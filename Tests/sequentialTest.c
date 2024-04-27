/*

Name: John Sabo
Id: 916454209
Homework #: 2

To Compile:
gcc -fopenmp hw2.c -o hw2.o

To Run:
./hw2.o <X-Dimension> <Y-Dimension> <Generations> <Thread-Count>

*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdbool.h>
#include <omp.h>
#include <unistd.h>

/*
 * Purpose: Return current time in unit of seconds
 * Parameters: none
 * Returns: a double of the current time in seconds
 */
double getTime(void) {
  struct timeval tval;

  gettimeofday(&tval, NULL);

  return( (double)tval.tv_sec + (double)tval.tv_usec/1000000.0 );
}

/*
 * Purpose: Allocate memory for a 2D array of dimension PxQ
 * Parameters: 
 *      P   -> first dimension of array
 *      Q   -> second dimension of array
 * Returns: The array
 */
double **allocArray(int P, int Q) {
  int i;
  double *p, **a;
  
  p = (double *)malloc(P*Q*sizeof(double));
  a = (double **)malloc(P*sizeof(double*));

  if (p == NULL || a == NULL) 
    printf("Error allocating memory\n");

  /* for row major storage */
  for (i = 0; i < P; i++)
    a[i] = &p[i*Q];
  
  return a;
}

/*
 * Purpose: Initialize entire array with specific value
 * Parameters: 
 *     a    -> Reference to array
 *   mrows  -> numbers of rows in the 2D array
 *   ncols  -> numbers of columns in the 2D array
 *   value  -> value that will be assigned to entire array
 * Returns: The array
 */
double **initArray(double **a, int mrows, int ncols, double value) {
  int i,j;

  for (i=0; i<mrows; i++)
    for (j=0; j<ncols; j++)
      // a[i][j] = drand48()*value;
      a[i][j] = value;
  
  return a;
}

/*
 * Purpose: Initialize array with same values every time
 *          ~ Will be used for comparing results with and without OpenMP
 * Parameters: 
 *     a    -> Reference to array
 *   mrows  -> numbers of rows in the 2D array
 *   ncols  -> numbers of columns in the 2D array
 * Returns: The array
 */
double **initTestArray(double **a, int mrows, int ncols) {
  int i,j;

  for (i=0; i<mrows; i++)
    for (j=0; j<ncols; j++)
      if(i % 2 == 0 && j % 2 == 0)
      {
        a[i][j] = 0.0;
      }
      else
      {
        a[i][j] = 1.0;
      }

  return a;
}

/*
 * Purpose: Print contents of the array
 * Parameters: 
 *     a    -> Reference to array
 *   mrows  -> numbers of rows in the 2D array
 *   ncols  -> numbers of columns in the 2D array
 * Returns: none
 */
void printArray(double **a, int mrows, int ncols) {
  int i,j;
  
  for (i=1; i<mrows-1; i++) {
    for (j=1; j<ncols-1; j++)
      printf("%.1f ", a[i][j]);
    printf("\n");
  }
}

/*
 * Purpose: Initalize array randomly with values 1.0 or 0.0 for each index
 * Parameters: 
 *     a    -> Reference to array
 *   mrows  -> numbers of rows in the 2D array
 *   ncols  -> numbers of columns in the 2D array
 * Returns: The array
 */
double **initRandomArray(double **a, int mrows, int ncols) {
  int i,j;

  for (i=0; i<mrows; i++)
    for (j=0; j<ncols; j++)
        a[i][j] = drand48() > 0.5 ? 1.0 : 0.0;
  
  return a;
}

/*
 * Purpose: Fills the border of the 2D array with values from opposite side of array
 *          ~ creates a border that mimics values to make fetching neighbors easier
 * Parameters: 
 *      a   -> Reference to array
 *      N   -> N dimension of 2D NxM array
 *      M   -> M dimension of 2D NxM array
 * Returns: none
 */
void ghostFill(double ***a, int N, int M)
{
    int i, j;

    for (i = 1; i < N - 1; i++) 
    {
        (*a)[i][M - 1] = (*a)[i][1];
        (*a)[i][0] = (*a)[i][M - 2];
    }

    for (j = 0; j < M; j++) 
    {
        (*a)[N - 1][j] = (*a)[1][j];
        (*a)[0][j] = (*a)[N - 2][j];
    }
}

/*
 * Purpose: Find sum of 8 neighors around cell, 
 *          ~ determine whether cell will be alive or dead next generation,
 *          ~ write to previous generation array with cell's new status,
 *          ~ effectively making (previous -> current) and (current -> previous)
 * Parameters: 
 *      i    -> row index of cell
 *      j    -> column index of cell
 *  current  -> reference to array of current generation
 *  previous -> reference to array of previous generation
 * Returns: boolean - whether this cell is alive or dead in next generation
 */
bool updateCell(int i, int j, double **current, double ***previous)
{
    double TL, TM, TR, ML, MM, MR, BL, BM, BR, total;

    TL = current[i - 1][j - 1];
    TM = current[i][j - 1];
    TR = current[i + 1][j - 1];

    ML = current[i - 1][j];
    MM = current[i][j];
    MR = current[i + 1][j];

    BL = current[i - 1][j + 1];
    BM = current[i][j + 1];
    BR = current[i + 1][j + 1];

    total = TL + TM + TR + ML + MR + BL + BM + BR;

    if(MM == 1.0)
    {
        if (total > 1.0 && total < 4.0) 
        { 
            (*previous)[i][j] = 1.0;
            return true;
        } 
        else 
        {
            (*previous)[i][j] = 0.0;
            return false;
        }
    }
    else
    {
        if (total == 3.0) 
        { 
            (*previous)[i][j] = 1.0;
            return true;
        }
        else
        {
            (*previous)[i][j] = 0.0;
            return false;
        }
    }
    return false;
}

/*
 * Purpose: Loop through all rows and columns of current generation array,
 *          ~ call $updateCell for each index of i and j
 * Parameters: 
 *  current  -> reference to array of current generation
 *  previous -> reference to array of previous generation
 *      N   -> N dimension of 2D NxM array
 *      M   -> M dimension of 2D NxM array
 * Returns: boolean - whether next generation will have changes
 */
bool step(double **current, double ***previous, int N, int M)
{
    int i, j;
    bool isAlive = false;
    bool isCell;
    
    for (i = 1; i < N - 1; i++) {
        for (j = 1; j < M - 1; j++) {
            isCell = updateCell(i, j, current, previous);

            if(isCell){
                isAlive = true;
            }
        }
    }
    return isAlive;
}


int main(int argc, char **argv) 
{
    int N, M, maxGen, threadCount;
    double **a=NULL, **b=NULL;
    double starttime, endtime;

    if (argc != 4) {
      printf("Usage: %s <X-Dimension> <Y-Dimension> <Generations>\n", argv[0]);
      exit(-1);
    }
    
    N = atoi(argv[1]);
    N += 2;
    M = atoi(argv[2]);
    M += 2;

    maxGen = atoi(argv[3]);
    maxGen -= 1;
    
    /* Allocate memory for matrices */
    a = allocArray(N, M);
    b = allocArray(N, M);
    
    /* Initialize the matrices */
    srand48(123456);
    a = initRandomArray(a, N, M);
    b = initArray(b, N, M, 0.0);
    ghostFill(&a, N, M);
    printArray(a, N, M);
    // a = initTestArray(a, N, M);

    /* Run Game of Life */
    // printf("Starting Game of Life\n");
    bool isAlive = true;
    int counter = 0;
    starttime = getTime();
    while(isAlive && counter < maxGen)
    {
        if(counter % 2 == 0)
        {
            ghostFill(&a, N, M);
            // printArray(a, N, M);
            isAlive = step(a, &b, N, M);
        }
        else
        {
            ghostFill(&b, N, M);
            // printArray(b, N, M);
            isAlive = step(b, &a, N, M);
        }
        counter++;
        // printf("Generation: %d Completed\n", counter);
    }
    endtime = getTime();

    // if(counter % 2 == 0)
    //   printArray(b, N, M);
    // else
    //   printArray(a, N, M);

    // printf("Time taken = %lf seconds\n", endtime-starttime);
    // printf("%lf\n", endtime-starttime);

    return 0;
}