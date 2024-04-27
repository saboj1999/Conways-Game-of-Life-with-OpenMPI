/*

Name: John Sabo
Id: 916454209
Homework #: 3

To Compile:
mpicc blocking_point_mpi.c -o blocking_point_mpi.o

To Run:
mpiexec -n 2 ./blocking_point_mpi.o <X-Dimension> <Y-Dimension> <Generations>

*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>
#include <time.h>

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
 * Purpose: Free memory used for a 2D array of dimension rowsxcols
 * Parameters: 
 *    array  -> The array
 * Returns: none
 */
void free2DArray(double **array) {
    free(array[0]);  // Free the continuous block of data
    free(array);     // Free the row pointers
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
 * Purpose: Print contents of the array
 * Parameters: 
 *     a    -> Reference to array
 *   mrows  -> numbers of rows in the 2D array
 *   ncols  -> numbers of columns in the 2D array
 * Returns: none
 */
void printArray(double **array, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.1f ", array[i][j]);
        }
        printf("\n");
    }
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

/*
 * Purpose: Distribute top and bottom rows to the above and below processes
 *        ~ to emulate the 'ghost' cells that would be there.
 * Parameters: 
 *   grid   -> reference to array of current generation
 *  numRows -> N dimension of 2D NxM array
 *  numCols -> M dimension of 2D NxM array
 *   rank   -> Rank of this process
 *   size   -> Total number of MPI processes
 * Returns: none
 */
void sendReceiveGhostRows(double **grid, int numRows, int numCols, int rank, int size) {
    MPI_Status status;
    int up, down;

    if (rank == 0) {  
        up = size - 1;  // No process above, use bottom process
        down = rank + 1;
    } else if (rank == size - 1) {  
        up = rank - 1;
        down = 0;  // No process below, use top process
    } else { 
        up = rank - 1;
        down = rank + 1;
    }

    // Send to up, receive from down
    MPI_Sendrecv(grid[1], numCols, MPI_DOUBLE, up, 0,
                 grid[numRows - 1], numCols, MPI_DOUBLE, down, 0,
                 MPI_COMM_WORLD, &status);

    // Send to down, receive from up
    MPI_Sendrecv(grid[numRows - 2], numCols, MPI_DOUBLE, down, 1,
                 grid[0], numCols, MPI_DOUBLE, up, 1,
                 MPI_COMM_WORLD, &status);

    for (int j = 0; j < numCols; ++j) {
        grid[0][j] = grid[numRows - 2][j];  // Fill left ghost column
        grid[numRows - 1][j] = grid[1][j];  // Fill right ghost column
    }
}


int main(int argc, char **argv) {
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int N, M, maxGen, localN;
  double **localGrid = NULL, **nextGrid = NULL;
  double *localGridFlat = NULL;
  double startTime, endTime;

  if (argc != 4) {
      if (rank == 0) {
          printf("Usage: %s <X-Dimension> <Y-Dimension> <Generations>\n", argv[0]);
      }
      MPI_Finalize();
      return -1;
  }

  N = atoi(argv[1]);
  M = atoi(argv[2]);
  maxGen = atoi(argv[3]);

  if (N % size != 0) {
      if (rank == 0) {
          printf("The number of rows must be divisible by the number of processes.\n");
      }
      MPI_Finalize();
      return -1;
  }

  localN = N / size;
  localGridFlat = malloc(localN * M * sizeof(double));
  localGrid = allocArray(localN, M);
  nextGrid = allocArray(localN, M);

  if (rank == 0) 
  {
      double **fullGrid = NULL;
      fullGrid = allocArray(N, M);

      srand48(123456); // Seed for random number generation
      fullGrid = initRandomArray(fullGrid, N, M);

      ghostFill(&fullGrid, N, M);
      // printf("Initial full grid:\n");
      printArray(fullGrid, N, M);

      startTime = MPI_Wtime();
      MPI_Scatter(&(fullGrid[0][0]), localN * M, MPI_DOUBLE, localGridFlat, localN * M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      free2DArray(fullGrid);
  } 
  else 
  {
      MPI_Scatter(NULL, localN * M, MPI_DOUBLE, localGridFlat, localN * M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  // Transform 1D localGridFlat into 2D localGrid array for each process
  for (int i = 0; i < localN; i++) 
  {
      for (int j = 0; j < M; j++) 
      {
          localGrid[i][j] = localGridFlat[i * M + j];
      }
  }

  // Run Game of Life
  bool changed = true;
  int gen = 0;
  while(gen < maxGen && changed) 
  {
    gen++;
    sendReceiveGhostRows(localGrid, localN, M, rank, size);
    changed = step(localGrid, &nextGrid, localN, M);
      if (changed) 
      {
        double **temp = localGrid;
        localGrid = nextGrid;
        nextGrid = temp;
      }
  }

  // Transform 2D localGrid into 1D localGridFlat array for gather method
  for (int i = 1; i < localN - 1; i++) 
  {
      for (int j = 1; j < M - 1; j++) 
      {
          localGridFlat[(i - 1) * M + (j - 1)] = localGrid[i][j];
      }
  }

  if (rank == 0) 
  {
      double **finalGrid = NULL;
      finalGrid = allocArray(N, M);

      MPI_Gather(localGridFlat, localN * M, MPI_DOUBLE, &(finalGrid[0][0]), localN * M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      endTime = MPI_Wtime();

      // printf("Final full grid:\n");
      // printArray(finalGrid, N, M);
      // printf("Time taken = %lf seconds\n", endTime - startTime);

      free2DArray(finalGrid);
  } 
  else 
  {
      MPI_Gather(localGridFlat, localN * M, MPI_DOUBLE, NULL, localN * M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  free2DArray(localGrid);
  free2DArray(nextGrid);
  free(localGridFlat);

  MPI_Finalize();
  return 0;
}
