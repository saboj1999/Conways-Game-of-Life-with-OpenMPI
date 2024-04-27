# Implementing OpenMPI in Conway's Game of Life
- https://github.com/saboj1999/Conways-Game-of-Life-with-OpenMPI.git

1. Compile 'hw3.c' with the following command: "mpicc hw3.c -o hw3.o"
2. Run the 'test' script inside of the 'Tests' folder with: "./test"
- This will compile and run the sequential and MPI versions of the script and test to make sure they have identical outputs for a series of parameters.
3. Run the time performance script with: "./runTimePerformance"
- This will run the command "mpiexec --oversubscribe -n x ./hw3.o 1024 1024 1000" with x as 1,2,4,8,16,32,64,128 and 256 processes. It will generate average times for 5 runs of each process-count and populate the file 'time_average_outputs.txt'.
4. Run the graphing script for time performance with: "python3 graphTimePerformance.py"
- This will generate 2 graphs for the MPI implemenation with 1024x1024 grid, 1000 generations. One for efficiency and one for speed-up.