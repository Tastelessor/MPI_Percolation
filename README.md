# MPP_cw

Run following commands to run the program:

```bash
make
mpirun -n <number of processes> ./percolate <seed> <rho> <map length> <map height>
```

Noted: 

If you input 3 arguments, the 3rd argument will be the length of map, and the map will be a square. 

However, if you input 4 arguments, the 3rd parameter will be the map width, the 4th will be the map height. The program runs well with any map size and any number of processes:)

## Header File

- percolate.h: Define dimension to be 2 for the decomposition and prototype of IO functions , functions for generating random number
- arralloc.h: Define prototype for allocating memory to arrays

## Source File

- percolate.c: Main function
- arralloc.c: Implementation of allocating functions
- percio.c: Implementation of function mapwrite()
- unirand.c: Implementation of functions that generate random number
- percolate_prospect.c: Implementation of prospect
