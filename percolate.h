/*
 *  Main header file for percolation code.
 */

/*
 * Number of dimension
 */

#define ndims 2

/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */

void mapwrite(char *percfile, int **map, int w, int h, int ncluster);

/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);
