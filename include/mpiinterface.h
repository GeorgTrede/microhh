#ifndef MPIINTERFACE
#define MPIINTERFACE

#include <mpi.h>
#include "input.h"
#include "grid.h"

class cmpi
{
  public:
    cmpi(cgrid *);
    ~cmpi();

    int readinifile(cinput *);
    int init();

    int boundary_cyclic(double *);
    int transposezx(double *, double *);
    int transposexz(double *, double *);
    int transposexy(double *, double *);
    int transposeyx(double *, double *);
    int transposeyz(double *, double *);

    int getmax(double *);
    int getsum(double *);

    int writefield3d(double *, char *);
    int readfield3d(double *, char *);

    int nprocs;
    int npx;
    int npy;
    int mpiid;
    int mpicoordx;
    int mpicoordy;

    int nnorth;
    int nsouth;
    int neast;
    int nwest;

  private:
    cgrid *grid;

    bool initialized;
    bool allocated;

    MPI_Comm commxy;
    MPI_Comm commx;
    MPI_Comm commy;

    MPI_Datatype eastwestedge;
    MPI_Datatype northsouthedge;
    MPI_Datatype transposez;
    MPI_Datatype transposex;
    MPI_Datatype transposex2;
    MPI_Datatype transposey;
    MPI_Datatype transposezgc;
    MPI_Datatype subarray;

    MPI_Request *reqs;
};
#endif
