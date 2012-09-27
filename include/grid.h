#ifndef GRID
#define GRID

#include <mpi.h>
#include "input.h"
#include "mpiinterface.h"

class cgrid
{
  public:
    // functions
    cgrid(cmpi *);
    ~cgrid();
    int readinifile(cinput *);
    int init();
    int create(cinput *);
    int calculate();
    int save();
    int load();

    // variables
    int itot;
    int jtot;
    int ktot;

    int imax;
    int jmax;
    int kmax;

    int iblock;
    int jblock;
    int kblock;

    int igc;
    int jgc;
    int kgc;

    int icells;
    int jcells;
    int kcells;
    int ncells;
    int istart;
    int jstart;
    int kstart;
    int iend;
    int jend;
    int kend;

    double xsize;
    double ysize;
    double zsize;

    double dx;
    double dy;
    double *dz;
    double *dzh;
    double *dzi;
    double *dzhi;
    double *dzi4;
    double *dzhi4;
    
    double *x;
    double *y;
    double *z;
    double *xh;
    double *yh;
    double *zh;

    // MPI functions
    int initmpi();
    int boundary_cyclic(double *);
    int transposezx(double *, double *);
    int transposexz(double *, double *);
    int transposexy(double *, double *);
    int transposeyx(double *, double *);
    int transposeyz(double *, double *);
    int transposezy(double *, double *);

    int getmax (double *);
    int getsum (double *);
    int getprof(double *);

    int savefield3d(double *, double *, double *, char *);
    int loadfield3d(double *, double *, double *, char *);

  private:
    cmpi *mpi;
    bool allocated;

    // MPI Datatypes
    MPI_Datatype eastwestedge;
    MPI_Datatype northsouthedge;
    MPI_Datatype transposez;
    MPI_Datatype transposez2;
    MPI_Datatype transposex;
    MPI_Datatype transposex2;
    MPI_Datatype transposey;
    MPI_Datatype transposey2;
    MPI_Datatype subi;
    MPI_Datatype subj;
    MPI_Datatype subarray;
    bool mpitypes;

    inline double interp4      (const double, const double, const double, const double);
    inline double interp4biasbot(const double, const double, const double, const double);
    inline double interp4biastop(const double, const double, const double, const double);
    inline double grad4x       (const double, const double, const double, const double);
    inline double grad4xbiasbot(const double, const double, const double, const double);
    inline double grad4xbiastop(const double, const double, const double, const double);
};
#endif
