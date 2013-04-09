#ifndef ADVEC_G42
#define ADVEC_G42

#include "grid.h"
#include "fields.h"
#include "advec.h"
#include "mpiinterface.h"

/**
 * Derived class for advection scheme with 2nd order interpolation in the vertical
 * and 4th order in the horizontal.
 */
class cadvec_g42 : public cadvec
{
  public:
    cadvec_g42(cgrid *, cfields *, cmpi *); ///< Constructor of the advection class.
    ~cadvec_g42();                          ///< Destructor of the advection class.

    unsigned long gettimelim(long unsigned int idt, double ifactor);
    double getcfl(double);
    int exec();

  private:
    cgrid   *grid;   ///< Pointer to grid class.
    cfields *fields; ///< Pointer to fields class.
    cmpi    *mpi;    ///< Pointer to mpi class.

    double calccfl(double *, double *, double *, double *, double);         ///< Calculate the CFL number.
    int advecu(double *, double *, double *, double *, double *);           ///< Calculate longitudinal velocity advection.
    int advecv(double *, double *, double *, double *, double *);           ///< Calculate latitudinal velocity advection.
    int advecw(double *, double *, double *, double *, double *);           ///< Calculate vertical velocity advection.
    int advecs(double *, double *, double *, double *, double *, double *); ///< Calculate scalar advection.

    inline double interp2(const double, const double);                                           ///< 2nd order interpolation.
    inline double interp4(const double, const double, const double, const double);               ///< 4th order interpolation.
    inline double grad4  (const double, const double, const double, const double, const double); ///< 4th order gradient.
};
#endif
