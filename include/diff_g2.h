#ifndef DIFF_G2
#define DIFF_G2

#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "mpiinterface.h"

class cdiff_g2 : public cdiff
{
  public:
    cdiff_g2(cgrid *, cfields *, cmpi *);
    ~cdiff_g2();

    int setvalues();
    int exec();

    unsigned long gettimelim(unsigned long, double);
    double getdn(double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double dnmul;

    int diffc(double *, double *, double *, double *, double);
    int diffw(double *, double *, double *, double *, double);
};
#endif
