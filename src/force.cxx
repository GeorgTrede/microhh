#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "force.h"
#include "defines.h"

cforce::cforce(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  std::printf("Creating instance of object force\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cforce::~cforce()
{
  std::printf("Destroying instance of object force\n");
}

int cforce::exec(double dt)
{
  flux((*fields->ut).data, (*fields->u).data, grid->dz, dt);

  return 0;
}

int cforce::flux(double * restrict ut, double * restrict u, double * restrict dz, double dt)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  
  double uavg, utavg;

  uavg  = 0.;
  utavg = 0.;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        uavg  = uavg  + u [ijk]*dz[k];
        utavg = utavg + ut[ijk]*dz[k];
      }

  mpi->getsum(&uavg);
  mpi->getsum(&utavg);

  uavg  = uavg  / (grid->itot*grid->jtot*grid->zsize);
  utavg = utavg / (grid->itot*grid->jtot*grid->zsize);

  const double uflux = 0.0282;

  double fbody; 
  fbody = (uflux - uavg) / dt - utavg;

  for(int n=0; n<grid->ncells; n++)
    ut[n] += fbody;

  return 0;
}

