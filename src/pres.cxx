#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "pres.h"
#include "defines.h"

cpres::cpres(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  // std::printf("Creating instance of object pres\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  pres_g2  = new cpres_g2 (grid, fields, mpi);
  pres_g42 = new cpres_g42(grid, fields, mpi);
  pres_g4  = new cpres_g4 (grid, fields, mpi);
}

cpres::~cpres()
{
  delete pres_g2;
  delete pres_g42;
  delete pres_g4;

  // std::printf("Destroying instance of object pres\n");
}

int cpres::readinifile(cinput *inputin)
{
  // input parameters
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&ipres, "physics", "ipres");

  // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

int cpres::init()
{
  if(ipres == 2)
    pres_g2->init();
  else if(ipres == 42)
    pres_g42->init();
  else if(ipres == 4)
    pres_g4->init();

  return 0;
}

int cpres::setvalues()
{
  if(ipres == 2)
    pres_g2->setvalues();
  else if(ipres == 42)
    pres_g42->setvalues();
  else if(ipres == 4)
    pres_g4->setvalues();

  return 0;
}

int cpres::exec(double dt)
{
  if(ipres == 0)
    return 0;

  if(ipres == 2)
  {
    // create the input for the pressure solver
    pres_g2->pres_in(fields->sd["p"]->data,
                     fields->u ->data, fields->v ->data, fields->w ->data,
                     fields->ut->data, fields->vt->data, fields->wt->data, 
                     grid->dzi, dt);

    // solve the system
    pres_g2->pres_solve(fields->sd["p"]->data, fields->sd["tmp1"]->data, fields->sd["tmp2"]->data, grid->dz,
                        grid->fftini, grid->fftouti, grid->fftinj, grid->fftoutj);

    // get the pressure tendencies from the pressure field
    pres_g2->pres_out(fields->ut->data, fields->vt->data, fields->wt->data, 
                      fields->sd["p"]->data, grid->dzhi);
  }
  else if(ipres == 42)
  {
    // create the input for the pressure solver
    pres_g42->pres_in(fields->sd["p"]->data,
                      fields->u ->data, fields->v ->data, fields->w ->data,
                      fields->ut->data, fields->vt->data, fields->wt->data, 
                      grid->dzi, dt);

    // solve the system
    pres_g42->pres_solve(fields->sd["p"]->data, fields->sd["tmp1"]->data, fields->sd["tmp2"]->data, grid->dz,
                         grid->fftini, grid->fftouti, grid->fftinj, grid->fftoutj);

    // get the pressure tendencies from the pressure field
    pres_g42->pres_out(fields->ut->data, fields->vt->data, fields->wt->data, 
                       fields->sd["p"]->data, grid->dzhi);
  }
  else if(ipres == 4)
  {
    // create the input for the pressure solver
    pres_g4->pres_in(fields->sd["p"]->data,
                     fields->u ->data, fields->v ->data, fields->w ->data,
                     fields->ut->data, fields->vt->data, fields->wt->data, 
                     grid->dzi4, dt);

    // solve the system
    pres_g4->pres_solve(fields->sd["p"]->data, fields->sd["tmp1"]->data, fields->sd["tmp2"]->data, grid->dz,
                        grid->fftini, grid->fftouti, grid->fftinj, grid->fftoutj);

    // get the pressure tendencies from the pressure field
    pres_g4->pres_out(fields->ut->data, fields->vt->data, fields->wt->data, 
                      fields->sd["p"]->data, grid->dzhi4);
  }

  return 0;
}

double cpres::check()
{
  double divmax = 0.;
  if(ipres == 2)
    divmax = pres_g2 ->calcdivergence((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi);
  else if(ipres == 42)
    divmax = pres_g42->calcdivergence((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi);
  else if(ipres == 4)
    divmax = pres_g4->calcdivergence((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi4);

  return divmax;
}
