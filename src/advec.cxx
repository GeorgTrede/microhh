#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "advec.h"
#include "defines.h"

cadvec::cadvec(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  advec_g2   = new cadvec_g2  (grid, fields, mpi);
  advec_g2i4 = new cadvec_g2i4(grid, fields, mpi);
  advec_g42  = new cadvec_g42 (grid, fields, mpi);
  advec_g4   = new cadvec_g4  (grid, fields, mpi);
  advec_g4m  = new cadvec_g4m (grid, fields, mpi);
}

cadvec::~cadvec()
{
  delete advec_g2;
  delete advec_g2i4;
  delete advec_g42;
  delete advec_g4;
  delete advec_g4m;
}

int cadvec::readinifile(cinput *inputin)
{
  // input parameters
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&iadvec, "physics", "iadvec");

  // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

double cadvec::getcfl(double dt)
{
  double cfl;

  // in case of no advection set cfl to a small number to avoid zero divisions
  if(iadvec == 0)
    cfl = dsmall;
  
  if(iadvec == 2)
    cfl = advec_g2->calccfl((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);
  else if(iadvec == 24)
    cfl = advec_g2i4->calccfl((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);
  else if(iadvec == 42)
    cfl = advec_g42->calccfl((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);
  else if(iadvec == 4)
    cfl = advec_g4->calccfl((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);
  else if(iadvec == 44)
    cfl = advec_g4m->calccfl((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);

  return cfl;
}

int cadvec::exec()
{
  if(iadvec == 0)
    return 0;

  if(iadvec == 2)
  {
    advec_g2->advecu((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
    advec_g2->advecv((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
    advec_g2->advecw((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzhi);

    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      advec_g2->advecs((*it->second).data, (*fields->s[it->first]).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi);
  }
  else if(iadvec == 24)
  {
    advec_g2i4->advecu((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
    advec_g2i4->advecv((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
    advec_g2i4->advecw((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzhi);

    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      advec_g2i4->advecs((*it->second).data, (*fields->s[it->first]).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi);
  }
  else if(iadvec == 42)
  {
    advec_g42->advecu((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
    advec_g42->advecv((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
    advec_g42->advecw((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzhi);

    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      advec_g42->advecs((*it->second).data, (*fields->s[it->first]).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi);
  }
  else if(iadvec == 4)
  {
    advec_g4->advecu((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi4 );
    advec_g4->advecv((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi4 );
    advec_g4->advecw((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzhi4);

    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      advec_g4->advecs((*it->second).data, (*fields->s[it->first]).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi4);
  }
  else if(iadvec == 44)
  {
    advec_g4m->advecu((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi4 );
    advec_g4m->advecv((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi4 );
    advec_g4m->advecw((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzhi4);

    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      advec_g4m->advecs((*it->second).data, (*fields->s[it->first]).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi4);
  }

  return 0;
}

