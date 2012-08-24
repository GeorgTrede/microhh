#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "advec_g4.h"
#include "defines.h"

cadvec_g4::cadvec_g4(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  std::printf("Creating instance of object advec_g4\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cadvec_g4::~cadvec_g4()
{
  std::printf("Destroying instance of object advec_g4\n");
}

double cadvec_g4::calccfl(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double dt)
{
  int    ijk;
  int    ii1,ii2,jj1,jj2,kk1,kk2;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;


  double cfl = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        cfl = std::max(cfl, std::abs(interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]))*dxi 
                          + std::abs(interp4(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]))*dyi 
                          + std::abs(interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k] );
      }

  grid->getmax(&cfl);

  cfl = cfl*dt;

  return cfl;
}

int cadvec_g4::advecu(double * restrict ut, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      ut[ijk] -= grad4(interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]) * interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]),
                       interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]),
                       interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]),
                       interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]) * interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]), dxi);

      ut[ijk] -= grad4(interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp4(u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ]),
                       interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]),
                       interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2]),
                       interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp4(u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3]), dyi);

      ut[ijk] -= grad4x(interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp4biasbot(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]),
                        interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp4       (u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]),
                        interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp4       (u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]),
                        interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp4       (u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3]))
                 * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        ut[ijk] -= grad4(interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]) * interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]),
                         interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]),
                         interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]),
                         interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]) * interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]), dxi);

        ut[ijk] -= grad4(interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp4(u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ]),
                         interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]),
                         interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2]),
                         interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp4(u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3]), dyi);

        ut[ijk] -= grad4x(interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp4(u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ]),
                          interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp4(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]),
                          interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp4(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]),
                          interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp4(u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3]))
                   * dzi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      ut[ijk] -= grad4(interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]) * interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]),
                       interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]),
                       interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]),
                       interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]) * interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]), dxi);

      ut[ijk] -= grad4(interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp4(u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ]),
                       interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]),
                       interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2]),
                       interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp4(u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3]), dyi);

      ut[ijk] -= grad4x(interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp4       (u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ]),
                        interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp4       (u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]),
                        interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp4       (u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]),
                        interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp4biastop(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]))
                 * dzi4[kend-1];
    }

  return 0;
}

int cadvec_g4::advecv(double * restrict vt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      vt[ijk] -= grad4(interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp4(v[ijk-ii3], v[ijk-ii2], v[ijk-ii1], v[ijk    ]),
                       interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]),
                       interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2]),
                       interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp4(v[ijk    ], v[ijk+ii1], v[ijk+ii2], v[ijk+ii3]), dxi);

      vt[ijk] -= grad4(interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]) * interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]),
                       interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]),
                       interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]),
                       interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]) * interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]), dyi);

      vt[ijk] -= grad4x(interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp4biasbot(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]),
                        interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp4       (v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]),
                        interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp4       (v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]),
                        interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp4       (v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3]))
                 * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        vt[ijk] -= grad4(interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp4(v[ijk-ii3], v[ijk-ii2], v[ijk-ii1], v[ijk    ]),
                         interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]),
                         interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2]),
                         interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp4(v[ijk    ], v[ijk+ii1], v[ijk+ii2], v[ijk+ii3]), dxi);

        vt[ijk] -= grad4(interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]) * interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]),
                         interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]),
                         interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]),
                         interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]) * interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]), dyi);

        vt[ijk] -= grad4x(interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp4(v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ]),
                          interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp4(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]),
                          interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp4(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]),
                          interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp4(v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3]))
                   * dzi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      vt[ijk] -= grad4(interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp4(v[ijk-ii3], v[ijk-ii2], v[ijk-ii1], v[ijk    ]),
                       interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]),
                       interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2]),
                       interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp4(v[ijk    ], v[ijk+ii1], v[ijk+ii2], v[ijk+ii3]), dxi);

      vt[ijk] -= grad4(interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]) * interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]),
                       interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]),
                       interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]),
                       interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]) * interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]), dyi);

      vt[ijk] -= grad4x(interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp4       (v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ]),
                        interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp4       (v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]),
                        interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp4       (v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]),
                        interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp4biastop(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]))
                 * dzi4[kend-1];
    }

  return 0;
}

int cadvec_g4::advecw(double * restrict wt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzhi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kstart+1)*kk1;
      wt[ijk] -= grad4(interp4(u[ijk-ii1-kk2], u[ijk-ii1-kk1], u[ijk-ii1], u[ijk-ii1+kk1]) * interp4(w[ijk-ii3], w[ijk-ii2], w[ijk-ii1], w[ijk    ]),
                       interp4(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]),
                       interp4(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2]),
                       interp4(u[ijk+ii2-kk2], u[ijk+ii2-kk1], u[ijk+ii2], u[ijk+ii2+kk1]) * interp4(w[ijk    ], w[ijk+ii1], w[ijk+ii2], w[ijk+ii3]), dxi);

      wt[ijk] -= grad4(interp4(v[ijk-jj1-kk2], v[ijk-jj1-kk1], v[ijk-jj1], v[ijk-jj1+kk1]) * interp4(w[ijk-jj3], w[ijk-jj2], w[ijk-jj1], w[ijk    ]),
                       interp4(v[ijk    -kk2], v[ijk    -kk1], v[ijk    ], v[ijk    +kk1]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]),
                       interp4(v[ijk+jj1-kk2], v[ijk+jj1-kk1], v[ijk+jj1], v[ijk+jj1+kk1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2]),
                       interp4(v[ijk+jj2-kk2], v[ijk+jj2-kk1], v[ijk+jj2], v[ijk+jj2+kk1]) * interp4(w[ijk    ], w[ijk+jj1], w[ijk+jj2], w[ijk+jj3]), dyi);

      wt[ijk] -= grad4x(interp4biasbot(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4biasbot(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]),
                        interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]),
                        interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]),
                        interp4       (w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]) * interp4       (w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]))
                 * dzhi4[kstart+1];
    }

  for(int k=grid->kstart+2; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        wt[ijk] -= grad4(interp4(u[ijk-ii1-kk2], u[ijk-ii1-kk1], u[ijk-ii1], u[ijk-ii1+kk1]) * interp4(w[ijk-ii3], w[ijk-ii2], w[ijk-ii1], w[ijk    ]),
                         interp4(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]),
                         interp4(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2]),
                         interp4(u[ijk+ii2-kk2], u[ijk+ii2-kk1], u[ijk+ii2], u[ijk+ii2+kk1]) * interp4(w[ijk    ], w[ijk+ii1], w[ijk+ii2], w[ijk+ii3]), dxi);

        wt[ijk] -= grad4(interp4(v[ijk-jj1-kk2], v[ijk-jj1-kk1], v[ijk-jj1], v[ijk-jj1+kk1]) * interp4(w[ijk-jj3], w[ijk-jj2], w[ijk-jj1], w[ijk    ]),
                         interp4(v[ijk    -kk2], v[ijk    -kk1], v[ijk    ], v[ijk    +kk1]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]),
                         interp4(v[ijk+jj1-kk2], v[ijk+jj1-kk1], v[ijk+jj1], v[ijk+jj1+kk1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2]),
                         interp4(v[ijk+jj2-kk2], v[ijk+jj2-kk1], v[ijk+jj2], v[ijk+jj2+kk1]) * interp4(w[ijk    ], w[ijk+jj1], w[ijk+jj2], w[ijk+jj3]), dyi);

        wt[ijk] -= grad4x(interp4(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]) * interp4(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]),
                          interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]),
                          interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]),
                          interp4(w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]) * interp4(w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]))
                   * dzhi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      wt[ijk] -= grad4(interp4(u[ijk-ii1-kk2], u[ijk-ii1-kk1], u[ijk-ii1], u[ijk-ii1+kk1]) * interp4(w[ijk-ii3], w[ijk-ii2], w[ijk-ii1], w[ijk    ]),
                       interp4(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]),
                       interp4(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2]),
                       interp4(u[ijk+ii2-kk2], u[ijk+ii2-kk1], u[ijk+ii2], u[ijk+ii2+kk1]) * interp4(w[ijk    ], w[ijk+ii1], w[ijk+ii2], w[ijk+ii3]), dxi);

      wt[ijk] -= grad4(interp4(v[ijk-jj1-kk2], v[ijk-jj1-kk1], v[ijk-jj1], v[ijk-jj1+kk1]) * interp4(w[ijk-jj3], w[ijk-jj2], w[ijk-jj1], w[ijk    ]),
                       interp4(v[ijk    -kk2], v[ijk    -kk1], v[ijk    ], v[ijk    +kk1]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]),
                       interp4(v[ijk+jj1-kk2], v[ijk+jj1-kk1], v[ijk+jj1], v[ijk+jj1+kk1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2]),
                       interp4(v[ijk+jj2-kk2], v[ijk+jj2-kk1], v[ijk+jj2], v[ijk+jj2+kk1]) * interp4(w[ijk    ], w[ijk+jj1], w[ijk+jj2], w[ijk+jj3]), dyi);

      wt[ijk] -= grad4x(interp4       (w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]) * interp4       (w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]),
                        interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]),
                        interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]),
                        interp4biastop(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4biastop(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]))
                 * dzhi4[kend-1];
    }

  return 0;
}

int cadvec_g4::advecs(double * restrict st, double * restrict s, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  kstart = grid->kstart;
  kend   = grid->kend;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      st[ijk] -= grad4(u[ijk-ii1] * interp4(s[ijk-ii3], s[ijk-ii2], s[ijk-ii1], s[ijk    ]),
                       u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]),
                       u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2]),
                       u[ijk+ii2] * interp4(s[ijk    ], s[ijk+ii1], s[ijk+ii2], s[ijk+ii3]), dxi);

      st[ijk] -= grad4(v[ijk-jj1] * interp4(s[ijk-jj3], s[ijk-jj2], s[ijk-jj1], s[ijk    ]),
                       v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]),
                       v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2]),
                       v[ijk+jj2] * interp4(s[ijk    ], s[ijk+jj1], s[ijk+jj2], s[ijk+jj3]), dyi);

      st[ijk] -= grad4x(w[ijk-kk1] * interp4biasbot(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]),
                        w[ijk    ] * interp4       (s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]),
                        w[ijk+kk1] * interp4       (s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]),
                        w[ijk+kk2] * interp4       (s[ijk    ], s[ijk+kk1], s[ijk+kk2], s[ijk+kk3]))
                 * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        st[ijk] -= grad4(u[ijk-ii1] * interp4(s[ijk-ii3], s[ijk-ii2], s[ijk-ii1], s[ijk    ]),
                         u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]),
                         u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2]),
                         u[ijk+ii2] * interp4(s[ijk    ], s[ijk+ii1], s[ijk+ii2], s[ijk+ii3]), dxi);

        st[ijk] -= grad4(v[ijk-jj1] * interp4(s[ijk-jj3], s[ijk-jj2], s[ijk-jj1], s[ijk    ]),
                         v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]),
                         v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2]),
                         v[ijk+jj2] * interp4(s[ijk    ], s[ijk+jj1], s[ijk+jj2], s[ijk+jj3]), dyi);

        st[ijk] -= grad4x(w[ijk-kk1] * interp4(s[ijk-kk3], s[ijk-kk2], s[ijk-kk1], s[ijk    ]),
                          w[ijk    ] * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]),
                          w[ijk+kk1] * interp4(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]),
                          w[ijk+kk2] * interp4(s[ijk    ], s[ijk+kk1], s[ijk+kk2], s[ijk+kk3]))
                   * dzi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      st[ijk] -= grad4(u[ijk-ii1] * interp4(s[ijk-ii3], s[ijk-ii2], s[ijk-ii1], s[ijk    ]),
                       u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]),
                       u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2]),
                       u[ijk+ii2] * interp4(s[ijk    ], s[ijk+ii1], s[ijk+ii2], s[ijk+ii3]), dxi);

      st[ijk] -= grad4(v[ijk-jj1] * interp4(s[ijk-jj3], s[ijk-jj2], s[ijk-jj1], s[ijk    ]),
                       v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]),
                       v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2]),
                       v[ijk+jj2] * interp4(s[ijk    ], s[ijk+jj1], s[ijk+jj2], s[ijk+jj3]), dyi);

      st[ijk] -= grad4x(w[ijk-kk1] * interp4       (s[ijk-kk3], s[ijk-kk2], s[ijk-kk1], s[ijk    ]),
                        w[ijk    ] * interp4       (s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]),
                        w[ijk+kk1] * interp4       (s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]),
                        w[ijk+kk2] * interp4biastop(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]))
                 * dzi4[kend-1];
    }

  return 0;
}

inline double cadvec_g4::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

inline double cadvec_g4::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

inline double cadvec_g4::interp4biasbot(const double a, const double b, const double c, const double d)
{
  return ((5./16.)*a + (15./16.)*b - (5./16.)*c + (1./16)*d);
}

inline double cadvec_g4::interp4biastop(const double a, const double b, const double c, const double d)
{
  return ((5./16.)*d + (15./16.)*c - (5./16.)*b + (1./16)*a);
}

inline double cadvec_g4::grad4(const double a, const double b, const double c, const double d, const double dxi)
{
  return ( -(1./24.)*(d-a) + (27./24.)*(c-b) ) * dxi;
}

inline double cadvec_g4::grad4x(const double a, const double b, const double c, const double d)
{
  return (-(d-a) + 27.*(c-b)); 
}

inline double cadvec_g4::grad4xbiasbot(const double a, const double b, const double c, const double d)
{
  return (-23.*a + 21.*b + 3.*c - d);
}

inline double cadvec_g4::grad4xbiastop(const double a, const double b, const double c, const double d)
{
  return ( 23.*d - 21.*c - 3.*b + a);
}

