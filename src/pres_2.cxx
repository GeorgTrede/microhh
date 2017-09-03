/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "pres_2.h"
#include "defines.h"
#include "model.h"

template<typename TF>
Pres_2<TF>::Pres_2(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) : 
    Pres<TF>(masterin, gridin, fieldsin, inputin)
{
    a = 0;
    c = 0;
    work2d = 0;
    bmati  = 0;
    bmatj  = 0;

#ifdef USECUDA
    a_g = 0;
    c_g = 0;
    work2d_g = 0;
    bmati_g  = 0;
    bmatj_g  = 0;
#endif
}

template<typename TF>
Pres_2<TF>::~Pres_2()
{
    delete[] a;
    delete[] c;
    delete[] work2d;

    delete[] bmati;
    delete[] bmatj;

#ifdef USECUDA
    clear_device();
#endif
}

#ifndef USECUDA
template<typename TF>
void Pres_2<TF>::exec(double dt)
{
    // create the input for the pressure solver
    input(fields.sd["p"]->data,
          fields.mp.at("u")->data, fields.mp.at("v")->data->data, fields.mp.at("w")->data->data,
          fields.mt.at("u")->data, fields.mt.at("v")->data->data, fields.mt.at("w")->data->data,
          grid->dzi, fields.rhoref, fields.rhorefh,
          dt);

    // solve the system
    solve(fields.sd.at("p")->data, fields.atmp.at("tmp1")->data, fields.atmp.at("tmp2")->data,
          grid->dz, fields.rhoref,
          grid->fftini, grid->fftouti, grid->fftinj, grid->fftoutj);

    // get the pressure tendencies from the pressure field
    output(fields.mt.at("u")->data, fields.mt.at("v")->data->data, fields.mt.at("w")->data->data, 
           fields.sd.at("p")->data, grid->dzhi);
}
#endif

#ifndef USECUDA
template<typename TF>
double Pres_2<TF>::check_divergence()
{
    return calc_divergence(fields.u->data, fields.v->data, fields.w->data, grid->dzi,
                           fields.rhoref, fields.rhorefh);
}
#endif

template<typename TF>
void Pres_2<TF>::init()
{
    const int itot = grid->itot;
    const int jtot = grid->jtot;
    const int imax = grid->imax;
    const int jmax = grid->jmax;
    const int kmax = grid->kmax;

    bmati = new double[itot];
    bmatj = new double[jtot];

    a = new double[kmax];
    c = new double[kmax];

    work2d = new double[imax*jmax];
}

template<typename TF>
void Pres_2<TF>::set_values()
{
    const int itot = grid->itot;
    const int jtot = grid->jtot;
    const int kmax = grid->kmax;
    const int kgc  = grid->kgc;

    // compute the modified wave numbers of the 2nd order scheme
    const double dxidxi = 1./(grid->dx*grid->dx);
    const double dyidyi = 1./(grid->dy*grid->dy);

    const double pi = std::acos(-1.);

    for (int j=0; j<jtot/2+1; j++)
        bmatj[j] = 2. * (std::cos(2.*pi*(double)j/(double)jtot)-1.) * dyidyi;

    for (int j=jtot/2+1; j<jtot; j++)
        bmatj[j] = bmatj[jtot-j];

    for (int i=0; i<itot/2+1; i++)
        bmati[i] = 2. * (std::cos(2.*pi*(double)i/(double)itot)-1.) * dxidxi;

    for (int i=itot/2+1; i<itot; i++)
        bmati[i] = bmati[itot-i];

    // create vectors that go into the tridiagonal matrix solver
    for (int k=0; k<kmax; k++)
    {
        a[k] = grid->dz[k+kgc] * fields.rhorefh[k+kgc  ]*grid->dzhi[k+kgc  ];
        c[k] = grid->dz[k+kgc] * fields.rhorefh[k+kgc+1]*grid->dzhi[k+kgc+1];
    }
}

template<typename TF>
void Pres_2<TF>::input(double* restrict p, 
                       double* restrict u , double* restrict v , double* restrict w ,
                       double* restrict ut, double* restrict vt, double* restrict wt,
                       double* restrict dzi, double* restrict rhoref, double* restrict rhorefh,
                       double dt)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const int jjp = grid->imax;
    const int kkp = grid->imax*grid->jmax;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;
    const double dti = 1./dt;

    const int igc = grid->igc;
    const int jgc = grid->jgc;
    const int kgc = grid->kgc;

    // set the cyclic boundary conditions for the tendencies
    grid->boundary_cyclic(ut, Edge::East_west_edge  );
    grid->boundary_cyclic(vt, Edge::North_south_edge);

    // write pressure as a 3d array without ghost cells
    for (int k=0; k<grid->kmax; k++)
        for (int j=0; j<grid->jmax; j++)
#pragma ivdep
            for (int i=0; i<grid->imax; i++)
            {
                const int ijkp = i + j*jjp + k*kkp;
                const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
                p[ijkp] = rhoref[k+kgc] * ( (ut[ijk+ii] + u[ijk+ii] * dti) - (ut[ijk] + u[ijk] * dti) ) * dxi
                        + rhoref[k+kgc] * ( (vt[ijk+jj] + v[ijk+jj] * dti) - (vt[ijk] + v[ijk] * dti) ) * dyi
                        + ( rhorefh[k+kgc+1] * (wt[ijk+kk] + w[ijk+kk] * dti) 
                          - rhorefh[k+kgc  ] * (wt[ijk   ] + w[ijk   ] * dti) ) * dzi[k+kgc];
            }
}

template<typename TF>
void Pres_2<TF>::solve(double* restrict p, double* restrict work3d, double* restrict b,
                   double* restrict dz, double* restrict rhoref,
                   double* restrict fftini, double* restrict fftouti, 
                   double* restrict fftinj, double* restrict fftoutj)
{
    const int imax   = grid->imax;
    const int jmax   = grid->jmax;
    const int kmax   = grid->kmax;
    const int iblock = grid->iblock;
    const int jblock = grid->jblock;
    const int igc    = grid->igc;
    const int jgc    = grid->jgc;
    const int kgc    = grid->kgc;

    int i,j,k,jj,kk,ijk;
    int iindex,jindex;

    grid->fft_forward(p, work3d, fftini, fftouti, fftinj, fftoutj);

    jj = iblock;
    kk = iblock*jblock;

    //for (int i=0; i<itot*jtot; ++i)
    //  printf("%i %e\n",i,p[i]);
    //exit(1);

    // solve the tridiagonal system
    // create vectors that go into the tridiagonal matrix solver
    for (k=0; k<kmax; k++)
        for (j=0; j<jblock; j++)
#pragma ivdep
            for (i=0; i<iblock; i++)
            {
                // swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes
                iindex = master->mpicoordy * iblock + i;
                jindex = master->mpicoordx * jblock + j;

                ijk  = i + j*jj + k*kk;
                b[ijk] = dz[k+kgc]*dz[k+kgc] * rhoref[k+kgc]*(bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
                p[ijk] = dz[k+kgc]*dz[k+kgc] * p[ijk];
            }

    for (j=0; j<jblock; j++)
#pragma ivdep
        for (i=0; i<iblock; i++)
        {
            iindex = master->mpicoordy * iblock + i;
            jindex = master->mpicoordx * jblock + j;

            // substitute BC's
            ijk = i + j*jj;
            b[ijk] += a[0];

            // for wave number 0, which contains average, set pressure at top to zero
            ijk  = i + j*jj + (kmax-1)*kk;
            if (iindex == 0 && jindex == 0)
                b[ijk] -= c[kmax-1];
            // set dp/dz at top to zero
            else
                b[ijk] += c[kmax-1];
        }

    // call tdma solver
    tdma(a, b, c, p, work2d, work3d);

    grid->fft_backward(p, work3d, fftini, fftouti, fftinj, fftoutj);

    jj = imax;
    kk = imax*jmax;

    int ijkp,jjp,kkp;
    jjp = grid->icells;
    kkp = grid->ijcells;

    // put the pressure back onto the original grid including ghost cells
    for (int k=0; k<grid->kmax; k++)
        for (int j=0; j<grid->jmax; j++)
#pragma ivdep
            for (int i=0; i<grid->imax; i++)
            {
                ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;
                ijk  = i + j*jj + k*kk;
                p[ijkp] = work3d[ijk];
            }

    // set the boundary conditions
    // set a zero gradient boundary at the bottom
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            ijk = i + j*jjp + grid->kstart*kkp;
            p[ijk-kkp] = p[ijk];
        }

    // set the cyclic boundary conditions
    grid->boundary_cyclic(p);
}

template<typename TF>
void Pres_2<TF>::output(double* restrict ut, double* restrict vt, double* restrict wt, 
                        double* restrict p , double* restrict dzhi)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    for (int k=grid->kstart; k<grid->kend; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ut[ijk] -= (p[ijk] - p[ijk-ii]) * dxi;
                vt[ijk] -= (p[ijk] - p[ijk-jj]) * dyi;
                wt[ijk] -= (p[ijk] - p[ijk-kk]) * dzhi[k];
            }
}

// tridiagonal matrix solver, taken from Numerical Recipes, Press
template<typename TF>
void Pres_2<TF>::tdma(double* restrict a, double* restrict b, double* restrict c, 
                      double* restrict p, double* restrict work2d, double* restrict work3d)

{
    int i,j,k,jj,kk,ijk,ij;
    int iblock,jblock,kmax;

    iblock = grid->iblock;
    jblock = grid->jblock;
    kmax = grid->kmax;

    jj = iblock;
    kk = iblock*jblock;

    for (j=0;j<jblock;j++)
#pragma ivdep
        for (i=0;i<iblock;i++)
        {
            ij = i + j*jj;
            work2d[ij] = b[ij];
        }

    for (j=0;j<jblock;j++)
#pragma ivdep
        for (i=0;i<iblock;i++)
        {
            ij = i + j*jj;
            p[ij] /= work2d[ij];
        }

    for (k=1; k<kmax; k++)
    {
        for (j=0;j<jblock;j++)
#pragma ivdep
            for (i=0;i<iblock;i++)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + k*kk;
                work3d[ijk] = c[k-1] / work2d[ij];
            }
        for (j=0;j<jblock;j++)
#pragma ivdep
            for (i=0;i<iblock;i++)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + k*kk;
                work2d[ij] = b[ijk] - a[k]*work3d[ijk];
            }
        for (j=0;j<jblock;j++)
#pragma ivdep
            for (i=0;i<iblock;i++)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + k*kk;
                p[ijk] -= a[k]*p[ijk-kk];
                p[ijk] /= work2d[ij];
            }
    }

    for (k=kmax-2; k>=0; k--)
        for (j=0;j<jblock;j++)
#pragma ivdep
            for (i=0;i<iblock;i++)
            {
                ijk = i + j*jj + k*kk;
                p[ijk] -= work3d[ijk+kk]*p[ijk+kk];
            }
}

#ifndef USECUDA
template<typename TF>
double Pres_2<TF>::calc_divergence(double* restrict u, double* restrict v, double* restrict w, double* restrict dzi,
                                   double* restrict rhoref, double* restrict rhorefh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    double div    = 0.;
    double divmax = 0.;

    for (int k=grid->kstart; k<grid->kend; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                div = rhoref[k]*((u[ijk+ii]-u[ijk])*dxi + (v[ijk+jj]-v[ijk])*dyi) 
                    + (rhorefh[k+1]*w[ijk+kk]-rhorefh[k]*w[ijk])*dzi[k];

                divmax = std::max(divmax, std::abs(div));
            }

    grid->get_max(&divmax);

    return divmax;
}
#endif
