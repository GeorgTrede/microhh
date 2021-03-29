/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "boundary_surface.h"
#include "boundary_surface_functions.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "master.h"
#include "cross.h"
#include "column.h"
#include "monin_obukhov.h"
#include "fast_math.h"

namespace
{
    // Make a shortcut in the file scope.
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace bsf = Boundary_surface_functions;

    // Size of the lookup table.
    const int nzL = 10000; // Size of the lookup table for MO iterations.

    template<typename TF>
    TF find_zL(
            const float* const restrict zL,
            const float* const restrict f,
            int& n, const float Ri)
    {
        // Determine search direction. All checks are at float accuracy.
        if ( (f[n]-Ri) > 0.f )
            while ( (f[n-1]-Ri) > 0.f && n > 0.f) { --n; }
        else
            while ( (f[n]-Ri) < 0.f && n < (nzL-1) ) { ++n; }

        const TF zL0 = (n == 0 || n == nzL-1) ? zL[n] : zL[n-1] + (Ri-f[n-1]) / (f[n]-f[n-1]) * (zL[n]-zL[n-1]);

        return zL0;
    }

    template<typename TF>
    TF calc_obuk_noslip_flux_lookup(
            const float* restrict zL, const float* restrict f,
            int& n,
            const TF du, const TF bfluxbot, const TF zsl)
    {
        // Calculate the appropriate Richardson number and reduce precision.
        const float Ri = -Constants::kappa<TF> * bfluxbot * zsl / fm::pow3(du);
        return zsl/find_zL<TF>(zL, f, n, Ri);
    }

    template<typename TF>
    TF calc_obuk_noslip_dirichlet_lookup(
            const float* restrict zL, const float* restrict f,
            int& n,
            const TF du, const TF db, const TF zsl)
    {
        // Calculate the appropriate Richardson number and reduce precision.
        const float Ri = Constants::kappa<TF> * db * zsl / fm::pow2(du);
        return zsl/find_zL<TF>(zL, f, n, Ri);
    }

    template<typename TF>
    void set_bc(
            TF* const restrict a, TF* const restrict agrad, TF* const restrict aflux,
            const Boundary_type sw, const TF aval, const TF visc, const TF offset,
            const int icells, const int jcells)
    {
        const int jj = icells;

        if (sw == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    a[ij] = aval - offset;
                }
        }
        else if (sw == Boundary_type::Neumann_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    agrad[ij] = aval;
                    aflux[ij] = -aval*visc;
                }
        }
        else if (sw == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    aflux[ij] = aval;
                    agrad[ij] = -aval/visc;
                }
        }
    }

    template<typename TF, bool sw_constant_z0>
    void stability(
            TF* const restrict ustar,
            TF* const restrict obuk,
            const TF* const restrict bfluxbot,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict b,
            const TF* const restrict ubot,
            const TF* const restrict vbot,
            const TF* const restrict bbot,
            TF* const restrict dutot,
            const TF* const restrict z,
            const TF* const restrict z0m,
            const TF* const restrict z0h,
            const float* const zL_sl,
            const float* const f_sl,
            int* const nobuk,
            const TF db_ref,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int jcells, const int kk,
            Boundary_type mbcbot, Boundary_type thermobc,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int ii = 1;
        const int jj = icells;

        // Calculate total wind.
        const TF minval = 1.e-1;

        // First, interpolate the wind to the scalar location.
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                const TF du2 = fm::pow2(TF(0.5)*(u[ijk] + u[ijk+ii]) - TF(0.5)*(ubot[ij] + ubot[ij+ii]))
                             + fm::pow2(TF(0.5)*(v[ijk] + v[ijk+jj]) - TF(0.5)*(vbot[ij] + vbot[ij+jj]));
                // Prevent the absolute wind gradient from reaching values less than 0.01 m/s,
                // otherwise evisc at k = kstart blows up
                dutot[ij] = std::max(std::pow(du2, TF(0.5)), minval);
            }

        boundary_cyclic.exec_2d(dutot);

        // Calculate Obukhov length
        // Case 1: fixed buoyancy flux and fixed ustar
        if (mbcbot == Boundary_type::Ustar_type && thermobc == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    obuk[ij] = -fm::pow3(ustar[ij]) / (Constants::kappa<TF>*bfluxbot[ij]);
                }
        }
        // Case 2: fixed buoyancy surface value and free ustar
        else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;

                    // Switch between the iterative and lookup solver
                    if (sw_constant_z0)
                        obuk[ij] = calc_obuk_noslip_flux_lookup(
                                zL_sl, f_sl, nobuk[ij], dutot[ij], bfluxbot[ij], z[kstart]);
                    else
                        obuk[ij] = bsf::calc_obuk_noslip_flux_iterative(
                                obuk[ij], dutot[ij], bfluxbot[ij], z[kstart], z0m[ij]);

                    ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m[ij], obuk[ij]);
                }
        }
        else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    const TF db = b[ijk] - bbot[ij] + db_ref;

                    // Switch between the iterative and lookup solver
                    if (sw_constant_z0)
                        obuk[ij] = calc_obuk_noslip_dirichlet_lookup(
                                zL_sl, f_sl, nobuk[ij], dutot[ij], db, z[kstart]);
                    else
                        obuk[ij] = bsf::calc_obuk_noslip_dirichlet_iterative(
                                obuk[ij], dutot[ij], db, z[kstart], z0m[ij], z0h[ij]);

                    ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m[ij], obuk[ij]);
                }
        }
    }

    template<typename TF>
    void stability_neutral(
            TF* const restrict ustar,
            TF* const restrict obuk,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict ubot,
            const TF* const restrict vbot,
            TF* const restrict dutot,
            const TF* const restrict z,
            const TF* const restrict z0m,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int jcells, const int kk,
            Boundary_type mbcbot,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int ii = 1;
        const int jj = icells;

        // calculate total wind
        TF du2;
        const TF minval = 1.e-1;

        // first, interpolate the wind to the scalar location
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                du2 = fm::pow2(TF(0.5)*(u[ijk] + u[ijk+ii]) - TF(0.5)*(ubot[ij] + ubot[ij+ii]))
                    + fm::pow2(TF(0.5)*(v[ijk] + v[ijk+jj]) - TF(0.5)*(vbot[ij] + vbot[ij+jj]));
                // prevent the absolute wind gradient from reaching values less than 0.01 m/s,
                // otherwise evisc at k = kstart blows up
                dutot[ij] = std::max(std::pow(du2, TF(0.5)), minval);
            }

        boundary_cyclic.exec_2d(dutot);

        // set the Obukhov length to a very large negative number
        // case 1: fixed buoyancy flux and fixed ustar
        if (mbcbot == Boundary_type::Ustar_type)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    obuk[ij] = -Constants::dbig;
                }
        }
        // case 2: free ustar
        else if (mbcbot == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    obuk [ij] = -Constants::dbig;
                    ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m[ij], obuk[ij]);
                }
        }
    }

    template<typename TF>
    void surfm(
            TF* const restrict ufluxbot,
            TF* const restrict vfluxbot,
            TF* const restrict ugradbot,
            TF* const restrict vgradbot,
            const TF* const restrict ustar,
            const TF* const restrict obuk,
            const TF* const restrict u, const TF* const restrict ubot,
            const TF* const restrict v, const TF* const restrict vbot,
            const TF* const restrict z0m,
            const TF zsl, const Boundary_type bcbot,
            const int istart, const int iend,
            const int jstart, const int jend, const int kstart,
            const int icells, const int jcells, const int kk,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int ii = 1;
        const int jj = icells;

        // the surface value is known, calculate the flux and gradient
        if (bcbot == Boundary_type::Dirichlet_type)
        {
            // first calculate the surface value
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;

                    // interpolate the whole stability function rather than ustar or obuk
                    ufluxbot[ij] = -(u[ijk]-ubot[ij])*TF(0.5)*
                        (ustar[ij-ii]*most::fm(zsl, z0m[ij-ii], obuk[ij-ii]) + ustar[ij]*most::fm(zsl, z0m[ij], obuk[ij]));
                    vfluxbot[ij] = -(v[ijk]-vbot[ij])*TF(0.5)*
                        (ustar[ij-jj]*most::fm(zsl, z0m[ij-jj], obuk[ij-jj]) + ustar[ij]*most::fm(zsl, z0m[ij], obuk[ij]));
                }

            boundary_cyclic.exec_2d(ufluxbot);
            boundary_cyclic.exec_2d(vfluxbot);
        }

        // the flux is known, calculate the surface value and gradient
        else if (bcbot == Boundary_type::Ustar_type)
        {
            // first redistribute ustar over the two flux components
            const TF minval = 1.e-2;

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;

                    const TF vonu2 = std::max(minval, TF(0.25)*(
                                fm::pow2(v[ijk-ii]-vbot[ij-ii]) + fm::pow2(v[ijk-ii+jj]-vbot[ij-ii+jj])
                              + fm::pow2(v[ijk   ]-vbot[ij   ]) + fm::pow2(v[ijk   +jj]-vbot[ij   +jj])) );
                    const TF uonv2 = std::max(minval, TF(0.25)*(
                                fm::pow2(u[ijk-jj]-ubot[ij-jj]) + fm::pow2(u[ijk+ii-jj]-ubot[ij+ii-jj])
                              + fm::pow2(u[ijk   ]-ubot[ij   ]) + fm::pow2(u[ijk+ii   ]-ubot[ij+ii   ])) );

                    const TF u2 = std::max(minval, fm::pow2(u[ijk]-ubot[ij]) );
                    const TF v2 = std::max(minval, fm::pow2(v[ijk]-vbot[ij]) );

                    const TF ustaronu4 = TF(0.5)*(fm::pow4(ustar[ij-ii]) + fm::pow4(ustar[ij]));
                    const TF ustaronv4 = TF(0.5)*(fm::pow4(ustar[ij-jj]) + fm::pow4(ustar[ij]));

                    ufluxbot[ij] = -copysign(TF(1), u[ijk]-ubot[ij]) * std::pow(ustaronu4 / (TF(1) + vonu2 / u2), TF(0.5));
                    vfluxbot[ij] = -copysign(TF(1), v[ijk]-vbot[ij]) * std::pow(ustaronv4 / (TF(1) + uonv2 / v2), TF(0.5));
                }

            boundary_cyclic.exec_2d(ufluxbot);
            boundary_cyclic.exec_2d(vfluxbot);

            // CvH: I think that the problem is not closed, since both the fluxes and the surface values
            // of u and v are unknown. You have to assume a no slip in order to get the fluxes and therefore
            // should not update the surface values with those that belong to the flux. This procedure needs
            // to be checked more carefully.
            /*
            // calculate the surface values
            for (int j=grid->jstart; j<grid->jend; ++j)
                #pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                ij  = i + j*jj;
                ijk = i + j*jj + kstart*kk;
                // interpolate the whole stability function rather than ustar or obuk
                ubot[ij] = 0.;// ufluxbot[ij] / (0.5*(ustar[ij-ii]*fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*fm(zsl, z0m, obuk[ij]))) + u[ijk];
                vbot[ij] = 0.;// vfluxbot[ij] / (0.5*(ustar[ij-jj]*fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*fm(zsl, z0m, obuk[ij]))) + v[ijk];
            }

            grid->boundary_cyclic_2d(ubot);
            grid->boundary_cyclic_2d(vbot);
            */
        }

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                // use the linearly interpolated grad, rather than the MO grad,
                // to prevent giving unresolvable gradients to advection schemes
                // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0m*ustar[ij]) * phih(zsl/obuk[ij]);
                ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
                vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
            }
    }

    template<typename TF>
    void surfs(
            TF* const restrict varbot,
            TF* const restrict vargradbot,
            TF* const restrict varfluxbot,
            const TF* const restrict ustar,
            const TF* const restrict obuk,
            const TF* const restrict var,
            const TF* const restrict z0h,
            const TF zsl, const Boundary_type bcbot,
            const int istart, const int iend,
            const int jstart, const int jend, const int kstart,
            const int icells, const int jcells, const int kk,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;

        // the surface value is known, calculate the flux and gradient
        if (bcbot == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*most::fh(zsl, z0h[ij], obuk[ij]);
                    // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
                    // use the linearly interpolated grad, rather than the MO grad,
                    // to prevent giving unresolvable gradients to advection schemes
                    vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
                }
        }
        else if (bcbot == Boundary_type::Flux_type)
        {
            // the flux is known, calculate the surface value and gradient
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    varbot[ij] = varfluxbot[ij] / (ustar[ij]*most::fh(zsl, z0h[ij], obuk[ij])) + var[ijk];
                    // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
                    // use the linearly interpolated grad, rather than the MO grad,
                    // to prevent giving unresolvable gradients to advection schemes
                    vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
                }
        }
    }

    template<typename TF>
    void calc_ra(
            TF* const restrict ra,
            const TF* const restrict ustar,
            const TF* const restrict obuk,
            const TF* const restrict z0h,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                ra[ij]  = TF(1) / (ustar[ij] * most::fh(zsl, z0h[ij], obuk[ij]));
            }
    }
}

template<typename TF>
Boundary_surface<TF>::Boundary_surface(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Boundary<TF>(masterin, gridin, fieldsin, inputin)
{
    swboundary = "surface";

    #ifdef USECUDA
    ustar_g = 0;
    obuk_g  = 0;
    nobuk_g = 0;
    zL_sl_g = 0;
    f_sl_g  = 0;
    #endif
}

template<typename TF>
Boundary_surface<TF>::~Boundary_surface()
{
    #ifdef USECUDA
    clear_device();
    #endif
}

template<typename TF>
void Boundary_surface<TF>::create(
        Input& input, Netcdf_handle& input_nc,
        Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross)
{
    const std::string group_name = "default";
    Boundary<TF>::process_time_dependent(input, input_nc);

    // add variables to the statistics
    if (stats.get_switch())
    {
        stats.add_time_series("ustar", "Surface friction velocity", "m s-1", group_name);
        stats.add_time_series("obuk", "Obukhov length", "m", group_name);
    }

    if (column.get_switch())
    {
        column.add_time_series("ustar", "Surface friction velocity", "m s-1");
        column.add_time_series("obuk", "Obukhov length", "m");
    }

    if (cross.get_switch())
    {
        const std::vector<std::string> allowed_crossvars = {"ustar", "obuk", "ra"};
        cross_list = cross.get_enabled_variables(allowed_crossvars);
    }
}

template<typename TF>
void Boundary_surface<TF>::init(Input& inputin, Thermo<TF>& thermo)
{
    // 1. Process the boundary conditions now all fields are registered.
    process_bcs(inputin);

    // 2. Read and check the boundary_surface specific settings.
    process_input(inputin, thermo);

    // 3. Allocate and initialize the 2D surface fields.
    init_surface(inputin);

    // 4. Initialize the boundary cyclic.
    boundary_cyclic.init();
}

template<typename TF>
void Boundary_surface<TF>::process_input(Input& inputin, Thermo<TF>& thermo)
{
    // Switch between heterogeneous and homogeneous z0's
    sw_constant_z0 = inputin.get_item<bool>("boundary", "sw_constant_z0", "", true);

    // crash in case fixed gradient is prescribed
    if (mbcbot == Boundary_type::Neumann_type)
    {
        std::string msg = "Neumann bc is not supported in surface model";
        throw std::runtime_error(msg);
    }
    // read the ustar value only if fixed fluxes are prescribed
    else if (mbcbot == Boundary_type::Ustar_type)
        ustarin = inputin.get_item<TF>("boundary", "ustar", "");

    // process the scalars
    for (auto& it : sbc)
    {
        // crash in case fixed gradient is prescribed
        if (it.second.bcbot == Boundary_type::Neumann_type)
        {
            std::string msg = "Fixed Gradient bc is not supported in surface model";
            throw std::runtime_error(msg);
        }

        // crash in case of fixed momentum flux and dirichlet bc for scalar
        if (it.second.bcbot == Boundary_type::Dirichlet_type && mbcbot == Boundary_type::Ustar_type)
        {
            std::string msg = "Fixed Ustar bc in combination with Dirichlet bc for scalars is not supported";
            throw std::runtime_error(msg);
        }
    }

    // check whether the prognostic thermo vars are of the same type
    std::vector<std::string> thermolist;
    thermo.get_prog_vars(thermolist);

    auto it = thermolist.begin();

    // save the bc of the first thermo field in case thermo is enabled
    if (it != thermolist.end())
        thermobc = sbc[*it].bcbot;
    else
        // Set the thermobc to Flux_type to avoid ininitialized errors.
        thermobc = Boundary_type::Flux_type;

    while (it != thermolist.end())
    {
        if (sbc[*it].bcbot != thermobc)
        {

            std::string msg = "All thermo variables need to have the same bc type";
            throw std::runtime_error(msg);
        }
        ++it;
    }
}

template<typename TF>
void Boundary_surface<TF>::init_surface(Input& input)
{
    auto& gd = grid.get_grid_data();

    obuk.resize(gd.ijcells);
    ustar.resize(gd.ijcells);

    if (sw_constant_z0)
        nobuk.resize(gd.ijcells);

    z0m.resize(gd.ijcells);
    z0h.resize(gd.ijcells);

    if (sw_constant_z0)
    {
        const TF z0m_hom = input.get_item<TF>("boundary", "z0m", "");
        const TF z0h_hom = input.get_item<TF>("boundary", "z0h", "");

        std::fill(z0m.begin(), z0m.end(), z0m_hom);
        std::fill(z0h.begin(), z0h.end(), z0h_hom);
    }

    // Initialize the obukhov length on a small number.
    std::fill(obuk.begin(),  obuk.end(), Constants::dsmall);
    if (sw_constant_z0)
        std::fill(nobuk.begin(), nobuk.end(), 0);
}

template<typename TF>
void Boundary_surface<TF>::load(const int iotime)
{
    // Obukhov length restart files are only needed for the iterative solver
    if (!sw_constant_z0)
    {
        auto tmp1 = fields.get_tmp();
        int nerror = 0;

        auto load_2d_field = [&](
                TF* const restrict field, const std::string& name, const int itime)
        {
            char filename[256];
            std::sprintf(filename, "%s.%07d", name.c_str(), itime);
            master.print_message("Loading \"%s\" ... ", filename);

            if (field3d_io.load_xy_slice(
                    field, tmp1->fld.data(),
                    filename))
            {
                master.print_message("FAILED\n");
                nerror += 1;
            }
            else
                master.print_message("OK\n");

            boundary_cyclic.exec_2d(field);
        };

        // Read Obukhov length
        load_2d_field(obuk.data(), "obuk", iotime);

        // Read spatial z0 fields
        load_2d_field(z0m.data(), "z0m", 0);
        load_2d_field(z0h.data(), "z0h", 0);

        master.sum(&nerror, 1);
        if (nerror)
            throw std::runtime_error("Error loading field(s)");

        fields.release_tmp(tmp1);
    }
}

template<typename TF>
void Boundary_surface<TF>::save(const int iotime)
{
    // Obukhov length restart files are only needed for the iterative solver
    if (!sw_constant_z0)
    {
        auto tmp1 = fields.get_tmp();

        char filename[256];
        std::sprintf(filename, "%s.%07d", "obuk", iotime);
        master.print_message("Saving \"%s\" ... ", filename);

        const int kslice = 0;
        int nerror = 0;
        if (field3d_io.save_xy_slice(
                obuk.data(), tmp1->fld.data(), filename, kslice))
        {
            master.print_message("FAILED\n");
            nerror += 1;
        }
        else
            master.print_message("OK\n");

        master.sum(&nerror, 1);
        if (nerror)
            throw std::runtime_error("Error saving field(s)");

        fields.release_tmp(tmp1);
    }
}

template<typename TF>
void Boundary_surface<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    auto& gd = grid.get_grid_data();
    auto tmp1 = fields.get_tmp();

    for (auto& it : cross_list)
    {
        if (it == "ustar")
            cross.cross_plane(ustar.data(), "ustar", iotime);
        else if (it == "obuk")
            cross.cross_plane(obuk.data(), "obuk", iotime);
        else if (it == "ra")
        {
            calc_ra(tmp1->flux_bot.data(), ustar.data(), obuk.data(),
                    z0h.data(), gd.z[gd.kstart], gd.istart,
                    gd.iend, gd.jstart, gd.jend, gd.icells);
            cross.cross_plane(tmp1->flux_bot.data(), "ra", iotime);
        }
    }

    fields.release_tmp(tmp1);
}

template<typename TF>
void Boundary_surface<TF>::exec_stats(Stats<TF>& stats)
{
    const TF no_offset = 0.;
    stats.calc_stats_2d("obuk", obuk, no_offset);
    stats.calc_stats_2d("ustar", ustar, no_offset);
}

#ifndef USECUDA
template<typename TF>
void Boundary_surface<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    column.calc_time_series("obuk", obuk.data(), no_offset);
    column.calc_time_series("ustar", ustar.data(), no_offset);
}
#endif

template<typename TF>
void Boundary_surface<TF>::set_values()
{
    auto& gd = grid.get_grid_data();

    // Call the base class function.
    Boundary<TF>::set_values();

    // Override the boundary settings in order to enforce dirichlet BC for surface model.
    set_bc<TF>(
            fields.mp.at("u")->fld_bot.data(),
            fields.mp.at("u")->grad_bot.data(),
            fields.mp.at("u")->flux_bot.data(),
            Boundary_type::Dirichlet_type, ubot,
            fields.visc, grid.utrans,
            gd.icells, gd.jcells);

    set_bc<TF>(
            fields.mp.at("v")->fld_bot.data(),
            fields.mp.at("v")->grad_bot.data(),
            fields.mp.at("v")->flux_bot.data(),
            Boundary_type::Dirichlet_type, vbot,
            fields.visc, grid.vtrans,
            gd.icells, gd.jcells);

    // in case the momentum has a fixed ustar, set the value to that of the input
    if (mbcbot == Boundary_type::Ustar_type)
        set_ustar();

    // Prepare the lookup table for the surface solver
    if (sw_constant_z0)
        init_solver();
}

template<typename TF>
void Boundary_surface<TF>::set_ustar()
{
    auto& gd = grid.get_grid_data();
    const int jj = gd.icells;

    set_bc<TF>(
            fields.mp.at("u")->fld_bot.data(),
            fields.mp.at("u")->grad_bot.data(),
            fields.mp.at("u")->flux_bot.data(),
            mbcbot, ubot, fields.visc, grid.utrans,
            gd.icells, gd.jcells);

    set_bc<TF>(
            fields.mp.at("v")->fld_bot.data(),
            fields.mp.at("v")->grad_bot.data(),
            fields.mp.at("v")->flux_bot.data(),
            mbcbot, vbot, fields.visc, grid.vtrans,
            gd.icells, gd.jcells);

    for (int j=0; j<gd.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<gd.icells; ++i)
        {
            const int ij = i + j*jj;
            // Limit ustar at 1e-4 to avoid zero divisions.
            ustar[ij] = std::max(static_cast<TF>(0.0001), ustarin);
        }
}

// Prepare the surface layer solver.
template<typename TF>
void Boundary_surface<TF>::init_solver()
{
    auto& gd = grid.get_grid_data();

    zL_sl.resize(nzL);
    f_sl.resize(nzL);

    std::vector<TF> zL_tmp(nzL);

    // BvS: TMP
    const TF z0m_hom = z0m[0];
    const TF z0h_hom = z0h[0];

    // Calculate the non-streched part between -5 to 10 z/L with 9/10 of the points,
    // and stretch up to -1e4 in the negative limit.
    // Alter next three values in case the range need to be changed.
    const TF zL_min = -1.e4;
    const TF zLrange_min = -5.;
    const TF zLrange_max = 10.;

    TF dzL = (zLrange_max - zLrange_min) / (9.*nzL/10.-1.);
    zL_tmp[0] = -zLrange_max;
    for (int n=1; n<9*nzL/10; ++n)
        zL_tmp[n] = zL_tmp[n-1] + dzL;

    // Stretch the remainder of the z/L values far down for free convection.
    const TF zLend = -(zL_min - zLrange_min);

    // Find stretching that ends up at the correct value using geometric progression.
    TF r  = 1.01;
    TF r0 = Constants::dhuge;
    while (std::abs( (r-r0)/r0 ) > 1.e-10)
    {
        r0 = r;
        r  = std::pow( 1. - (zLend/dzL)*(1.-r), (1./ (nzL/10.) ) );
    }

    for (int n=9*nzL/10; n<nzL; ++n)
    {
        zL_tmp[n] = zL_tmp[n-1] + dzL;
        dzL *= r;
    }

    // Calculate the final array and delete the temporary array.
    for (int n=0; n<nzL; ++n)
        zL_sl[n] = -zL_tmp[nzL-n-1];

    // Calculate the evaluation function.
    if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Flux_type)
    {
        const TF zsl = gd.z[gd.kstart];
        for (int n=0; n<nzL; ++n)
            f_sl[n] = zL_sl[n] * std::pow(most::fm(zsl, z0m_hom, zsl/zL_sl[n]), 3);
    }
    else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Dirichlet_type)
    {
        const TF zsl = gd.z[gd.kstart];
        for (int n=0; n<nzL; ++n)
            f_sl[n] = zL_sl[n] * std::pow(most::fm(zsl, z0m_hom, zsl/zL_sl[n]), 2) / most::fh(zsl, z0h_hom, zsl/zL_sl[n]);
    }
}

#ifndef USECUDA
template<typename TF>
void Boundary_surface<TF>::calc_mo_stability(
        Thermo<TF>& thermo, Land_surface<TF>& lsm)
{
    auto& gd = grid.get_grid_data();

    // Start with retrieving the stability information.
    if (thermo.get_switch() == "0")
    {
        auto dutot = fields.get_tmp();
        stability_neutral(
                ustar.data(), obuk.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                fields.mp.at("u")->fld_bot.data(), fields.mp.at("v")->fld_bot.data(),
                dutot->fld.data(), gd.z.data(), z0m.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend, gd.kstart,
                gd.icells, gd.jcells, gd.ijcells,
                mbcbot, boundary_cyclic);
        fields.release_tmp(dutot);
    }
    else
    {
        auto buoy = fields.get_tmp();
        auto tmp = fields.get_tmp();

        thermo.get_buoyancy_surf(*buoy, false);
        const TF db_ref = thermo.get_db_ref();

        if (sw_constant_z0)
            stability<TF, true>(
                    ustar.data(), obuk.data(),
                    buoy->flux_bot.data(),
                    fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                    buoy->fld.data(), fields.mp.at("u")->fld_bot.data(),
                    fields.mp.at("v")->fld_bot.data(), buoy->fld_bot.data(),
                    tmp->fld.data(), gd.z.data(),
                    z0m.data(), z0h.data(),
                    zL_sl.data(), f_sl.data(), nobuk.data(), db_ref,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend, gd.kstart,
                    gd.icells, gd.jcells, gd.ijcells,
                    mbcbot, thermobc, boundary_cyclic);
        else
            stability<TF, false>(
                    ustar.data(), obuk.data(),
                    buoy->flux_bot.data(),
                    fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), buoy->fld.data(),
                    fields.mp.at("u")->fld_bot.data(), fields.mp.at("v")->fld_bot.data(), buoy->fld_bot.data(),
                    tmp->fld.data(), gd.z.data(),
                    z0m.data(), z0h.data(),
                    zL_sl.data(), f_sl.data(), nobuk.data(), db_ref,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart,
                    gd.icells, gd.jcells, gd.ijcells,
                    mbcbot, thermobc, boundary_cyclic);

        fields.release_tmp(buoy);
        fields.release_tmp(tmp);
    }
}

template<typename TF>
void Boundary_surface<TF>::calc_mo_bcs_momentum(
        Thermo<TF>& thermo, Land_surface<TF>& lsm)
{
    auto& gd = grid.get_grid_data();

    // Calculate the surface value, gradient and flux depending on the chosen boundary condition.
    surfm(fields.mp.at("u")->flux_bot.data(),
          fields.mp.at("v")->flux_bot.data(),
          fields.mp.at("u")->grad_bot.data(),
          fields.mp.at("v")->grad_bot.data(),
          ustar.data(), obuk.data(),
          fields.mp.at("u")->fld.data(), fields.mp.at("u")->fld_bot.data(),
          fields.mp.at("v")->fld.data(), fields.mp.at("v")->fld_bot.data(),
          z0m.data(), gd.z[gd.kstart], mbcbot,
          gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart,
          gd.icells, gd.jcells, gd.ijcells,
          boundary_cyclic);
}

template<typename TF>
void Boundary_surface<TF>::calc_mo_bcs_scalars(
        Thermo<TF>& thermo, Land_surface<TF>& lsm)
{
    auto& gd = grid.get_grid_data();

    for (auto& it : fields.sp)
        surfs(it.second->fld_bot.data(),
              it.second->grad_bot.data(),
              it.second->flux_bot.data(),
              ustar.data(), obuk.data(),
              it.second->fld.data(), z0h.data(),
              gd.z[gd.kstart], sbc.at(it.first).bcbot,
              gd.istart, gd.iend,
              gd.jstart, gd.jend, gd.kstart,
              gd.icells, gd.jcells, gd.ijcells,
              boundary_cyclic);
}

template<typename TF>
void Boundary_surface<TF>::get_ra(Field3d<TF>& fld)
{
    auto& gd = grid.get_grid_data();

    calc_ra(fld.flux_bot.data(),
            ustar.data(),
            obuk.data(),
            z0h.data(),
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
}

template<typename TF>
void Boundary_surface<TF>::get_duvdz(
        std::vector<TF>& dudz, std::vector<TF>& dvdz)
{
    auto& gd = grid.get_grid_data();

    Boundary_surface_functions::calc_duvdz(
            dudz.data(), dvdz.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("u")->fld_bot.data(),
            fields.mp.at("v")->fld_bot.data(),
            fields.mp.at("u")->flux_bot.data(),
            fields.mp.at("v")->flux_bot.data(),
            ustar.data(), obuk.data(), z0m.data(),
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart,
            gd.icells, gd.ijcells);
}

template<typename TF>
void Boundary_surface<TF>::get_dbdz(
        std::vector<TF>& dbdz, std::vector<TF>& bfluxbot)
{
    auto& gd = grid.get_grid_data();

    Boundary_surface_functions::calc_dbdz(
            dbdz.data(), bfluxbot.data(),
            ustar.data(), obuk.data(),
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
}
#endif

template<typename TF>
void Boundary_surface<TF>::update_slave_bcs()
{
    // This function does nothing when the surface model is enabled, because
    // the fields are computed by the surface model in update_bcs.
}

template class Boundary_surface<double>;
template class Boundary_surface<float>;
