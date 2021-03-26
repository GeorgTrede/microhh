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
#include "grid.h"
#include "fields.h"
#include "input.h"

#include "boundary_surface_tiled.h"
#include "boundary_surface_functions.h"

#include "land_surface.h"
#include "fast_math.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "monin_obukhov.h"

#include "cross.h"
#include "column.h"

namespace
{
    // Make a shortcut in the file scope.
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace bsf = Boundary_surface_functions;

    template<typename TF>
    void init_tile(MO_surface_tile<TF>& tile, const int ijcells)
    {
        tile.obuk.resize(ijcells);
        tile.ustar.resize(ijcells);
        tile.z0m.resize(ijcells);
        tile.z0h.resize(ijcells);
    }

    template<typename TF>
    void calc_tiled_mean(
            TF* const restrict fld,
            const TF* const restrict fld_veg,
            const TF* const restrict fld_soil,
            const TF* const restrict fld_wet,
            const TF* const restrict c_veg,
            const TF* const restrict c_soil,
            const TF* const restrict c_wet,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells, const int jcells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;

                fld[ij] = c_veg [ij] * fld_veg [ij] +
                          c_soil[ij] * fld_soil[ij] +
                          c_wet [ij] * fld_wet [ij];
            }
    }

    template<typename TF>
    void calc_bulk_obuk(
            TF* const restrict obuk,
            const TF* const restrict bfluxbot,
            const TF* const restrict ustar,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                obuk[ij] = -fm::pow3(ustar[ij]) / (Constants::kappa<TF> * bfluxbot[ij]);
            }
    }

    template<typename TF>
    void calc_dutot(
            TF* const restrict dutot,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict ubot,
            const TF* const restrict vbot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int jcells,
            const int ijcells)
    {
        const TF minval = 1.e-1;
        const int ii = 1;
        const int jj = icells;

        // Interpolate the wind to the scalar location.
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kstart*ijcells;

                const TF du2 = fm::pow2(TF(0.5)*(u[ijk] + u[ijk+ii]) - TF(0.5)*(ubot[ij] + ubot[ij+ii]))
                             + fm::pow2(TF(0.5)*(v[ijk] + v[ijk+jj]) - TF(0.5)*(vbot[ij] + vbot[ij+jj]));

                // Prevent the absolute wind gradient from reaching values less than 0.01 m/s,
                // otherwise evisc at k = kstart blows up
                dutot[ij] = std::max(std::pow(du2, TF(0.5)), minval);
            }
    }

    template<typename TF>
    void stability(
            TF* const restrict ustar,
            TF* const restrict obuk,
            const TF* const restrict dutot,
            const TF* const restrict b,
            const TF* const restrict bbot,
            const TF* const restrict z0m,
            const TF* const restrict z0h,
            const TF db_ref,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int jcells,
            const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                const TF db = b[ijk] - bbot[ij] + db_ref;

                obuk[ij] = bsf::calc_obuk_noslip_dirichlet_iterative(
                        obuk[ij], dutot[ij], db, zsl, z0m[ij], z0h[ij]);
                ustar[ij] = dutot[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
            }
    }

    template<typename TF>
    void surfm(
            TF* const restrict ufluxbot,
            TF* const restrict vfluxbot,
            TF* const restrict ugradbot,
            TF* const restrict vgradbot,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict ubot,
            const TF* const restrict vbot,
            const TF* const restrict ustar,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int icells,
            const int jcells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

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

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                // Use the linearly interpolated grad, rather than the MO grad,
                // to prevent giving unresolvable gradients to advection schemes
                ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
                vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
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

    template<typename TF>
    void calc_effective_ra(
            TF* const restrict ra,
            const TF* const restrict thl,
            const TF* const restrict qt,
            const TF* const restrict thl_bot,
            const TF* const restrict qt_bot,
            const TF* const restrict thl_fluxbot,
            const TF* const restrict qt_fluxbot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kstart*ijcells;

                ra[ij] = (thl_bot[ij] - thl[ijk]) / thl_fluxbot[ij];
            }
    }

//    template<typename TF>
//    void surfs(
//            TF* const restrict varbot,
//            TF* const restrict vargradbot,
//            TF* const restrict varfluxbot,
//            const TF* const restrict ustar,
//            const TF* const restrict obuk,
//            const TF* const restrict var,
//            const TF* const restrict z0h,
//            const TF zsl, const Boundary_type bcbot,
//            const int istart, const int iend,
//            const int jstart, const int jend, const int kstart,
//            const int icells, const int jcells, const int kk,
//            Boundary_cyclic<TF>& boundary_cyclic)
//    {
//        const int jj = icells;
//
//        // the surface value is known, calculate the flux and gradient
//        if (bcbot == Boundary_type::Dirichlet_type)
//        {
//            for (int j=0; j<jcells; ++j)
//                #pragma ivdep
//                for (int i=0; i<icells; ++i)
//                {
//                    const int ij  = i + j*jj;
//                    const int ijk = i + j*jj + kstart*kk;
//                    varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*most::fh(zsl, z0h[ij], obuk[ij]);
//                    // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
//                    // use the linearly interpolated grad, rather than the MO grad,
//                    // to prevent giving unresolvable gradients to advection schemes
//                    vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
//                }
//        }
//        else if (bcbot == Boundary_type::Flux_type)
//        {
//            // the flux is known, calculate the surface value and gradient
//            for (int j=0; j<jcells; ++j)
//                #pragma ivdep
//                for (int i=0; i<icells; ++i)
//                {
//                    const int ij  = i + j*jj;
//                    const int ijk = i + j*jj + kstart*kk;
//                    varbot[ij] = varfluxbot[ij] / (ustar[ij]*most::fh(zsl, z0h[ij], obuk[ij])) + var[ijk];
//                    // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
//                    // use the linearly interpolated grad, rather than the MO grad,
//                    // to prevent giving unresolvable gradients to advection schemes
//                    vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
//                }
//        }
//    }
}

template<typename TF>
Boundary_surface_tiled<TF>::Boundary_surface_tiled(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        Boundary<TF>(masterin, gridin, fieldsin, inputin)
{
    swboundary = "surface_tiled";

    // Create the land-surface tiles
    mo_tiles.emplace("veg",  MO_surface_tile<TF>{});
    mo_tiles.emplace("soil", MO_surface_tile<TF>{});
    mo_tiles.emplace("wet",  MO_surface_tile<TF>{});
}

template<typename TF>
Boundary_surface_tiled<TF>::~Boundary_surface_tiled()
{
}

template<typename TF>
void Boundary_surface_tiled<TF>::create(
        Input& input, Netcdf_handle& input_nc,
        Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross)
{
}

template<typename TF>
void Boundary_surface_tiled<TF>::init(Input& inputin, Thermo<TF>& thermo)
{
    // 1. Process the boundary conditions now all fields are registered.
    process_bcs(inputin);

    // 2. Check settings (surface_tiled is really only meant to be used with the LSM)
    check_settings();

    // 3. Allocate and initialize the 2D surface fields.
    init_surface(inputin);

    // 4. Initialize the boundary cyclic.
    boundary_cyclic.init();
}

template<typename TF>
void Boundary_surface_tiled<TF>::check_settings()
{
    // Checks on input
    if (mbcbot != Boundary_type::Dirichlet_type)
    {
        std::string error = "Boundary_surface_tiled requires mbcbot=noslip";
        throw std::runtime_error(error);
    }

    for (auto& it : fields.sp)
    {
        if(sbc.at(it.first).bcbot != Boundary_type::Dirichlet_type)
        {
            std::string error = "Boundary_surface_tiled requires sbcbot=dirichlet";
            throw std::runtime_error(error);
        }
    }
}

template<typename TF>
void Boundary_surface_tiled<TF>::init_surface(Input& input)
{
    auto& gd = grid.get_grid_data();

    // Allocate the surface tiles
    for (auto& tile : mo_tiles)
        init_tile(tile.second, gd.ijcells);

    // Grid point mean quantities
    obuk.resize(gd.ijcells);
    ustar.resize(gd.ijcells);
    z0m.resize(gd.ijcells);

    // Switch between constant (from .ini) or spatial varying z0's
    sw_constant_z0 = input.get_item<TF>("boundary", "sw_constant_z0", "", true);

    // Init roughness lengths
    if (sw_constant_z0)
        for (auto& tile : mo_tiles)
        {
            const TF z0m = input.get_item<TF>("boundary", "z0m", tile.first);
            const TF z0h = input.get_item<TF>("boundary", "z0h", tile.first);

            std::fill(tile.second.z0m.begin(), tile.second.z0m.end(), z0m);
            std::fill(tile.second.z0h.begin(), tile.second.z0h.end(), z0m);
        }

    // Initialize the obukhov lengths on a small number.
    for (auto& tile : mo_tiles)
        std::fill(tile.second.obuk.begin(), tile.second.obuk.end(), Constants::dsmall);
}

template<typename TF>
void Boundary_surface_tiled<TF>::load(const int iotime)
{
}

template<typename TF>
void Boundary_surface_tiled<TF>::save(const int iotime)
{
}

template<typename TF>
void Boundary_surface_tiled<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
}

template<typename TF>
void Boundary_surface_tiled<TF>::exec_stats(Stats<TF>& stats)
{
}

#ifndef USECUDA
template<typename TF>
void Boundary_surface_tiled<TF>::exec_column(Column<TF>& column)
{
}
#endif

template<typename TF>
void Boundary_surface_tiled<TF>::set_values()
{
}

#ifndef USECUDA
template<typename TF>
void Boundary_surface_tiled<TF>::calc_mo_stability(
        Thermo<TF>& thermo, Land_surface<TF>& lsm)
{
    auto& gd = grid.get_grid_data();

    auto buoy  = fields.get_tmp();
    auto dutot = fields.get_tmp();

    // Get buoyancy at surface and first model level
    thermo.get_buoyancy_surf(*buoy, false);
    const TF db_ref = thermo.get_db_ref();

    // Get tile information from LSM
    Tile_map<TF>& lsm_tiles = lsm.get_tiles();

    // Calculate absolute wind speed difference with
    // surface, limited at some minimum value.
    calc_dutot(
            dutot->fld_bot.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("u")->fld_bot.data(),
            fields.mp.at("v")->fld_bot.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart,
            gd.icells, gd.jcells, gd.ijcells);

    boundary_cyclic.exec_2d(dutot->fld_bot.data());

    for (auto& tile : mo_tiles)
    {
        // Calculate surface buoyancy tile
        thermo.get_buoyancy_surf(
                buoy->fld_bot,
                lsm_tiles.at(tile.first).thl_bot,
                lsm_tiles.at(tile.first).qt_bot);

        // Calculate ustar and Obukhov length
        stability(
                tile.second.ustar.data(),
                tile.second.obuk.data(),
                dutot->fld_bot.data(),
                buoy->fld.data(),
                buoy->fld_bot.data(),
                tile.second.z0m.data(),
                tile.second.z0h.data(),
                db_ref, gd.z[gd.kstart],
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart,
                gd.icells, gd.jcells, gd.ijcells);
    }

    // Calculate tile fraction averaged friction velocity
    calc_tiled_mean(
            ustar.data(),
            mo_tiles.at("veg").ustar.data(),
            mo_tiles.at("soil").ustar.data(),
            mo_tiles.at("wet").ustar.data(),
            lsm_tiles.at("veg").fraction.data(),
            lsm_tiles.at("soil").fraction.data(),
            lsm_tiles.at("wet").fraction.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells, gd.jcells);

    // Calculate Obukhov length from mean buoyancy
    // flux and mean friction velocity
    calc_bulk_obuk(
            obuk.data(),
            buoy->flux_bot.data(),
            ustar.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    boundary_cyclic.exec_2d(ustar.data());
    boundary_cyclic.exec_2d(obuk.data());

    fields.release_tmp(buoy);
    fields.release_tmp(dutot);
}

template<typename TF>
void Boundary_surface_tiled<TF>::calc_mo_bcs_momentum(
        Thermo<TF>& thermo, Land_surface<TF>& lsm)
{
    auto& gd = grid.get_grid_data();

    // Calculate the surface fluxes and gradients
    surfm(fields.mp.at("u")->flux_bot.data(),
          fields.mp.at("v")->flux_bot.data(),
          fields.mp.at("u")->grad_bot.data(),
          fields.mp.at("v")->grad_bot.data(),
          fields.mp.at("u")->fld.data(),
          fields.mp.at("v")->fld.data(),
          fields.mp.at("u")->fld_bot.data(),
          fields.mp.at("v")->fld_bot.data(),
          ustar.data(),
          gd.z[gd.kstart],
          gd.istart, gd.iend,
          gd.jstart, gd.jend,
          gd.kstart, gd.icells,
          gd.jcells, gd.ijcells);

    boundary_cyclic.exec_2d(fields.mp.at("u")->flux_bot.data());
    boundary_cyclic.exec_2d(fields.mp.at("v")->flux_bot.data());
}

template<typename TF>
void Boundary_surface_tiled<TF>::calc_mo_bcs_scalars(
        Thermo<TF>& thermo, Land_surface<TF>& lsm)
{
    auto& gd = grid.get_grid_data();
    auto tmp  = fields.get_tmp();

    // Calculate "effective" aerodynamic resistance
    calc_effective_ra(
            tmp->fld_bot.data(),
            fields.sp.at("thl")->fld.data(),
            fields.sp.at("qt")->fld.data(),
            fields.sp.at("thl")->fld_bot.data(),
            fields.sp.at("qt")->fld_bot.data(),
            fields.sp.at("thl")->flux_bot.data(),
            fields.sp.at("qt")->flux_bot.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart,
            gd.icells, gd.ijcells);

    throw 1;

    //for (auto& it : fields.sp)
    //    surfs(it.second->fld_bot.data(),
    //          it.second->grad_bot.data(),
    //          it.second->flux_bot.data(),
    //          ustar.data(), obuk.data(),
    //          it.second->fld.data(), z0h.data(),
    //          gd.z[gd.kstart], sbc.at(it.first).bcbot,
    //          gd.istart, gd.iend,
    //          gd.jstart, gd.jend, gd.kstart,
    //          gd.icells, gd.jcells, gd.ijcells,
    //          boundary_cyclic);

    fields.release_tmp(tmp);
}

template<typename TF>
void Boundary_surface_tiled<TF>::get_ra(Field3d<TF>& fld)
{
}

template<typename TF>
void Boundary_surface_tiled<TF>::get_ra(Field3d<TF>& fld, std::string tile)
{
    auto& gd = grid.get_grid_data();

    calc_ra(fld.flux_bot.data(),
            mo_tiles.at(tile).ustar.data(),
            mo_tiles.at(tile).obuk.data(),
            mo_tiles.at(tile).z0h.data(),
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
}

template<typename TF>
void Boundary_surface_tiled<TF>::get_duvdz(
        std::vector<TF>& dudz, std::vector<TF>& dvdz)
{
}

template<typename TF>
void Boundary_surface_tiled<TF>::get_dbdz(
        std::vector<TF>& dbdz, std::vector<TF>& bfluxbot)
{
}
#endif

template<typename TF>
void Boundary_surface_tiled<TF>::update_slave_bcs()
{
    // This function does nothing when the surface model is enabled, because
    // the fields are computed by the surface model in update_bcs.
}

template class Boundary_surface_tiled<double>;
template class Boundary_surface_tiled<float>;
