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
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <netcdf>       // C++
#include <netcdf.h>     // C, for sync() using older netCDF-C++ versions
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
//#include "thermo_moist.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
//#include "diff_smag2.h"
#include "timeloop.h"
#include "diff.h"

using namespace netCDF;
using namespace netCDF::exceptions;
using namespace Constants;

namespace
{
    // Help functions to switch between the different NetCDF data types
    template<typename TF> NcType netcdf_fp_type();
    template<> NcType netcdf_fp_type<double>() { return ncDouble; }
    template<> NcType netcdf_fp_type<float>()  { return ncFloat; }

    template<typename TF> TF netcdf_fp_fillvalue();
    template<> double netcdf_fp_fillvalue<double>() { return NC_FILL_DOUBLE; }
    template<> float  netcdf_fp_fillvalue<float>()  { return NC_FILL_FLOAT; }


    template<typename TF, Stats_mask_type mode>
    TF is_false(const TF value, const TF threshold)
    {
        if (mode == Stats_mask_type::Plus)
            return (value <= threshold);
        else if (mode == Stats_mask_type::Min)
            return (value > threshold);
    }

    template<typename TF, Stats_mask_type mode>
    void calc_mask_thres(unsigned int* const restrict mfield, unsigned int* const restrict mfield_bot, const unsigned int flag, const unsigned int flagh,
                        const TF* const restrict fld, const TF* const restrict fldh,
                        const TF* const restrict fld_bot, const TF threshold,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend,
                        const int icells, const int ijcells)
    {

        #pragma omp parallel for
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    mfield[ijk] -= (mfield[ijk] & flag)  * is_false<TF, mode>(fld[ijk], threshold);
                    mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk], threshold);
                }

        // Set the mask for surface projected quantities
        #pragma omp parallel for
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*icells;
                // const int ij2 = i + j*icells + (kend+1)*ijcells;
                mfield_bot[ij] -= (mfield_bot[ij] & flag)  * is_false<TF, mode>(fld_bot[ij], threshold);
                // mfield[ij2]    -= (mfield[ij2] & flagh) * is_false<TF, mode>(fldh[ij2], threshold);
            }
    }


    template<typename TF, Stats_mask_type mode>
    void calc_mask_thres_pert(unsigned int* const restrict mfield, unsigned int* const restrict mfield_bot, const unsigned int flag, const unsigned int flagh,
                        const TF* const restrict fld, const TF* const restrict fld_mean, const TF* const restrict fldh,
                        const TF* const restrict fldh_mean, const TF* const restrict fld_bot, const TF threshold,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend,
                        const int icells, const int ijcells)
    {

        #pragma omp parallel for
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    mfield[ijk] -=  (mfield[ijk] & flag)  * is_false<TF, mode>(fld[ijk]-fld_mean[k], threshold);
                    mfield[ijk] -=  (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk]-fldh_mean[k], threshold);
                }
        // Set the mask for surface projected quantities
        #pragma omp parallel for
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*icells;
                // const int ij2 = i + j*icells + (kend+1)*ijcells;
                mfield_bot[ij] -= (mfield_bot[ij] & flag)  * is_false<TF, mode>(fld_bot[ij]-fld_mean[kstart], threshold);
                // mfield[ij2]    -= (mfield[ij2] & flagh) * is_false<TF, mode>(fldh[ij2]-fld_mean[kend+1], threshold);
            }

    }

    template<typename TF>
    void calc_area(TF* const restrict area, const int loc[3], const int* const restrict nmask, const int kstart, const int kend, const int ijtot)
    {
        for (int k=kstart; k<kend+loc[2]; k++)
        {
            if (nmask[k])
                area[k] = static_cast<TF>(nmask[k]) / static_cast<TF>(ijtot);
            else
                area[k] = 0.;
        }
    }

    // Sets all the mask values to one (non-masked field)
    template<typename TF>
    void calc_nmask(int* restrict nmask_full, int* restrict nmask_half, int& nmask_bottom,
                    const unsigned int* const mfield, const unsigned int* const mfield_bot,const unsigned int flag,const unsigned int flagh,
                   const int istart, const int iend, const int jstart, const int jend,
                   const int kstart, const int kend, const int icells, const int ijcells, const int kcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            nmask_full[k] = 0;
            nmask_half[k] = 0;

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk  = i + j*icells + k*ijcells;
                    nmask_full[k]+=((mfield[ijk] & flag)>0);
                    nmask_half[k]+=((mfield[ijk] & flagh)>0);
                }
        }
        nmask_bottom = 0;
        for (int j=jstart; j<jend; ++j)
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                nmask_bottom+=((mfield_bot[ij] & flag)>0);
            }
    }

        template<typename TF>
        void calc_mean(TF* const restrict prof, const TF* const restrict fld, const TF offset,
                        const unsigned int* const mask, const unsigned int flag, const int* const nmask,
                        const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
        {
            #pragma omp parallel for
            for (int k=kstart; k<kend+1; k++)
            {
                if (nmask[k])
                {
                    prof[k] = 0.;
                    for (int j=jstart; j<jend; j++)
                        #pragma ivdep
                        for (int i=istart; i<iend; i++)
                        {
                            const int ijk  = i + j*icells + k*ijcells;
                            prof[k] += static_cast<TF>((mask[ijk] & flag)>0)*(fld[ijk] + offset);
                        }
                    prof[k] /= static_cast<TF>(nmask[k]);
                }
                else
                    prof[k] = netcdf_fp_fillvalue<TF>();
            }
        }


        template<typename TF>
        void calc_mean_2d(TF& out, const TF* const restrict fld, const TF offset,
                        const int istart, const int iend, const int jstart, const int jend, const int icells, const int itot, const int jtot)
        {
                    out = 0.;
                    for (int j=jstart; j<jend; j++)
                        #pragma ivdep
                        for (int i=istart; i<iend; i++)
                        {
                            const int ij  = i + j*icells;
                            out = fld[ij] + offset;
                        }
                    out /= static_cast<TF>(itot*jtot);
        }

    template<typename TF>
    void calc_moment(TF* const restrict prof, const TF* const restrict fld, const TF* const restrict fld_mean, const TF offset,
                    const unsigned int* const mask, const unsigned int flag, const int* const nmask, const int power,
                    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; k++)
        {
            if (nmask[k])
            {
                prof[k] = 0.;
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk  = i + j*icells + k*ijcells;
                        prof[k] += static_cast<TF>(mask[ijk] & flag)*std::pow(fld[ijk] - fld_mean[k] + offset, power);
                    }
                prof[k] /= static_cast<TF>(nmask[k]);
            }
            else
                prof[k] = netcdf_fp_fillvalue<TF>();
        }
    }

    template<typename TF>
    void add_fluxes(TF* const restrict flux, const TF* const restrict turb, const TF* const restrict diff, const TF fillvalue, const int kstart, const int kend)
    {
        for (int k=kstart; k<kend+1; ++k)
        {
            if (turb[k] == fillvalue || diff[k] == fillvalue)
                flux[k] = fillvalue;
            else
                flux[k] = turb[k] + diff[k];
        }
    }


    template<typename TF>
    void calc_frac(TF* const restrict prof, const TF* const restrict fld, const TF offset, const TF threshold,
                    const unsigned int* const mask, const unsigned int flag, const int* const nmask,
                    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; k++)
        {
            if (nmask[k])
            {
                prof[k] = 0.;
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk  = i + j*icells + k*ijcells;
                        prof[k] += static_cast<TF>((mask[ijk] & flag) > 0)*(fld[ijk] + offset > threshold);
                    }
                prof[k] /= static_cast<TF>(nmask[k]);
            }
            else
                prof[k] = netcdf_fp_fillvalue<TF>();
        }
    }

    template<typename TF>
    void calc_path(TF& path, const TF* const restrict data, const TF* const restrict dz, const TF* const restrict rho, const int* const restrict nmask, const int kstart, const int kend)
    {
        path = 0.;
        for (int k=kstart; k<kend; k++)
        {
            if(nmask[k])
            {
                path += data[k]*rho[k]*dz[k];
            }
        }
        path /= static_cast<TF>(nmask[kstart]);
    }

    template<typename TF>
    void calc_cover(TF& cover, const TF* const restrict fld, const TF offset, const TF threshold,
                    const unsigned int* const mask, const unsigned int flag, const int* const nmask,
                    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
    {
        cover  = 0.;
        TF nmaskcover = 0.;
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                for (int k=kstart; k<kend+1; k++)
                {
                    const int ijk  = i + j*icells + k*ijcells;
                    if((mask[ijk] & flag) > 0)
                    {
                        ++nmaskcover;
                        if((fld[ijk] + offset > threshold))
                        {
                            ++cover;
                            break;
                        }
                    }
                }

            }
        cover /= nmaskcover;
    }


    bool has_only_digits(const std::string s)
    {
        return s.find_first_not_of( "23456789" ) == std::string::npos;
    }
}

template<typename TF>
Stats<TF>::Stats(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin):
    master(masterin), grid(gridin), fields(fieldsin),  boundary_cyclic(master, grid)

{
    swstats = inputin.get_item<bool>("stats", "swstats", "", false);

    if (swstats)
    {
        sampletime = inputin.get_item<double>("stats", "sampletime", "");
        masklist   = inputin.get_list<std::string>("stats", "masklist", "", std::vector<std::string>());
        masklist.push_back("default");  // Add the default mask, which calculates the domain mean without sampling.
        std::vector<std::string> whitelistin = inputin.get_list<std::string>("stats", "whitelist", "", std::vector<std::string>());

        for (auto& it : whitelistin)
        {
            std::regex re(it);
            whitelist.push_back(re);
        }
        std::vector<std::string> blacklistin = inputin.get_list<std::string>("stats", "blacklist", "", std::vector<std::string>());
        for (auto& it : blacklistin)
        {
            std::regex re(it);
            blacklist.push_back(re);
        }

    }
}

template<typename TF>
Stats<TF>::~Stats()
{
}

template<typename TF>
void Stats<TF>::init(double ifactor)
{
    if (!swstats)
        return;

    auto& gd = grid.get_grid_data();

    boundary_cyclic.init();

    isampletime = static_cast<unsigned long>(ifactor * sampletime);
    statistics_counter = 0;

    // Vectors which hold the amount of grid points sampled on each model level.
    mfield.resize(gd.ncells);
    mfield_bot.resize(gd.ijcells);
}

template<typename TF>
void Stats<TF>::create(int iotime, std::string sim_name)
{
    // Do not create statistics file if stats is disabled.
    if (!swstats)
        return;

    int nerror = 0;
    auto& gd = grid.get_grid_data();

    // Create a NetCDF file for each of the masks.
    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        if (master.get_mpiid() == 0)
        {
            std::stringstream filename;
            filename << sim_name << "." << m.name << "." << std::setfill('0') << std::setw(7) << iotime << ".nc";

            // Create new NetCDF file, and catch any exceptions locally, to be able
            // to communicate them to the other processes.
            try
            {
                m.data_file = new NcFile(filename.str(), NcFile::newFile);
            }
            catch(NcException& e)
            {
                master.print_error("NetCDF exception: %s\n",e.what());
                ++nerror;
            }
        }

        // Crash on all processes in case the file could not be written.
        master.broadcast(&nerror, 1);
        if (nerror)
            throw 1;

        // Create dimensions.
        if (master.get_mpiid() == 0)
        {
            m.z_dim  = m.data_file->addDim("z" , gd.kmax);
            m.zh_dim = m.data_file->addDim("zh", gd.kmax+1);
            m.t_dim  = m.data_file->addDim("t");

            NcVar z_var;
            NcVar zh_var;

            // Create variables belonging to dimensions.
            m.iter_var = m.data_file->addVar("iter", ncInt, m.t_dim);
            m.iter_var.putAtt("units", "-");
            m.iter_var.putAtt("long_name", "Iteration number");

            m.t_var = m.data_file->addVar("time", ncDouble, m.t_dim);
            m.t_var.putAtt("units", "s");
            m.t_var.putAtt("long_name", "Time");

            z_var = m.data_file->addVar("z", netcdf_fp_type<TF>(), m.z_dim);
            z_var.putAtt("units", "m");
            z_var.putAtt("long_name", "Full level height");

            zh_var = m.data_file->addVar("zh", netcdf_fp_type<TF>(), m.zh_dim);
            zh_var.putAtt("units", "m");
            zh_var.putAtt("long_name", "Half level height");

            // Save the grid variables.
            z_var .putVar(&gd.z [gd.kstart]);
            zh_var.putVar(&gd.zh[gd.kstart]);

            // Synchronize the NetCDF file.
            // BvS: only the last netCDF4-c++ includes the NcFile->sync()
            //      for now use sync() from the netCDF-C library to support older NetCDF4-c++ versions
            //m.data_file->sync();
            nc_sync(m.data_file->getId());
        }

        m.nmask. resize(gd.kcells);
        m.nmaskh.resize(gd.kcells);
    }

    // For each mask, add the area as a variable.
    add_prof("area" , "Fractional area contained in mask", "-", "z" );
    add_prof("areah", "Fractional area contained in mask", "-", "zh");

}

template<typename TF>
unsigned long Stats<TF>::get_time_limit(unsigned long itime)
{
    // If statistics is disabled, return large (... huge!) value.
    if (!swstats)
        return Constants::ulhuge;

    unsigned long idtlim = isampletime - itime % isampletime;
    return idtlim;
}

template<typename TF>
bool Stats<TF>::do_statistics(unsigned long itime)
{
    // Check if stats are enabled.
    if (!swstats)
        return false;

    // Check if time for execution.
    if (itime % isampletime != 0)
        return false;

    // Return true such that stats are computed.
    return true;
}

template<typename TF>
void Stats<TF>::exec(int iteration, double time, unsigned long itime)
{
    if (!swstats)
        return;

    auto& gd = grid.get_grid_data();

    // check if time for execution
    // BvS: why was this used? This function is only called after a stats->do_statistics(), which already checks the sampletime...
    //if (itime % isampletime != 0)
    //    return;

    // Write message in case stats is triggered
    master.print_message("Saving statistics for time %f\n", time);

    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        // Put the data into the NetCDF file
        if (master.get_mpiid() == 0)
        {
            const std::vector<size_t> time_index = {static_cast<size_t>(statistics_counter)};

            // Write the time and iteration number
            m.t_var   .putVar(time_index, &time     );
            m.iter_var.putVar(time_index, &iteration);

            const std::vector<size_t> time_height_index = {static_cast<size_t>(statistics_counter), 0};
            std::vector<size_t> time_height_size  = {1, 0};

            //for (Prof_map::const_iterator it=m.profs.begin(); it!=m.profs.end(); ++it)
            for (auto& p : m.profs)
            {
                time_height_size[1] = m.profs[p.first].ncvar.getDim(1).getSize();
                m.profs[p.first].ncvar.putVar(time_height_index, time_height_size, &m.profs[p.first].data.data()[gd.kstart]);
            }

            for (auto& ts: m.tseries)
                m.tseries[ts.first].ncvar.putVar(time_index, &m.tseries[ts.first].data);

            // Synchronize the NetCDF file
            // BvS: only the last netCDF4-c++ includes the NcFile->sync()
            //      for now use sync() from the netCDF-C library to support older NetCDF4-c++ versions
            //m.dataFile->sync();
            nc_sync(m.data_file->getId());
        }
    }
    wmean_set = false;
    // Increment the statistics index
    ++statistics_counter;
}

// Retrieve the user input list of requested masks
template<typename TF>
const std::vector<std::string>& Stats<TF>::get_mask_list()
{
    return masklist;
}

// Add a new mask to the mask map
template<typename TF>
void Stats<TF>::add_mask(const std::string maskname)
{
    auto& gd = grid.get_grid_data();

    masks[maskname].name = maskname;
    masks[maskname].data_file = 0;
    int nmasks = masks.size();
    masks[maskname].flag = (1 << (2 * (nmasks - 1)));
    masks[maskname].flagh = (1 << (2 * (nmasks-1) + 1));
}

// Add a new profile to each of the NetCDF files
template<typename TF>
void Stats<TF>::add_prof(std::string name, std::string longname, std::string unit, std::string zloc)
{
    auto& gd = grid.get_grid_data();

    //Check whether variable is part of whitelist/blacklist;
    if(is_blacklisted(name))
        return;

    // Add profile to all the NetCDF files
    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        // Create the NetCDF variable
        if (master.get_mpiid() == 0)
        {
            std::vector<NcDim> dim_vector = {m.t_dim};

            if (zloc == "z")
            {
                dim_vector.push_back(m.z_dim);
                m.profs[name].ncvar = m.data_file->addVar(name, netcdf_fp_type<TF>(), dim_vector);
                //m.profs[name].data = NULL;
            }
            else if (zloc == "zh")
            {
                dim_vector.push_back(m.zh_dim);
                m.profs[name].ncvar = m.data_file->addVar(name.c_str(), netcdf_fp_type<TF>(), dim_vector);
                //m.profs[name].data = NULL;
            }
            m.profs[name].ncvar.putAtt("units", unit.c_str());
            m.profs[name].ncvar.putAtt("long_name", longname.c_str());
            m.profs[name].ncvar.putAtt("_FillValue", netcdf_fp_type<TF>(), netcdf_fp_fillvalue<TF>());

            nc_sync(m.data_file->getId());
        }

        // Resize the vector holding the data at all processes
        m.profs[name].data.resize(gd.kcells);

        varlist.push_back(name);
    }
}

template<typename TF>
void Stats<TF>::add_fixed_prof(std::string name, std::string longname, std::string unit, std::string zloc, TF* restrict prof)
{
    auto& gd = grid.get_grid_data();

    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        // Create the NetCDF variable
        if (master.get_mpiid() == 0)
        {
           NcVar var;
           if (zloc == "z")
               var = m.data_file->addVar(name, netcdf_fp_type<TF>(), m.z_dim);
           else if (zloc == "zh")
               var = m.data_file->addVar(name, netcdf_fp_type<TF>(), m.zh_dim);
           var.putAtt("units", unit.c_str());
           var.putAtt("long_name", longname.c_str());
           var.putAtt("_FillValue", netcdf_fp_type<TF>(), netcdf_fp_fillvalue<TF>());

           const std::vector<size_t> index = {0};
           if (zloc == "z")
           {
               const std::vector<size_t> size  = {static_cast<size_t>(gd.kmax)};
               var.putVar(index, size, &prof[gd.kstart]);
           }
           else if (zloc == "zh")
           {
               const std::vector<size_t> size  = {static_cast<size_t>(gd.kmax+1)};
               var.putVar(index, size, &prof[gd.kstart]);
           }
       }
   }
}
//
template<typename TF>
void Stats<TF>::add_time_series(const std::string name, const std::string longname, const std::string unit)
{
    //Check whether variable is part of whitelist/blacklist;
    if(is_blacklisted(name))
        return;

    // add the series to all files
    for (auto& mask : masks)
    {
        // shortcut
        Mask<TF>& m = mask.second;

        // create the NetCDF variable
        if (master.get_mpiid() == 0)
        {
            m.tseries[name].ncvar = m.data_file->addVar(name.c_str(), netcdf_fp_type<TF>(), m.t_dim);
            m.tseries[name].ncvar.putAtt("units", unit.c_str());
            m.tseries[name].ncvar.putAtt("long_name", longname.c_str());
            m.tseries[name].ncvar.putAtt("_FillValue", netcdf_fp_type<TF>(), netcdf_fp_fillvalue<TF>());
        }

        // Initialize at zero
        m.tseries[name].data = 0.;
    }
    varlist.push_back(name);
}

template<typename TF>
bool Stats<TF>::is_blacklisted(const std::string name)
{
    for (const auto& it : whitelist)
    {
        if(std::regex_match(name, it))
            return false;
    }

    for (const auto& it : blacklist)
    {
        if(std::regex_match(name, it))
        {
            master.print_message("NOT doing statistics for  %s\n", name.c_str());
            return true;
        }
    }
    return false;
}

template<typename TF>
void Stats<TF>::initialize_masks()
{
    auto& gd = grid.get_grid_data();
    unsigned int flagmax = 0;
    for (auto& it : masks)
    {
        flagmax+=it.second.flag+it.second.flagh;
    }
    for (int n=0; n<gd.ncells; ++n)
        mfield[n] = flagmax;
    for (int n=0; n<gd.ijcells; ++n)
        mfield_bot[n] = flagmax;
}


template<typename TF>
void Stats<TF>::finalize_masks()
{
    auto& gd = grid.get_grid_data();
    const int sloc[] = {0,0,0};
    const int wloc[] = {0,0,1};

    boundary_cyclic.exec(mfield.data());
    boundary_cyclic.exec_2d(mfield_bot.data());
    for (auto& it : masks)
    {
        calc_nmask<TF>(it.second.nmask.data(), it.second.nmaskh.data(), it.second.nmask_bot,
                       mfield.data(), mfield_bot.data(), it.second.flag, it.second.flagh,
                       gd.istart, gd.iend, gd.jstart, gd.jend, 0, gd.kcells,
                       gd.icells, gd.ijcells, gd.kcells);
        master.sum(it.second.nmask.data() , gd.kcells);
        master.sum(it.second.nmaskh.data(), gd.kcells);
        it.second.nmask_bot = it.second.nmaskh[gd.kstart];
        calc_area(it.second.profs["area" ].data.data(), sloc, it.second.nmask .data(), gd.kstart, gd.kend, gd.itot*gd.jtot);
        calc_area(it.second.profs["areah"].data.data(), wloc, it.second.nmaskh.data(), gd.kstart, gd.kend, gd.itot*gd.jtot);
    }
}

template<typename TF>
void Stats<TF>::set_mask_thres(std::string mask_name, Field3d<TF>& fld, Field3d<TF>& fldh, TF threshold, Stats_mask_type mode)
{
    auto& gd = grid.get_grid_data();
    unsigned int flag, flagh;
    bool found_mask = false;

    for(auto& it : masks)
    {
        if(it.second.name == mask_name)
        {
            found_mask = true;
            flag = it.second.flag;
            flagh = it.second.flagh;
        }
    }
    if(!found_mask)
        throw std::runtime_error("Invalid mask name in set_mask_thres()");

    if (mode == Stats_mask_type::Plus)
        calc_mask_thres<TF, Stats_mask_type::Plus>(mfield.data(), mfield_bot.data(), flag, flagh,
            fld.fld.data(), fldh.fld.data(), fldh.fld_bot.data(), threshold,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
    else if (mode == Stats_mask_type::Min)
        calc_mask_thres<TF, Stats_mask_type::Min>(mfield.data(), mfield_bot.data(), flag, flagh,
            fld.fld.data(), fldh.fld.data(), fldh.fld_bot.data(), threshold,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
    else
        throw std::runtime_error("Invalid mask type in set_mask_thres()");
}

template<typename TF>
void Stats<TF>::set_mask_thres_pert(std::string mask_name, Field3d<TF>& fld, Field3d<TF>& fldh, TF threshold, Stats_mask_type mode)
{

    auto& gd = grid.get_grid_data();
    unsigned int flag, flagh;
    bool found_mask = false;

    for(auto& it : masks)
    {
        if(it.second.name == mask_name)
        {
            found_mask = true;
            flag  = it.second.flag;
            flagh = it.second.flagh;
        }
    }
    if(!found_mask)
        throw std::runtime_error("Invalid mask name in set_mask_thres()");

    if (mode == Stats_mask_type::Plus)
        calc_mask_thres_pert<TF, Stats_mask_type::Plus>(mfield.data(), mfield_bot.data(), flag, flagh,
            fld.fld.data(), fld.fld_mean.data(), fldh.fld.data(), fldh.fld_mean.data(), fldh.fld_bot.data(), threshold,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
    else if (mode == Stats_mask_type::Min)
        calc_mask_thres_pert<TF, Stats_mask_type::Min>(mfield.data(), mfield_bot.data(), flag, flagh,
            fld.fld.data(), fld.fld_mean.data(), fldh.fld.data(), fldh.fld_mean.data(), fldh.fld_bot.data(), threshold,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
    else
        throw std::runtime_error("Invalid mask type in set_mask_thres_pert()");

}


template<typename TF>
void Stats<TF>::set_prof(const std::string varname, const std::vector<TF> prof)
{
    for (auto& it : masks)
    {
        it.second.profs.at(varname).data = prof;
    }
}

template<typename TF>
void Stats<TF>::calc_stats(const std::string varname, const Field3d<TF>& fld, const int* loc, const TF offset, const TF threshold, std::vector<std::string> operations, Diff<TF>& diff)
{
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    std::string name;

    sanatize_operations_vector(operations);
    // Process mean first
    auto it = std::find(operations.begin(), operations.end(), "mean");
    if (it != operations.end())
    {
        auto it1 = std::find(varlist.begin(), varlist.end(), varname);
        if (it1 != varlist.end())
        {
            for (auto& m : masks)
            {
                if(loc[2]==0)
                    flag = m.second.flag;
                else
                    flag = m.second.flagh;

                calc_mean(m.second.profs.at(varname).data.data(), fld.fld.data(), offset, mfield.data(), flag, m.second.nmask.data(),
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
                master.sum(m.second.profs.at(varname).data.data(), gd.kcells);
            }
            if(varname == "w")
                wmean_set = true;
            operations.erase(it);
        }
    }

    //Loop over all other operations.
    for(auto& it : operations)
    {
        name = varname+it;
        auto it1 = std::find(varlist.begin(), varlist.end(), name);
        if (it1 == varlist.end())
        {
        }
        else if(has_only_digits(it))
        {
            int power = std::stoi(it);
            for (auto& m : masks)
            {

                if(loc[2]==0)
                    flag = m.second.flag;
                else
                    flag = m.second.flagh;
                calc_moment(m.second.profs.at(name).data.data(), fld.fld.data(), m.second.profs.at(varname).data.data(), offset, mfield.data(), flag, m.second.nmask.data(),
                        power, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);

                master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            }

        }
        else if (it == "w")
        {
            if(!wmean_set)
                throw std::runtime_error("W mean not calculated in stat - needed for flux");

            auto tmp = fields.get_tmp();
            for (auto& m : masks)
            {
                if(loc[2]==0)
                    flag = m.second.flagh;
                else
                    flag = m.second.flag;
                if (grid.get_spatial_order() == Grid_order::Second)
                {
                    calc_flux_2nd(m.second.profs.at(name).data.data(), fld.fld.data(), m.second.profs.at(varname).data.data(), fields.mp["w"]->fld.data(), m.second.profs.at("w").data.data(),
                        tmp->fld.data(), loc, mfield.data(), flag, m.second.nmask.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
                }
                else if (grid.get_spatial_order() == Grid_order::Fourth)
                {
                    calc_flux_2nd(m.second.profs.at(name).data.data(), fld.fld.data(), m.second.profs.at(varname).data.data(), fields.mp["w"]->fld.data(), m.second.profs.at("w").data.data(),
                        tmp->fld.data(), loc, mfield.data(), flag, m.second.nmask.data(), gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
                }
            master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            }
            fields.release_tmp(tmp);

        }
        else if (it == "diff")
        {
            auto diffusion = fields.get_tmp();
            diff.diff_flux(*diffusion, fld, loc);

            for (auto& m : masks)
            {
                if(loc[2]==0)
                    flag = m.second.flagh;
                else
                    flag = m.second.flag;

                calc_mean(m.second.profs.at(name).data.data(), diffusion->fld.data(), offset, mfield.data(), flag, m.second.nmask.data(),
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
                master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            }
            fields.release_tmp(diffusion);

        }
        else if (it == "flux")
        {
            for (auto& m : masks)
            {
                add_fluxes(m.second.profs.at(name).data.data(), m.second.profs.at(varname+"w").data.data(), m.second.profs.at(varname+"diff").data.data(),
                        netcdf_fp_fillvalue<TF>(), gd.kstart, gd.kend);
                master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            }

        }
        else if (it == "grad")
        {
            for (auto& m : masks)
            {
                if(loc[2]==0)
                    flag = m.second.flagh;
                else
                    flag = m.second.flag;
                if (grid.get_spatial_order() == Grid_order::Second)
                {
                    calc_grad_2nd(m.second.profs.at(name).data.data(), fld.fld.data(), gd.dzhi.data(), mfield.data(), flag, m.second.nmask.data(),
                            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
                }
                else if (grid.get_spatial_order() == Grid_order::Fourth)
                {
                    calc_grad_4th(m.second.profs.at(name).data.data(), fld.fld.data(), gd.dzhi4.data(), mfield.data(), flag, m.second.nmask.data(),
                            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
                }
                master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            }
        }
        else if (it == "path")
        {
            for (auto& m : masks)
            {
                calc_path(m.second.tseries.at(name).data, m.second.profs.at(varname).data.data(), gd.dz.data(), fields.rhoref.data(), m.second.nmask.data(), gd.kstart, gd.kend);
                master.sum(&m.second.tseries.at(name).data, 1);
            }
        }
        else if (it == "cover")
        {
            for (auto& m : masks)
            {
                if(loc[2]==0)
                    flag = m.second.flag;
                else
                    flag = m.second.flagh;

                calc_cover(m.second.tseries.at(name).data,fld.fld.data(), offset, threshold, mfield.data(), flag, m.second.nmask.data(),
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
                master.sum(&m.second.tseries.at(name).data, 1);
            }
        }
        else if (it == "frac")
        {
            for (auto& m : masks)
            {
                if(loc[2]==0)
                    flag = m.second.flag;
                else
                    flag = m.second.flagh;

                calc_frac(m.second.profs.at(name).data.data(), fld.fld.data(), offset, threshold, mfield.data(), flag, m.second.nmask.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
                master.sum(m.second.profs.at(varname).data.data(), gd.kcells);

            }
        }
        else
        {
            throw std::runtime_error("Invalid operations in stat.");
        }
    }
}

template<typename TF>
void Stats<TF>::calc_stats_2d(const std::string varname, const std::vector<TF>& fld, const TF offset, std::vector<std::string> operations)
{
    auto& gd = grid.get_grid_data();

    // Process mean first
    auto it = std::find(operations.begin(), operations.end(), "mean");
    if (it != operations.end())
    {
        auto it1 = std::find(varlist.begin(), varlist.end(), varname);
        if (it1 != varlist.end())
        {
            for (auto& m : masks)
            {

                calc_mean_2d(m.second.tseries.at(varname).data, fld.data(), offset,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.icells, gd.itot, gd.jtot);
                master.sum(&m.second.tseries.at(varname).data, 1);
            }
        }
    }
}


template<typename TF>
void Stats<TF>::calc_flux_2nd(TF* const restrict prof, const TF* const restrict data, const TF* const restrict fld_mean, TF* const restrict w, const TF* const restrict wmean,
                              TF* restrict tmp1, const int loc[3], const unsigned int* const mask, const unsigned int flag, const int* const nmask,
                              const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
{

    // set a pointer to the field that contains w, either interpolated or the original
    TF* restrict calcw = w;

    // define the locations
    const int wloc [3] = {0,0,1};
    const int uwloc[3] = {1,0,1};
    const int vwloc[3] = {0,1,1};

    if (loc[0] == 1)
    {
        grid.interpolate_2nd(tmp1, w, wloc, uwloc);
        calcw = tmp1;
    }
    else if (loc[1] == 1)
    {
        grid.interpolate_2nd(tmp1, w, wloc, vwloc);
        calcw = tmp1;
    }

    #pragma omp parallel for
    for (int k=kstart; k<kend+1; ++k)
    {
        if (nmask[k] && fld_mean[k-1] != netcdf_fp_fillvalue<TF>() && fld_mean[k] != netcdf_fp_fillvalue<TF>())
        {
            prof[k] = 0.;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk  = i + j*icells + k*ijcells;
                    prof[k] += static_cast<TF>(mask[ijk] & flag)*(0.5*(data[ijk-ijcells]+data[ijk])-0.5*(fld_mean[k-1]+fld_mean[k]))*(calcw[ijk]-wmean[k]);
                }
            prof[k] /= static_cast<TF>(nmask[k]);
        }
        else
            prof[k] = netcdf_fp_fillvalue<TF>();
    }
}


template<typename TF>
void Stats<TF>::calc_flux_4th(TF* const restrict prof, const TF* const restrict data, const TF* const restrict fld_mean, TF* const restrict w, const TF* const restrict wmean,
                              TF* restrict tmp1, const int loc[3], const unsigned int* const mask, const unsigned int flag, const int* const nmask,
                              const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
{

    using namespace Finite_difference::O4;

    auto& gd = grid.get_grid_data();

    const int jj  = 1*icells;
    const int kk1 = 1*ijcells;
    const int kk2 = 2*ijcells;

    // set a pointer to the field that contains w, either interpolated or the original
    TF* restrict calcw = w;

    // define the locations
    const int wloc [3] = {0,0,1};
    const int uwloc[3] = {1,0,1};
    const int vwloc[3] = {0,1,1};

    if (loc[0] == 1)
    {
        grid.interpolate_4th(tmp1, w, wloc, uwloc);
        calcw = tmp1;
    }
    else if (loc[1] == 1)
    {
        grid.interpolate_4th(tmp1, w, wloc, vwloc);
        calcw = tmp1;
    }

    #pragma omp parallel for
    for (int k=kstart; k<kend+1; ++k)
    {
        if (nmask[k] && fld_mean[k-1] != netcdf_fp_fillvalue<TF>() && fld_mean[k] != netcdf_fp_fillvalue<TF>())
        {
            prof[k] = 0.;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk  = i + j*icells + k*ijcells;
                    prof[k] += static_cast<TF>(mask[ijk] & flag)*(ci0<TF>*data[ijk-kk2] + ci1<TF>*data[ijk-kk1] + ci2<TF>*data[ijk] + ci3<TF>*data[ijk+kk1])*calcw[ijk];
                }
            prof[k] /= static_cast<TF>(nmask[k]);
        }
        else
            prof[k] = netcdf_fp_fillvalue<TF>();
    }
}


template<typename TF>
void Stats<TF>::calc_grad_2nd(TF* const restrict prof, const TF* const restrict data, const TF* const restrict dzhi,
                    const unsigned int* const mask, const unsigned int flag, const int* const nmask,
                    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
{
    #pragma omp parallel for
    for (int k=kstart; k<kend+1; ++k)
    {
        if (nmask[k])
        {
            prof[k] = 0.;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    prof[k] +=static_cast<TF>(mask[ijk] & flag)*(data[ijk]-data[ijk-ijcells])*dzhi[k];
                }
            prof[k] /= static_cast<TF>(nmask[k]);
        }
        else
            prof[k] = netcdf_fp_fillvalue<TF>();

    }

}

template<typename TF>
void Stats<TF>::calc_grad_4th(
        TF* const restrict prof, const TF* const restrict data, const TF* const restrict dzhi4,
        const unsigned int* const mask, const unsigned int flag, const int* const nmask,
        const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
{
    using namespace Finite_difference::O4;

    auto& gd = grid.get_grid_data();

    const int jj  = 1*icells;
    const int kk1 = 1*ijcells;
    const int kk2 = 2*ijcells;

    #pragma omp parallel for
    for (int k=kstart; k<kend+1; ++k)
    {
        if (nmask[k])
        {
            prof[k] = 0.;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk1;
                    prof[k] +=static_cast<TF>(mask[ijk] & flag)*(cg0<TF>*data[ijk-kk2] + cg1<TF>*data[ijk-kk1] + cg2<TF>*data[ijk] + cg3<TF>*data[ijk+kk1])*dzhi4[k];
                }
            prof[k] /= static_cast<TF>(nmask[k]);
        }
        else
            prof[k] = netcdf_fp_fillvalue<TF>();

    }
}


template<typename TF>
void Stats<TF>::sanatize_operations_vector(std::vector<std::string> operations)
{
    // Sanatize the operations vector:
    //find instances that need a mean ({2,3,4,5}); if so, add it to the vector if necessary
    std::vector<std::string> tmpvec = operations;
    for(auto it : tmpvec)
    {
        if(it == "flux")
        {
            operations.push_back("diff");
            operations.push_back("w");
        }
    }
    for(auto it : tmpvec)
    {
        if(has_only_digits(it) || (it == "w") || (it == "path"))
        {
            operations.push_back("mean");
        }
    }
    // Check for duplicates
    std::sort( operations.begin(), operations.end() );
    operations.erase( unique( operations.begin(), operations.end() ), operations.end() );
    // Make sure that flux goes at the end
    for(auto& it : operations)
    {
        if(it == "flux" )
        {
            std::swap(it, operations.back());
            break;
        }
    }
}
template class Stats<double>;
template class Stats<float>;
