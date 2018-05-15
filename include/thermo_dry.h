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

#ifndef THERMO_DRY
#define THERMO_DRY

#include "thermo.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
class Data_block;


/**
 * Class for the dry thermodynamics.
 * This class is responsible for the computation of the right hand side term related to
 * the acceleration by buoyancy. In the dry thermodynamics temperature and buoyancy are
 * equivalent and no complex buoyancy function is required.
 */

enum class Basestate_type {anelastic, boussinesq};

template<typename TF>
class Thermo_dry : public Thermo<TF>
{
    public:
        Thermo_dry(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the dry thermodynamics class.
        virtual ~Thermo_dry(); ///< Destructor of the dry thermodynamics class.

        void init();
        void create(Input&, Data_block&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&);
        void exec(); ///< Add the tendencies belonging to the buoyancy.
        unsigned long get_time_limit(unsigned long, double); ///< Compute the time limit (n/a for thermo_dry)

        void exec_stats(Stats<TF>&, std::string, Field3d<TF>&, Field3d<TF>&, const Diff<TF>&);
        void exec_cross(Cross<TF>&, unsigned long);
        void exec_dump(Dump<TF>&, unsigned long);
        void exec_column(Column<TF>&);

        bool check_field_exists(std::string name);
        struct background_state
        {
            Basestate_type swbasestate;

            TF pbot;   ///< Surface pressure.
            TF thref0; ///< Reference potential temperature in case of Boussinesq

            std::vector<TF> thref;
            std::vector<TF> threfh;
            std::vector<TF> pref;
            std::vector<TF> prefh;
            std::vector<TF> exnref;
            std::vector<TF> exnrefh;

            // GPU functions and variables
            TF*  thref_g;
            TF*  threfh_g;
            TF*  pref_g;
            TF*  prefh_g;
            TF*  exnref_g;
            TF*  exnrefh_g;
        };
        background_state bs;
        background_state bs_stats;
        void get_thermo_field(Field3d<TF>&, std::string, bool, background_state);
        void get_buoyancy_surf(Field3d<TF>&, background_state);
        void get_buoyancy_fluxbot(Field3d<TF>&, background_state);
        void get_T_bot(Field3d<TF>&, background_state);
        const std::vector<TF>& get_p_vector() const;
        const std::vector<TF>& get_ph_vector() const;

        void get_prog_vars(std::vector<std::string>&); ///< Retrieve a list of prognostic variables.
        TF get_buoyancy_diffusivity();

        #ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        void forward_device();
        void backward_device();
        void get_thermo_field_g(Field3d<TF>&, std::string, bool);
        void get_buoyancy_surf_g(Field3d<TF>&);
        void get_buoyancy_fluxbot_g(Field3d<TF>&);
        #endif

        // Empty functions that are allowed to pass.
        void get_mask(Field3d<TF>&, Field3d<TF>&, Stats<TF>&, std::string) {};
        void update_time_dependent() {};

    private:
        using Thermo<TF>::swthermo;
        using Thermo<TF>::master;
        using Thermo<TF>::grid;
        using Thermo<TF>::fields;

        Boundary_cyclic<TF> boundary_cyclic;

        // cross sections
        std::vector<std::string> crosslist;        ///< List with all crosses from ini file
        std::vector<std::string> allowedcrossvars; ///< List with allowed cross variables
        std::vector<std::string> dumplist;         ///< List with all 3d dumps from the ini file.

        void create_stats(Stats<TF>&);   ///< Initialization of the statistics.
        void create_column(Column<TF>&); ///< Initialization of the single column output.
        void create_dump(Dump<TF>&);     ///< Initialization of the single column output.
        void create_cross(Cross<TF>&);   ///< Initialization of the single column output.


};
#endif
