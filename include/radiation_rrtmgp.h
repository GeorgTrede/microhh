/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#ifndef RADIATION_RRTMGP_H
#define RADIATION_RRTMGP_H

#include "radiation.h"
#include "field3d_operators.h"

#include "Gas_optics.h"
#include "Gas_concs.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
template<typename> class Thermo;
template<typename> class Timeloop;

// RRTMGP classes.
template<typename> class Gas_concs;
template<typename> class Gas_optics;

template<typename TF>
class Radiation_rrtmgp : public Radiation<TF>
{
	public:
		Radiation_rrtmgp(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Radiation_rrtmgp() {}

		bool check_field_exists(std::string name)
        { throw std::runtime_error("Not implemented"); }

        void init();
        void create(Thermo<TF>&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&);
        void exec(Thermo<TF>&, double, Timeloop<TF>&, Stats<TF>&);

		// Empty functions that should throw
		// Stats/dump not implemented in rrmtg now
        void get_radiation_field(Field3d<TF>&, std::string, Thermo<TF>&, Timeloop<TF>&)
        { throw std::runtime_error("Not implemented"); }

        void exec_stats(Stats<TF>&, Thermo<TF>&, Timeloop<TF>&){};
        void exec_cross(Cross<TF>&, unsigned long, Thermo<TF>&, Timeloop<TF>&){};
        void exec_dump(Dump<TF>&, unsigned long, Thermo<TF>&, Timeloop<TF>&){};
        void exec_column(Column<TF>&, Thermo<TF>&, Timeloop<TF>&){};

	private:
		using Radiation<TF>::swradiation;
		using Radiation<TF>::master;
		using Radiation<TF>::grid;
		using Radiation<TF>::fields;
		using Radiation<TF>::field3d_operators;

        const std::string tend_name = "rad";
        const std::string tend_longname = "Radiation";

        Gas_concs<TF> gas_concs;
        std::unique_ptr<Gas_optics<TF>> kdist_lw;
        std::unique_ptr<Gas_optics<TF>> kdist_sw;
};
#endif
