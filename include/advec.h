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

#ifndef ADVEC
#define ADVEC

#include <string>
#include <memory>
#include "field3d_operators.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;

/**
 * Base class for the advection scheme. This class is abstract and only
 * derived classes can be instantiated. Derived classes are
 * implemented that handle different advection schemes.
 */
template<typename TF>
class Advec
{
    public:
        Advec(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the advection class.
        virtual ~Advec(); ///< Destructor of the advection class.

        static std::shared_ptr<Advec> factory(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Factory function for advection class generation.

        std::string get_switch();

        // Pure virtual functions that have to be implemented in derived class.
        virtual void exec() = 0; ///< Execute the advection scheme.
        virtual unsigned long get_time_limit(unsigned long, double) = 0; ///< Get the maximum time step imposed by advection scheme
        virtual double get_cfl(double) = 0; ///< Retrieve the CFL number.

    protected:
        Master& master; ///< Pointer to master class.
        Grid<TF>& grid; ///< Pointer to grid class.
        Fields<TF>& fields; ///< Pointer to fields class.
        Field3d_operators<TF> field3d_operators; ///< Instance of the field3d_operators

        double cflmax; ///< Maximum allowed value for the CFL criterion.
        const double cflmin; ///< Minimum value for CFL used to avoid overflows.

        std::string swadvec;
};
#endif
