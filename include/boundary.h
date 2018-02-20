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

#ifndef BOUNDARY
#define BOUNDARY

class Master;
template<typename> class Grid;
template<typename> class Fields;
class Input;
// struct Mask;

enum class Boundary_type   {Dirichlet_type, Neumann_type, Flux_type, Ustar_type};
enum class Boundary_w_type {Normal_type, Conservation_type};

/**
 * Base class for the boundary scheme.
 * This class handles the case when the boundary is turned off. Derived classes are
 * implemented that handle different boundary schemes.
 */
template<typename TF>
class Boundary
{
    public:
        Boundary(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constuctor of the boundary class.
        virtual ~Boundary();      ///< Destructor of the boundary class.

        static std::shared_ptr<Boundary> factory(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Factory function for boundary class generation.

        virtual void init(Input&);   ///< Initialize the fields.
        virtual void create(Input&); ///< Create the fields.

        // virtual void update_time_dependent(); ///< Update the time dependent parameters.

        virtual void set_values(); ///< Set all 2d fields to the prober BC value.

        virtual void exec(); ///< Update the boundary conditions.
        virtual void set_ghost_cells_w(Boundary_w_type); ///< Update the boundary conditions.

        // virtual void exec_stats(Mask*); ///< Execute statistics of surface
        // virtual void exec_cross();       ///< Execute cross sections of surface

        // virtual void get_mask(Field3d*, Field3d*, Mask*); ///< Calculate statistics mask
        // virtual void get_surface_mask(Field3d*);          ///< Calculate surface mask

        std::string get_switch();

        // GPU functions and variables
        // virtual void prepare_device();
        // virtual void forward_device();
        // virtual void backward_device();

    protected:
        Master& master; ///< Pointer to master class.
        Grid<TF>& grid;   ///< Pointer to grid class.
        Fields<TF>& fields; ///< Pointer to fields class.

        std::string swboundary;

        Boundary_type mbcbot;
        Boundary_type mbctop;

        TF ubot;
        TF utop;
        TF vbot;
        TF vtop;

        /**
         * Structure containing the boundary options and values per 3d field.
         */
        struct Field3dBc
        {
            TF bot; ///< Value of the bottom boundary.
            TF top; ///< Value of the top boundary.
            Boundary_type bcbot; ///< Switch for the bottom boundary.
            Boundary_type bctop; ///< Switch for the top boundary.
        };

        typedef std::map<std::string, Field3dBc> BcMap;
        BcMap sbc;

        // Variables to handle time dependency.
        std::string swtimedep;
        // std::vector<double> timedeptime;
        std::vector<std::string> timedeplist;
        // std::map<std::string, double*> timedepdata;

        void process_bcs(Input&); ///< Process the boundary condition settings from the ini file.

        // void process_time_dependent(Input *); ///< Process the time dependent settings from the ini file.

        // void set_bc(double*, double*, double*, Boundary_type, double, double, double); ///< Set the values for the boundary fields.

        // GPU functions and variables
        void set_bc_g(TF*, TF*, TF*, Boundary_type, TF, TF, TF); ///< Set the values for the boundary fields.

    private:
        virtual void update_bcs();       ///< Update the boundary values.
        virtual void update_slave_bcs(); ///< Update the slave boundary values.
};
#endif
