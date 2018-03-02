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

#ifndef FIELD3D_IO
#define FIELD3D_IO

class Master;
template<typename> class Grid;

template<typename TF>
class Field3d_io
{
    public:
        Field3d_io(Master&, Grid<TF>&);
        ~Field3d_io();

        int save_field3d(TF*, TF*, TF*, char*, const TF); // Saves a full 3d field.
        int load_field3d(TF*, TF*, TF*, char*, const TF); // Loads a full 3d field.

        int save_xz_slice(TF*, TF*, char*, int);           // Saves a xz-slice from a 3d field.
        int save_yz_slice(TF*, TF*, char*, int);           // Saves a yz-slice from a 3d field.
        int save_xy_slice(TF*, TF*, char*, int kslice=-1); // Saves a xy-slice from a 3d field.
        int load_xy_slice(TF*, TF*, char*, int kslice=-1); // Loads a xy-slice.

    private:
        Master& master;
        Grid<TF>& grid;
};
#endif

