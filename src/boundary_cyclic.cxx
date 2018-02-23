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

#include "master.h"
#include "grid.h"
#include "boundary_cyclic.h"

template<typename TF>
Boundary_cyclic<TF>::Boundary_cyclic(Master& masterin, Grid<TF>& gridin)
    master(masterin),
    grid(gridin)
{
}

template<typename TF>
Boundary_cyclic<TF>::~Boundary_cyclic()
{
}

template<typename TF>
void Boundary_cyclic<TF>::init()
{
}

template<typename TF>
void Boundary_cyclic<TF>::create()
{
}

