/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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

#ifndef DIFF_G2
#define DIFF_G2

#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "mpiinterface.h"

class cdiff_g2 : public cdiff
{
  public:
    cdiff_g2(cgrid *, cfields *, cmpi *);
    ~cdiff_g2();

    int setvalues();
    int exec();

    unsigned long gettimelim(unsigned long, double);
    double getdn(double);

  private:
    double dnmul;

    int diffc(double *, double *, double *, double *, double);
    int diffw(double *, double *, double *, double *, double);
};
#endif
