/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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

#ifndef PRES
#define PRES

// forward declarations to speed up build time
class cmodel;
class cgrid;
class cfields;
class cmaster;

class cpres
{
  public:
    cpres(cmodel *, cinput *);
    virtual ~cpres();
    static cpres* factory(cmaster *, cinput *, cmodel *, const std::string); ///< Factory function for pres class generation.

    virtual void init();
    virtual void setvalues();

    virtual void exec(double);
    virtual double check();

    virtual int prepareGPU();

  protected:
    cmaster *master;
    cmodel  *model;
    cgrid   *grid;
    cfields *fields;
};
#endif
