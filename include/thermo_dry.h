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

#ifndef THERMO_DRY
#define THERMO_DRY

#include "thermo.h"

// forward declarations to speed up build time
class Master;
class Grid;
class Fields;
class Stats;

/**
 * Class for the dry thermodynamics.
 * This class is responsible for the computation of the right hand side term related to
 * the acceleration by buoyancy. In the dry thermodynamics temperature and buoyancy are
 * equivalent and no complex buoyancy function is required.
 */
class ThermoDry : public Thermo
{
  public:
    ThermoDry(Model *, Input *); ///< Constructor of the dry thermodynamics class.
    virtual ~ThermoDry();                  ///< Destructor of the dry thermodynamics class.

    virtual void init();
    virtual void create(Input *);
    virtual void exec();                ///< Add the tendencies belonging to the buoyancy.
    virtual void execStats(Mask *);
    virtual void execCross();

    virtual bool checkThermoField(std::string name);
    virtual void getThermoField(Field3d *, Field3d *, std::string name);
    virtual void getBuoyancySurf(Field3d *);              ///< Compute the near-surface and bottom buoyancy for usage in another routine.
    virtual void getBuoyancyFluxbot(Field3d *);           ///< Compute the bottom buoyancy flux for usage in another routine.
    virtual void getProgVars(std::vector<std::string> *); ///< Retrieve a list of prognostic variables.

    #ifdef USECUDA
    // GPU functions and variables
    virtual void prepareDevice();
    virtual void clearDevice();
    #endif

  private:
    int calcbuoyancy(double *, double *, double *);     ///< Calculation of the buoyancy.
    int calcN2(double *, double *, double *, double *); ///< Calculation of the Brunt-Vaissala frequency.
    
    // cross sections
    std::vector<std::string> crosslist;        ///< List with all crosses from ini file
    std::vector<std::string> allowedcrossvars; ///< List with allowed cross variables

    int calcbuoyancybot(double *, double *,
                        double *, double *,
                        double *, double *);                ///< Calculation of the near-surface and surface buoyancy.
    int calcbuoyancyfluxbot(double *, double *, double *);  ///< Calculation of the buoyancy flux at the bottom.
    int calcbuoyancytend_2nd(double *, double *, double *); ///< Calculation of the buoyancy tendency with 2nd order accuracy.
    int calcbuoyancytend_4th(double *, double *, double *); ///< Calculation of the buoyancy tendency with 4th order accuracy.

    Stats *stats;

    double pbot;   ///< Surface pressure.
    double thref0; ///< Reference potential temperature in case of Boussinesq

    double *thref;
    double *threfh;
    double *pref;
    double *prefh;
    double *exner;
    double *exnerh;

    // GPU functions and variables
    double *thref_g;
    double *threfh_g;
    double *pref_g;
    double *prefh_g;
    double *exner_g;
    double *exnerh_g;

};
#endif
