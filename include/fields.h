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

#ifndef FIELDS
#define FIELDS
#include <map>
#include "field3d.h"

// forward declarations to reduce compilation time
class cmaster;
class cmodel;
class cgrid;
class cstats;

typedef std::map<std::string, cfield3d *> fieldmap;

class cfields
{
  public:
    // functions
    cfields(cmodel *);
    ~cfields();
    int readinifile(cinput *);
    int init();
    int create(cinput *);
    int exec();
    int statsexec();

    int initmomfld(cfield3d*&, cfield3d*&, std::string, std::string, std::string);
    int initpfld(std::string, std::string, std::string);
    int initdfld(std::string, std::string, std::string);
    
    int save(int);
    int load(int);

    double checkmom ();
    double checktke ();
    double checkmass();

    int setcalcprofs(bool);
    int execcross();

    // 3d fields for momentum
    cfield3d *u;
    cfield3d *v;
    cfield3d *w;

    cfield3d *ut;
    cfield3d *vt;
    cfield3d *wt;

    // maps of 3d fields
    fieldmap a;
    fieldmap ap;
    fieldmap at;

    fieldmap mp;
    fieldmap mt;

    fieldmap s;
    fieldmap sd;
    fieldmap sp;
    fieldmap st;

    // TODO remove these to and bring them to diffusion model
    double visc;

    double *rhoref;

  private:
    // variables
    cmodel  *model;
    cgrid   *grid;
    cmaster *master;
    cstats  *stats;

    bool allocated;
    bool calcprofs;

    // cross sections
    std::vector<std::string> crosslist;      // List with all crosses from ini file
    std::vector<std::string> csimple;        // Cross sections split per type
    std::vector<std::string> clngrad;        //  ""      ""      ""   
    std::vector<std::string> cbot;
    std::vector<std::string> ctop;
    std::vector<std::string> cfluxbot;
    std::vector<std::string> cfluxtop;
    int checkaddcross(std::string, std::string, std::vector<std::string> *, std::vector<std::string> *);

    // perturbations
    double rndamp;
    double rndz;
    double rndexp;
    double vortexamp;
    int vortexnpair;
    std::string vortexaxis;
    
    // functions
    double calcmom_2nd(double *, double *, double *, double *);
    double calctke_2nd(double *, double *, double *, double *);
    int addmeanprofile(cinput *, std::string, double *, double);
    int randomnize(cinput *, std::string, double *);
    int addvortexpair(cinput* inputin);
    double calcmass(double *, double *);
    inline double interp2(const double, const double);

    // statistics
    double *umodel;
    double *vmodel;
};
#endif
