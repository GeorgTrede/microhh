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

#ifndef FIELDS
#define FIELDS
#include <map>
#include <vector>
#include "field3d.h"

// forward declarations to reduce compilation time
class Master;
class Input;
class Model;
class Grid;
class Stats;
struct Mask;

typedef std::map<std::string, Field3d *> fieldmap;

class Fields
{
  public:
    // functions
    Fields(Model *, Input *); ///< Constructor of the fields class.
    ~Fields();                ///< Destructor of the fields class

    void init();          ///< Initialization of the field arrays
    void create(Input *); ///< Initialization of the fields (random perturbations, vortices)  

    int exec();
    int getMask(Field3d *, Field3d *, Mask *);
    int execStats(Mask *);

    int initmomfld(Field3d*&, Field3d*&, std::string, std::string, std::string);
    int initpfld(std::string, std::string, std::string);
    int initdfld(std::string, std::string, std::string);
    int inittmpfld(std::string, std::string, std::string);
    
    void save(int);
    void load(int);

    double checkmom ();
    double checktke ();
    double checkmass();

    int setcalcprofs(bool);

    void execCross();

    Field3d *u; ///< Field3d instance of x velocity component
    Field3d *v; ///< Field3d instance of y velocity component
    Field3d *w; ///< Field3d instance of vertical velocity component

    Field3d *ut; ///< Field3d instance of x velocity component tendency 
    Field3d *vt; ///< Field3d instance of y velocity component tendency
    Field3d *wt; ///< Field3d instance of vertical velocity component tendency 

    fieldmap a;  ///< Map containing all field3d instances
    fieldmap ap; ///< Map containing all prognostic field3d instances
    fieldmap at; ///< Map containing all tendency field3d instances

    fieldmap mp; ///< Map containing all momentum field3d instances
    fieldmap mt; ///< Map containing all momentum tendency field3d instances

    fieldmap sd; ///< Map containing all diagnostic scalar field3d instances
    fieldmap sp; ///< Map containing all prognostic scalar field3d instances
    fieldmap st; ///< Map containing all prognostic scalar tendency field3d instances

    fieldmap atmp; ///< fieldmap containing all temporary field3d instances

    double *rhoref;  ///< Reference density at full levels 
    double *rhorefh; ///< Reference density at half levels

    // TODO remove these to and bring them to diffusion model
    double visc;

    /* 
     *Device (GPU) functions and variables
     */
    enum OffsetType {Offset, NoOffset};

    void prepareDevice();  ///< Allocation of all fields at device 
    void forwardDevice();  ///< Copy of all fields from host to device
    void backwardDevice(); ///< Copy of all fields required for statistics and output from device to host
    void clearDevice();    ///< Deallocation of all fields at device

    void forward3DFieldDevice (double *, double *, OffsetType); ///< Copy of a single 3d field from host to device
    void forward2DFieldDevice (double *, double *, OffsetType); ///< Copy of a single 2d field from host to device
    void forward1DFieldDevice (double *, double *, int);        ///< Copy of a single array from host to device
    void backward3DFieldDevice(double *, double *, OffsetType); ///< Copy of a single 3d field from device to host
    void backward2DFieldDevice(double *, double *, OffsetType); ///< Copy of a single 2d field from device to host
    void backward1DFieldDevice(double *, double *, int);        ///< Copy of a single array from device to host

    double *rhoref_g;  ///< Reference density at full levels at device
    double *rhorefh_g; ///< Reference density at half levels at device
    
  private:
    // variables
    Model  *model;
    Grid   *grid;
    Master *master;
    Stats  *stats;

    bool calcprofs;

    // cross sections
    std::vector<std::string> crosslist;      // List with all crosses from ini file
    // Cross sections split per type
    std::vector<std::string> crosssimple;
    std::vector<std::string> crosslngrad;   
    std::vector<std::string> crossbot;
    std::vector<std::string> crosstop;
    std::vector<std::string> crossfluxbot;
    std::vector<std::string> crossfluxtop;
    int checkaddcross(std::string, std::string, std::vector<std::string> *, std::vector<std::string> *);

    // masks
    int calcmaskwplus(double *, double *, double *, int *, int *, int *, double *);
    int calcmaskwmin (double *, double *, double *, int *, int *, int *, double *);

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
    int addmeanprofile(Input *, std::string, double *, double);
    int randomnize(Input *, std::string, double *);
    int addvortexpair(Input* inputin);
    double calcmass(double *, double *);

    // statistics
    double *umodel;
    double *vmodel;

    /* 
     *Device (GPU) functions and variables
     */
    void forwardField3dDevice(Field3d *);  ///< Copy of a complete Field3d instance from host to device
    void backwardField3dDevice(Field3d *); ///< Copy of a complete Field3d instance from device to host
};
#endif
