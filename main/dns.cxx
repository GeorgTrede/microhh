#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "timeloop.h"
#include "advec.h"
#include "diff.h"
#include "force.h"
#include "pres.h"

int main()
{
  // read the input data
  cinput input;
  if(input.readinifile())
    return 1;

  cmpi mpi;
  if(input.readinifile())
    return 1;
  
  // create the objects, read the inputdata
  cgrid grid;
  if(grid.readinifile(&input))
    return 1;

  cfields fields(&grid);
  if(fields.readinifile(&input))
    return 1;

  // initialize the objects, allocate the required memory
  mpi.init();
  grid.initgrid();
  fields.initfields();

  // create the model and the operators
  ctimeloop timeloop(&grid, &fields);
  cadvec    advec   (&grid, &fields);
  cdiff     diff    (&grid, &fields);
  cpres     pres    (&grid, &fields);
  cforce    force   (&grid, &fields);

  // read the inputdata
  if(timeloop.readinifile(&input))
    return 1;

  // fill the fields with data
  if(grid.load())
    return 1;
  if(pres.load())
    return 1;
  if(timeloop.load(timeloop.iteration))
    return 1;
  if(fields.load(timeloop.iteration))
    return 1;

  // initialize the diffusion to get the time step requirement
  diff.init();

  // initialize the pressure solver
  pres.init();

  // initialize the check variables
  int    iter;
  double time, dt, cputime;
  double mom, tke, mass;
  double div;
  double cfl;
  
  std::printf("%8s  %12s  %10s  %10s  %8s  %13s  %13s  %13s  %13s\n", 
    "ITER", "TIME", "CPUDT", "DT", "CFL", "DIV", "MOM", "TKE", "MASS");

  // set the boundary conditions
  fields.boundary();

  // start the time loop
  while(timeloop.loop)
  {
    // determine the time step
    if(!timeloop.insubstep())
    {
      cfl = advec.getcfl(timeloop.dt);
      timeloop.settimestep(cfl);
    }

    // advection
    advec.exec();
    // diffusion
    diff.exec();
    // large scale forcings
    force.exec(timeloop.getsubdt());
    // pressure
    pres.exec(timeloop.getsubdt());
    // perform the timestepping substep
    timeloop.exec();

    // boundary conditions
    fields.boundary();

    if(!timeloop.insubstep())
      timeloop.timestep();

    if(timeloop.docheck() && !timeloop.insubstep())
    {
      iter    = timeloop.iteration;
      time    = timeloop.time;
      dt      = timeloop.dt;
      cputime = timeloop.check();
      div     = pres.check();
      mom     = fields.check(0);
      tke     = fields.check(1);
      mass    = fields.check(2);

      std::printf("%8d  %12.3f  %10.4f  %10.4f  %8.4f  %13.5E  %13.5E  %13.5E  %13.5E\n", 
          iter, time, cputime, dt, cfl, div, mom, tke, mass);
    }

    if(timeloop.dosave() && !timeloop.insubstep())
    {
      timeloop.save(timeloop.iteration);
      fields.save(timeloop.iteration);
    }
  }
  
  return 0;
}

