# Input and Restart Databases in IBAMR

## Overview

Most SAMRAI objects (and, hence, most IBAMR objects) are constructed with an
input database. For example, `IBExplicitHierarchyIntegrator` has a single
constructor with the following arguments:

```cpp
IBExplicitHierarchyIntegrator(std::string object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                              SAMRAI::tbox::Pointer<IBStrategy> ib_method_ops,
                              SAMRAI::tbox::Pointer<INSHierarchyIntegrator> ins_hier_integrator,
                              bool register_for_restart = true);
```

Of these, the second argument configures the object every time it is constructed
and the first, second, and last arguments are used when the application is
restarted.

## Initial configuration is done with input databases

Here the second argument is the input database, which contains all the
parameters that configure that class. For example, the constructor of this class
decides whether or not to use a structure predictor if `use_structure_predictor`
is provided as a boolean to the input database.

These databases are specified in the input file. For example, one would create
an `IBExplicitHierarchyIntegrator` in the following way:
```cpp
Pointer<IBHierarchyIntegrator> time_integrator =
    new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                      app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                      ib_method_ops,
                                      navier_stokes_integrator);
```
Here `app_initializer->getComponentDatabase()` gets a database from the input
file named `IBHierarchyIntegrator`. This is spelt out in the input file provided
to the application. Hence, finally, the input file might contain the following
input database:
```
IBHierarchyIntegrator {
   start_time          = START_TIME
   end_time            = END_TIME
   grow_dt             = GROW_DT
   num_cycles          = NUM_CYCLES
   regrid_cfl_interval = REGRID_CFL_INTERVAL
   dt_max              = DT
   error_on_dt_change  = ERROR_ON_DT_CHANGE
   enable_logging      = ENABLE_LOGGING

   use_structure_predictor = FALSE
}
```

which will completely configure the `IBExplicitHierarchyIntegrator` object. Some
options here are used by the base class. IBAMR uses the convention that
constants common to multiple classes are written in all caps and are defined
at the top level of the input database.

## Restart configuration is done with a restart database

IBAMR applications tend to be very computationally expensive and can require
days or weeks to finish. Hence all IBAMR objects can store their current state
to a file so that the computation can be stopped and started several times.

Many objects store state that changes over the duration of a computation (such
as the current time, the center of mass of a structure, or a velocity field):
the state of the object at any point in time can be saved to an output file with
the `putToDatabase()` and `putToDatabaseSpecialized()` methods. SAMRAI
associates objects with the restart database with their object names (which is
typically the first argument to the constructor): i.e., data is saved or loaded
from a part of the restart file which is identified by the object name. Objects
don't have to be registered with the restart database: in the case of
`IBExplicitHierarchyIntegrator` the last argument to the constructor is a
boolean indicating whether or not the object should be placed in the restart
database.

IBAMR uses SAMRAI's `RestartManager` to store information for restarting a
computation. Hence, most IBAMR object constructors contain the following code
after the input database is read in:

```cpp
// Initialize object with data read from the input and restart databases.
bool from_restart = RestartManager::getManager()->isFromRestart();
if (from_restart) getFromRestart();
```

Hence, object initialization is done in several stages:
1. Set up everything provided by the constructor that is not in the input database.
2. Set up everything provided by the input database.
3. If the current application is a restart (i.e., there is something in the
   restart database), then read the current state of the object from the restart
   file.

Similarly, an object should save the information that describes its current
state (not its initial configuration) into the restart database, since the same
input file will always be used for the initial configuration of the object.
