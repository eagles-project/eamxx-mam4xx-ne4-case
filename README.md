# Eamxx-Mam4xx Nucleation Demo

This repo is a set of utility scripts and working notes for assembling a simulation that demonstrates the
nucleation of H2SO4 gas into SO4 aerosol in the aitken mode within EAMxx.

## E3SM Tests

These might be good to run to make sure you've installed E3SM correctly, and should work under EAMxx/SCREAM builds configured with CIME.

* `SMS_D.ne4pg2_oQU480.F2010`

Compsets

* `SMS_D_Ld1_P96x1.ne4pg2_ne4pg2.F2010-SCREAMv1`

Here's what all the symbols mean:

* `D` indicates a Debug build, and can be removed for optimized builds
* `L` indicates the length of the simulation (e.g. Ld1 means "run for 1 day"). The default simulation
  length is 5 days.
* `P` indicates the number of processors and threads to use, delimited by x. So P96x1 means 96 MPI
  processes and a single thread per process. Note that the ne4 grid has 96 = 4x4x6 spectral elements.
* `ne4pg2_ne4pg2` here indicates that the ne4 grid is to be used for dynamics and the pg2 grid for
  physics, in contrast to `ne4_ne4`, which uses the ne4 grid for both.
* `F2010-SCREAMv1` is the compset (or set of components) corresponding to the SCREAM atmospheric model
  with the other components (ocean, land, land/sea ice, etc) supplied by data files.

## Files of interest:

* `config_machines.xml`: a little auxiliary machine file for Jeff's Linux workstation
* `add_aerosol_ics`: a Python script that adds aerosol tracers to the default F-case
  initial conditions file (requires the Python `netCDF4` module)
* `create_nuc_case`: a Bash script that generates and builds the nucleation demo case

## How to Build and Run the Demo

These commands download any needed simulation data and then run the case:

```
./create_nuc_case`
./case.submit --no-batch
```

## Resources

* [CIME Documentation](http://esmci.github.io/cime/versions/master/html/index.html)
* [How to Run SCREAMV1](https://acme-climate.atlassian.net/wiki/spaces/NGDNA/pages/3386015745/How+To+Run+SCREAMv1)
