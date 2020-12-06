# MD

Performs a classical propagation of a system on a potential energy surface (PES)
given its initial conditions. Time propagation is perform with leapfrog algorithm.

## Input file
A file named `input_file.py` is requiered with all the information to perform the 
molecular dynamcics. A template can be found in `input_file_example`. It has the 
following general form:
```python
input_file={
        'nat' : 4,                                   # Number of particles in the system.
        'mintraj' : 1,                               # Initial trajectory.
        'maxtraj' : 10,                              # Final trajectory.
        'initial_cond' : 'initial-conditions-A.dat', # File with initial conditions.
        'mass' : [12e0, 15e0, 1e0, 1e0],             # Masses of the particles (amu).
        'tfin' : 100,                                # Final propagation time.
        'tstep' : 0.1,                               # Time step for leapfrog propagation.
        'it_print' : 10,                             # Printing interval.
        'atoms' : ['C', 'O', 'H', 'H'],              # List with the atomic symbol of the particles.
        'random_init_cond' : True,                   # True = Random indexes for selecting initial conditions.
                                                     # False = Sequencial indexes.
        }

def check_end(t, XP):
  '''
  Function to decide whether a trayectory ended or not.
  It is usefull if the trajectory will be stopped if some distance breaks, for example.
  '''
  
  end = False
  return end
```
## Initial conditions
For a system of N particles, each initial condition is expressed as a 2\*3\*N vector, where the 
first 3\*N elements are the positions of each particles and the following 3\*N elements their 
momentum.

Initial conditions must be provided in atomic units.

## Potential energy surface
Your PES is called through a function named `potxyz` where a vector of positions (in atomic units) 
is input and the energy and derivatives with respect to each coordinate are returned (all of them in atomic units).

If your PES is a Fortran code, compile it through f2py with a module name pes (`-m pes`). 

# Execution
To run the molecular dynamics trajectories execute:
```
python MD.py
```

Several files will be created:
* end-conditions: File with position and momenta of the system at the end of the propagation.
* traj_\<id\>.xyz: File where `<id>` trajectory is printed if `it_print != 0`.
* Standard output: Prints the trajectory currently propagating as well as initial and final energy to check its 
conservation.
