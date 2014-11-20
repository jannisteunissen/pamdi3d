Particle simulations of discharges in 3D, with adaptive particle control
(super-particles) and adaptive grid refinement.

## Getting the code on your machine

First, clone the repository (in this case into a folder pamdi3d)
```
$ git clone https://github.com/jannisteunissen/pamdi3d pamdi3d
```

The code uses git submodules, so for the first compilation you should initialize
and update the submodules.

```
$ cd pamdi3d
$ git submodule init
$ git submodule update
```

After that, you can compile with just

```
$ make
```

## Updating the code

Depending on your git preferences you can use

```
$ git fetch
```

or

```
$ git pull; git merge origin/master
```

The submodules might be updated, so after that, run

```
$ git submodule update
```

To do a fresh compilation:

```
$ make clean
$ make
```

## Running the code

Running (sequential):

```
$ ./pampi3d my_config_file.txt
```

Running (parallel, N = number of tasks):

```
$ mpirun -n N ./pampi3d my_config_file.txt
```

You can also specify multiple configuration files, like:

```
$ ./pampi3d cfg_base.txt cfg_1.txt
```

In each configuration file, you can specify a different "sim_name" variable.
These names will be appended to each other.

## Visualizing the output

I'd highly recommend
[Visit](https://wci.llnl.gov/simulation/computer-codes/visit/downloads) to
visualize the output. Output files have the
[Silo](https://wci.llnl.gov/simulation/computer-codes/silo) format.
