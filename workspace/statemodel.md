# Statemodel

- ChemNonEq1T

Sets the current state of the mixture.  Variable sets can be the
 following:

0: conserved variables (species densities, total energy density)

1: primitive set 1 (species densities, mixture temperature)

2: primitive set 2 (species mass fractions, {P, T} array)

- ChemNonTTv

Sets the current state of the mixture.  Variable sets can be the following:

​    0: {species densities}, {total energy density, vib. energy density}

​    1: {species densities} and {T, Tv}    

- Equil

Sets the mixture state by computing the equilibrium composition at the

given mixture conditions based on the given variable set. Variable sets
 can be the following:

0: conserved variables (mixture density, static energy density)

1: primitive set 1 (pressure, temperature)

2: primitive set 2 (elemental mole fractions,  {P, T} array)

- EquilTP

1:Sets the state of the equilibrium mixture given temperature and pressure.