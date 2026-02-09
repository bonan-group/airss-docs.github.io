---
title: "Building Structures"
layout: single
classes: wide
lang: en
lang-ref: building-structures
sidebar:
  nav: "docs"
---

# How to Build Structures

This guide covers the practical workflow for building random structures with AIRSS, from creating seed files to launching searches.

## Typical Workflow - Single Composition Search

### Step 1: Generate a Base Cell with `gencell`

The `gencell` tool generates seed files for a given composition. The syntax is:

```console
$ gencell <volume> <units> <species1> <count1> [<species2> <count2> ...]
```

Where:
- `<volume>` - Target volume in Å³ per formula unit
- `<units>` - Number of formula units (typically 1)
- `<species> <count>` - Element symbol and count per formula unit

**Example: SrTiO₃**

```console
$ gencell 60 1 Sr 1 Ti 1 O 3
```

This creates three files:
- `SrTiO3.cell` - The seed file for buildcell
- `SrTiO3.param` - CASTEP parameters for DFT relaxation
- `SrTiO3.par` - Parameters for RASH (relax-and-shake)

The generated cell file contains a cubic lattice with edge length = volume^(1/3) ≈ 3.91 Å, with all atoms placed at the origin (0, 0, 0). The `buildcell` program will randomize these positions. 

The exact coordinates and *shape* of the lattice does not matter here: they will be randomized anyway.

The generated `SrTiO3.cell` file looks like this:

```
%BLOCK LATTICE_CART
3.914867641167553 0 0
0 3.914867641167553 0
0 0 3.914867641167553
%ENDBLOCK LATTICE_CART

#VARVOL=60

%BLOCK POSITIONS_FRAC
Sr 0.0 0.0 0.0 # Sr1 % NUM=1
Ti 0.0 0.0 0.0 # Ti1 % NUM=1
O 0.0 0.0 0.0 # O1 % NUM=1
O 0.0 0.0 0.0 # O2 % NUM=1
O 0.0 0.0 0.0 # O3 % NUM=1
%ENDBLOCK POSITIONS_FRAC

##SPECIES=Sr,Ti,O
##NATOM=3-9
##FOCUS=3

#SYMMOPS=2-4
#NFORM=1
#MINSEP=1-3 AUTO
#COMPACT
...
```

You can also create seed files manually. A minimal seed file looks like this:

```
%BLOCK LATTICE_CART
2 0 0
0 2 0
0 0 2
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Al 0.0 0.0 0.0 # Al1 % NUM=8
%ENDBLOCK POSITIONS_FRAC

#MINSEP=1.5
```

Key elements:
- **LATTICE_CART** or **LATTICE_ABC**: Defines the unit cell volume - shape does not matter here
- **POSITIONS_FRAC** or **POSITIONS_ABS**: Atom positions (all at origin is fine - buildcell randomizes them)
- **NUM=n**: Creates n copies of this atom in the structure
- **#MINSEP**: Minimum allowed distance between atoms (in Å)

For multi-species systems, you should specify pair-specific separations:

```
%BLOCK LATTICE_CART
4 0 0
0 4 0
0 0 4
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Sr 0.0 0.0 0.0 # Sr1 % NUM=1
Ti 0.0 0.0 0.0 # Ti1 % NUM=1
O  0.0 0.0 0.0 # O1  % NUM=3
%ENDBLOCK POSITIONS_FRAC

#VARVOL=60
#MINSEP=1.5 Sr-Sr=3.5 Sr-Ti=3.0 Sr-O=2.4 Ti-O=1.8 O-O=2.5
#NFORM=1
#SYMMOPS=2-4
```

### Step 2: Generate Random Structures

Use `buildcell` to generate random structures from the seed:

```console
$ buildcell < SrTiO3.cell > random.cell
```

To generate multiple structures for inspection:

```bash
for i in $(seq 1 10); do
  buildcell < SrTiO3.cell > SrTiO3-$i.cell
done
```

### Step 3: Visualize and Validate

Before launching a full search, generate a few structures and check them visually using tools like VESTA, ASE, or OVITO.

**Checklist:**
1. **Atomic distances** - Are atoms at physically reasonable separations?
2. **Spatial distribution** - Are there holes or clustering of species?
3. **Build success** - Does buildcell generate structures without hanging?

### Common Pitfalls and Solutions

| Problem | Symptoms | Solution |
|---------|----------|----------|
| Atoms too close | Short bonds, overlapping atoms | Add pair-specific MINSEP values for the problematic species pairs |
| Structure won't build | buildcell hangs or fails repeatedly | VARVOL is too small for the MINSEP constraints. Increase VARVOL until structures build successfully |
| No symmetry | All structures are P1 | Check that `#SYMMOPS` and `#NFORM` are set |
| Too much symmetry | Missing low-symmetry polymorphs | Use `#SYMMOPS=1-4` to include low-symmetry structures |

### Step 4: Launch the Search

Once you're satisfied with the seed file, launch the search:

```console
$ airss.pl -castep -max 100 -seed SrTiO3
```

Or with a pair potential for faster testing:

```console
$ airss.pl -pp3 -max 100 -seed SrTiO3
```

> **Note:** Start with small formula units (`#NFORM=1` or `#NFORM=1-2`) to validate the search setup before scaling up.
> **Note:** In HPC environment, you should launch as series of `airss.pl` tasks to in the search in parallel. The easiest way to achieve this is to use **job array**. You should, however, wait for the appearance of a few `.res` file to reassure that the search is process before submit a large amount of jobs as job array to avoid wasting computing resources on errorous parameters.



### Step 5: Monitor the search

Analyze search results with:

```console
$ ca -r                  # Rank by energy
$ ca -u 0.01 -r          # Unify similar structures
$ ca -s -cl -r           # Summary with clustering
```

You should stop the search if the lowest energies structure have been found multiple times. The number of identical structures found is the last column when the `-u` flags is used for unifying similar structures. Otherwise this columns is always `1`.

---

## Variable Composition Searches

AIRSS can explore multiple compositions simultaneously using the `#SPECIES` tag and related directives.

### The `{}` Expansion Syntax

The `{a,b,c,...}` syntax allows random selection from a set of values. Each value has equal probability, but you can weight the selection by repeating values:

```
#NFORM={2,2,2,4,4,6}
```

This selects 2 with probability 3/6, 4 with probability 2/6, and 6 with probability 1/6.

Compare with range syntax:

```
#NFORM=2-6
```

This gives uniform probability to integers 2, 3, 4, 5, and 6.

The `{}` expansion works for any directive and is processed before other parsing.

### The SPECIES Tag

The `#SPECIES` tag enables variable composition searching:

```
#VARVOL=15
#SPECIES=A,B
#NATOM=2-8
#MINSEP=1.5
```

This explores all A-B compositions with 2-8 total atoms, including A₂, AB, A₂B, AB₂, A₃B, etc.

### Fixed Stoichiometry with NUM

To fix the ratio between species, use `%NUM=` modifiers:

```
#SPECIES=Si%NUM=1,O%NUM=2
#NFORM=4
```

This generates SiO₂ structures with 4 formula units (4 Si + 8 O = 12 atoms).

### The FOCUS Tag

When running variable composition searches, use `#FOCUS=n` to restrict generation to a specific composition:

```
#SPECIES=A,B,C
#NATOM=4-8
#FOCUS=1
```

This only generates structures with the first composition variant. Useful for parallelizing searches across compositions.

### Example: Ternary Variable Composition

```
#VARVOL=15
#SPECIES=A,B,C
#NATOM=3-8
#MINSEP=1.5
```

This explores all ternary compositions from 3-8 atoms.

> **Warning:** Variable composition searches are computationally expensive. The number of possible compositions grows rapidly with the number of species and atom count. Start with small ranges and expand gradually.

---

## Tips and Tricks

### Weighted Distributions

Use `{}` expansion for non-uniform sampling:

```
#SYMMOPS={1,1,2,2,2,3,4}    # Favor low symmetry
#NFORM={1,1,2,2,4}          # Favor small cells
```

### Pair-Specific MINSEP for Diversity

Always specify MINSEP for different species pairs:

```
#MINSEP=1.0 Si-Si=3.00 Si-O=1.60 O-O=2.58
```

Benefits:
- **Physical relevance**: Different atom pairs have different equilibrium distances
- **Better sampling**: Prevents unphysical configurations early
- **Improved diversity**: Allows closer approach for some pairs while preventing others

### Visualization Tools

| Tool | Usage | Notes |
|------|-------|-------|
| VESTA | GUI application | Excellent for crystal structures |
| ASE | `ase gui structure.cell` | Python-based, scriptable |
| OVITO | GUI application | Good for large systems |

To convert results to extended XYZ format:

```console
$ cabal res xyz
```

### Using Range Syntax

Many directives support range syntax `min-max`:

```
#NFORM=1-4        # 1, 2, 3, or 4 formula units
#SYMMOPS=2-6      # 2 to 6 symmetry operations
#VARVOL=50-80     # Volume between 50-80 Å³
#MINSEP=1.5-2.0   # Global minimum separation varies
```

---

## Guessing Starting Points for VARVOL and MINSEP

### Understanding `VARVOL`

The `#VARVOL` tag specifies the target volume in Å³ (cubic Angstroms) for the **unexpanded** cell - that is, the cell before any expansion from `#NFORM` or `NUM` tags.

```
#VARVOL=60
```

When you use `#NFORM` or `NUM` to create multiple formula units, the volume is automatically scaled. For example:

- `#VARVOL=60` with `#NFORM=1` → target volume of 60 Å³
- `#VARVOL=60` with `#NFORM=2` → target volume of 120 Å³ (automatically scaled)
- `#VARVOL=60` with `#NFORM=4` → target volume of 240 Å³ (automatically scaled)

This means you should set `#VARVOL` based on your **per-formula-unit** volume, regardless of how many formula units you plan to generate. The same seed file can then be used with different `#NFORM` values without adjusting `#VARVOL`.


In fact, the scaling is done based on the number of atoms in the unit cell. Consider following example:

```
%BLOCK POSITIONS_FRAC
Si 0 0 0 # Si % NUM=1-3
O 0 0 0 #  O % NUM=1-3
%ENDBLOCK POSITIONS_FRAC

#VARVOL=30
```



### VARVOL Guidelines

The target volume per atom depends on the material type:

| Material Type | Volume (Å³/atom) | Examples |
|--------------|-----------------|----------|
| Dense metals, high pressure | 5-10 | Fe, Al at high P |
| Typical metals | 10-15 | Cu, Ag at ambient |
| Oxides, semiconductors | 10-15 | SiO₂, TiO₂ |
| Molecular crystals | 15-25 | Ice, organic crystals |
| Open frameworks | 20-30+ | Zeolites, MOFs |

If unsure, start with **15 Å³/atom** and adjust based on whether structures build successfully.

> **Note:** Without `#VARVOL`, buildcell uses a default 5% volume variation around the cell volume defined in LATTICE_CART/LATTICE_ABC.

### MINSEP from Known Bond Lengths

Use typical bond lengths as a guide, setting MINSEP slightly shorter than equilibrium values to allow flexibility:

| Bond Type | Typical Length (Å) | Suggested MINSEP |
|-----------|-------------------|------------------|
| C-C (single) | 1.54 | 1.3 |
| C-C (double) | 1.34 | 1.2 |
| C-C (aromatic) | 1.40 | 1.2 |
| C-O | 1.43 | 1.2 |
| C-H | 1.09 | 0.9 |
| O-H | 0.96 | 0.8 |
| Si-O | 1.61 | 1.5 |
| Ti-O | 1.95 | 1.8 |
| Metal-O (general) | 1.8-2.4 | 1.6-2.2 |
| O-O (non-bonded) | 2.8+ | 2.5 |
| Metal-Metal | varies | 2.0-3.0 |

For non-bonded pairs (like O-O in oxides), use larger MINSEP values (2.5+ Å) to prevent unphysical clustering.

### Quick Reference for Common Systems

**Carbon (high pressure):**
```
#VARVOL=5
#MINSEP=1.3
```

**SiO₂:**
```
#VARVOL=45
#MINSEP=1.0 Si-Si=3.00 Si-O=1.60 O-O=2.58
```

**Metal oxide (e.g., TiO₂):**
```
#VARVOL=35
#MINSEP=1.5 Ti-Ti=2.8 Ti-O=1.8 O-O=2.5
```

**Perovskite (e.g., SrTiO₃):**
```
#VARVOL=60
#MINSEP=1.5 Sr-Sr=3.5 Sr-Ti=3.0 Sr-O=2.4 Ti-Ti=2.8 Ti-O=1.8 O-O=2.5
```
