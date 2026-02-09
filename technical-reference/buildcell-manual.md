---
title: "Buildcell Manual"
layout: single
classes: wide
lang: en
lang-ref: buildcell-manual
sidebar:
  nav: "docs"
toc: true
toc_sticky: true
---

The construction of reasonable, or _sensible_, random structures is central to AIRSS. The Fortran `buildcell` tool is provided in the AIRSS package for this purpose. It can build structures from scratch, or modify structures specified using the CASTEP `.cell` format. The random structures generated are output in the CASTEP `.cell` format.

`buildcell` reads from standard input (`stdin`), and writes to standard output (`stdout`). Additional information is reported to standard error (`stderr`).

## CASTEP Cell File Format

Buildcell uses the CASTEP `.cell` file format for both input and output. Understanding this format is essential for writing buildcell input files.

### Basic Structure

A CASTEP cell file consists of:
- **Block sections** enclosed by `%BLOCK name` and `%ENDBLOCK name`
- **Keywords** specified as `KEYWORD = value` or `KEYWORD : value`
- **Comments** starting with `#`, `!`, or `;`

> **Important:** In standard CASTEP, lines beginning with `#` are treated as comments and ignored. However, `buildcell` repurposes the `#` character to specify generation tags. This means `#MINSEP=1.5` is a comment to CASTEP but an instruction to buildcell.

### Lattice Specification

The unit cell can be specified in two ways:

**Cartesian vectors** (`LATTICE_CART`):
```
%BLOCK LATTICE_CART
ang              ! Optional: units (ang, bohr). Default: ang
  4.0   0.0   0.0
  0.0   4.0   0.0
  0.0   0.0   4.0
%ENDBLOCK LATTICE_CART
```

**Lattice parameters** (`LATTICE_ABC`):
```
%BLOCK LATTICE_ABC
ang              ! Optional: units
  4.0  4.0  4.0  ! a, b, c lengths
 90.0 90.0 90.0  ! α, β, γ angles in degrees
%ENDBLOCK LATTICE_ABC
```

### Atomic Positions

Positions can be given in fractional or absolute coordinates:

**Fractional coordinates** (`POSITIONS_FRAC`):
```
%BLOCK POSITIONS_FRAC
Si  0.0   0.0   0.0
Si  0.25  0.25  0.25
%ENDBLOCK POSITIONS_FRAC
```

**Absolute coordinates** (`POSITIONS_ABS`):
```
%BLOCK POSITIONS_ABS
ang              ! Optional: units (ang, bohr)
Si  0.0  0.0  0.0
Si  1.0  1.0  1.0
%ENDBLOCK POSITIONS_ABS
```

### Comments and Buildcell Tags

The comment character serves different purposes:

| Syntax | CASTEP | Buildcell |
|--------|--------|-----------|
| `# text` | Ignored (comment) | **Parsed as generation tag** if valid format |
| `## text` | Ignored (comment) | Ignored (comment) |
| `! text` | Ignored (comment) | **Not stripped** - avoid using with buildcell tags |

This design allows buildcell input files to remain valid CASTEP cell files. When CASTEP reads a buildcell-enhanced file, it simply ignores the `#TAG=value` lines as comments.

> **Note:** Unlike CASTEP, buildcell does **not** recognize `!` as a comment character. Text after `!` is passed to the parser and may cause errors. Use `##` for comments in buildcell input files.

**Example showing comment styles:**
```
%BLOCK POSITIONS_FRAC
## This is a comment (ignored by both CASTEP and buildcell)
Si 0.0 0.0 0.0  # Si1 % NUM=4
%ENDBLOCK POSITIONS_FRAC

## Use double-hash for comments in buildcell files
#MINSEP=1.5
##MINSEP=2.0    ## This line is ignored by buildcell
```

### Species and Labels

Each position line has the format:
```
Species  x  y  z  [# label] [% per-atom-tags]
```

- **Species**: Element symbol (e.g., `Si`, `O`, `Fe`)
- **Coordinates**: Three numbers (fractional or absolute depending on block type)
- **Label** (optional): After `#`, used for identification
- **Per-atom tags** (optional): After `%`, control per-atom generation behavior

---

## Basic Example

```console
$ cat Al.cell

%BLOCK LATTICE_CART
2 0 0
0 2 0
0 0 2
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Al 0.0 0.0 0.0 # Al1 % NUM=8
%ENDBLOCK POSITIONS_FRAC

#MINSEP=1.5

$ buildcell < Al.cell > Al-rand.cell
```

## Tag Format Reference

Buildcell tags control structure generation. They are specified using the `#` prefix in the input file.

### Tag Formats

| Format | Example | Description |
|--------|---------|-------------|
| Boolean flag | `#CLUSTER` | Presence enables the feature |
| Single value | `#MINSEP=1.5` | Sets a parameter to a specific value |
| Range | `#SYMMOPS=2-4` | Randomly selects from range (inclusive) |
| List | `#SYMMNO=1,2,5-10,15` | Randomly selects from list/ranges |
| Species pairs | `#MINSEP=1.5 Si-O=1.8 O-O=2.0` | Default value with species-specific overrides |
| Approximate | `#SYMMOPS=~4` | Prefix `~` allows approximate matching |

### Per-Atom Tags

Per-atom tags appear after the `%` symbol in position lines:

```
Al 0.0 0.0 0.0 # label % NUM=8 POSAMP=2.0 FIX
```

### Comments

Lines beginning with `##` are treated as comments and ignored.

---

## Global Tags by Category

### Cell Geometry

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `FIX` | flag | false | Fix the entire unit cell. Place in LATTICE block. |
| `ABFIX` | flag | false | Fix the a- and b-axes. Place in LATTICE block. |
| `CFIX` | flag | false | Fix the c-axis only. Place in LATTICE block. |
| `CELLAMP` | float | -1 | Amplitude for random variation of supplied cell. Negative = generate from scratch. |
| `CELLCON` | 6 floats | none | Cell constraints vector (a, b, c, α, β, γ). Use -1 for free, specific value to fix. |
| `SUPERCELL` | int(s) | identity | Supercell transformation. 1 value = isotropic, 3 = diagonal, 9 = full matrix. |
| `SYSTEM` | string | none | Enforce crystal system: Tric, Mono, Orth, Tetr, Hexa, Rhom, Cubi. Conflicts with CELLCON. |
| `COMPACT` | flag | auto | Force Niggli reduction. Default true unless cell is fixed or cluster. |
| `NOCOMPACT` | flag | false | Disable Niggli reduction. |
| `CONS` | float | 0.4 | Cell shape constraint (0 = total freedom, 1 = cubic only). Controls aspect ratio. |
| `ACONS` | float | 0.5 | Angle constraint. Rejects cells that are too flat. Higher values favor 3D cells. |

**Examples:**

```
#CELLCON=-1 -1 -1 90 90 90   ! Cubic cell (a,b,c free, angles fixed to 90)
#CELLCON=-1 -1 -1 -1 -1 -1   ! Rhombohedral cell (all equal)
#SYSTEM=Cubi                  ! Enforce cubic system
#SUPERCELL=2 2 1              ! 2x2x1 supercell
```

### Volume and Packing

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `TARGVOL` | float or range | auto | Target volume per formula unit (ų). Supports ranges e.g., `40-50`. |
| `VOL` | float or range | auto | Target cell volume (ų). Alternative to TARGVOL. Supports ranges. |
| `VARVOL` | float | from TARGVOL | Variable target volume. Overrides TARGVOL if set. |
| `PACKING` | float | 0.3401 | Packing fraction for volume estimation. Default is diamond packing. |
| `VACUUM` | float | 0.0 | Vacuum padding in Ångstroms (for clusters, surfaces). |

**Examples:**

```
#TARGVOL=35                   ! 35 ų per formula unit
#TARGVOL=30-40                ! Random volume between 30-40 ų
#VOL=100                      ! Exact cell volume of 100 ų
#VACUUM=10.0                  ! 10 Å vacuum padding
```

### Composition Control

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `SEED` | string | system time | Random number seed for reproducibility. |
| `FORMULA` | string | from atoms | Chemical formula (e.g., Si2O4). Cannot use with SPECIES. |
| `NATOM` | int or range | from seed | Total number of atoms. Supports ranges e.g., `5-10`. |
| `SPECIES` | string | from atoms | Comma-separated species with optional modifiers (e.g., `Fe%NUM=2,O%NUM=3`). Cannot use with POSITIONS block. |
| `NFORM` | int or range | -1 | Number of formula units. -1 (default) = auto-determine based on symmetry. |
| `FOCUS` | int | 0 | Focus on n-ary compositions: 1=elements, 2=binaries, 3=ternaries, etc. |

**Examples:**

```
#FORMULA=Si2O4                ! Silicon dioxide with 4 oxygen
#NATOM=8                      ! Exactly 8 atoms
#NATOM=6-12                   ! 6 to 12 atoms
#SPECIES=Fe,O                 ! Iron and oxygen
#SPECIES=Si%NUM=2,O%NUM=4     ! 2 Si and 4 O atoms
#NFORM=2                      ! 2 formula units
#FOCUS=2                      ! Focus on binary compositions
```

### Symmetry Control

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `SYMM` | string | P1 | Space group by symbol (e.g., `Fm-3m`). Prefix with `~` for approximate. |
| `SYMMNO` | int, range, or list | 1 | Space group number (1-230). Supports ranges and comma-separated lists. |
| `SYMMOPS` | int or range | 1 | Number of symmetry operations. Crystal: 1,2,3,4,6,8,12,16,24,48. Cluster: 1-12,24. |
| `SGRANK` | int | 230 | Maximum allowed space group rank. Lower values restrict to simpler groups. |
| `SYMMORPHIC` | flag | false | Restrict to symmorphic space groups only. |
| `CHIRAL` | flag | false | Restrict to chiral (Sohncke) space groups only. |
| `ADJGEN` | int or range | 0 | Adjust general positions. 0 = maximum general positions, higher = allow special positions. |

**Examples:**

```
#SYMMOPS=4                    ! Exactly 4 symmetry operations
#SYMMOPS=2-8                  ! 2 to 8 operations
#SYMMOPS=~4                   ! Approximately 4 operations
#SYMMNO=225                   ! Fm-3m (space group 225)
#SYMMNO=1-15,75-80            ! Select from these space groups
#SYMM=Fm-3m                   ! Fm-3m by name
#SGRANK=50                    ! Only allow rank ≤ 50
#ADJGEN=0-1                   ! Allow some special positions
```

### Position Generation

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `POSAMP` | float | -1 (or sphere radius) | Position randomization amplitude. Negative = spherical distribution filling cell. |
| `MINAMP` | float | 0.0 | Minimum distance from origin for position randomization. |
| `XAMP` | float | -1 | X-direction amplitude. Negative = use spherical POSAMP. |
| `YAMP` | float | -1 | Y-direction amplitude. Negative = use spherical POSAMP. |
| `ZAMP` | float | -1 | Z-direction amplitude. Negative = use spherical POSAMP. |
| `ANGAMP` | float | -1 | Angular rotation amplitude in degrees for units/molecules. Negative = full rotation (360°). |
| `BREAKAMP` | float | -1 | Random displacement amplitude to break symmetry. Negative = no breaking. |

**Examples:**

```
#POSAMP=2.0                   ! Random positions within 2 Å sphere
#ZAMP=0.5                     ! Only ±0.5 Å variation in Z
#XAMP=1.0                     ! ±1 Å in X
#YAMP=1.0                     ! ±1 Å in Y
#ANGAMP=30                    ! ±30° rotation for molecular units
```

### Separation Constraints

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `MINSEP` | float, range, pairs | auto | Minimum atomic separation. Supports species pairs and AUTO keyword. |
| `OVERLAP` | float | -999.9 | Maximum allowed atomic overlap. Negative = disabled. |
| `COORD` | int or range | -1 | Target coordination number. Supports ranges. -1 = no constraint. |
| `MINBANGLE` | float | 0.0 | Minimum bond angle in degrees for coordination constraints. |
| `MAXBANGLE` | float | 180.0 | Maximum bond angle in degrees for coordination constraints. |
| `RAD` | float | 0.0 | Global hard-sphere radius for all atoms. |
| `NFAILS` | int | 0 | Number of constraint failures to tolerate before rejection. |

**MINSEP Format:**

```
#MINSEP=1.5                           ! Global 1.5 Å separation
#MINSEP=1.5-2.0                       ! Random between 1.5-2.0 Å
#MINSEP=1.5 Si-O=1.6 O-O=2.0          ! Species-specific separations
#MINSEP=1.5-2.0 Si-O=1.6-1.8          ! Ranges for species pairs
#MINSEP=AUTO                          ! Use comp2minsep utility
#MINSEP=AUTOVOL                       ! Auto volume only
```

**Examples:**

```
#MINSEP=1.8 Fe-O=1.9 O-O=2.2
#COORD=4                              ! 4-coordinated
#COORD=3-6                            ! 3 to 6 coordination
#MINBANGLE=90                         ! No angles < 90°
```

### Structure Optimization

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `NOPUSH` | flag | false | Disable PUSH (PUt-and-SHake) optimization algorithm. |
| `PUSHSTEP` | float | 0.25 | Step size for atomic displacement in PUSH algorithm. |
| `PUSHMAX` | int | 100 | Maximum number of PUSH iterations. |
| `RASH` | flag | false | Enable RASH (Relax-And-SHake) algorithm. |
| `RASH_POSAMP` | float | 1.0 | Position amplitude for RASH perturbations. |
| `RASH_ANGAMP` | float | 30.0 | Angular amplitude for RASH perturbations (degrees). |
| `CELLADAPT` | flag | false | Allow cell shape changes during distance constraint optimization. |
| `TIGHT` | flag | false | Enable tight packing mode. |

**Examples:**

```
#OVERLAP=0.1                          ! Allow 0.1 Å overlap
#PUSHMAX=200                          ! More optimization steps
#RASH                                 ! Enable RASH
#RASH_POSAMP=0.5                      ! Smaller RASH perturbations
```

### Geometry Constraints

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `CLUSTER` | flag | false | Build cluster (non-periodic) geometry. |
| `SLAB` | flag | false | Build slab geometry. |
| `SURFACE` | flag | false | Build surface model (implies SLAB). |
| `SPHERE` | float | -1 | Confining sphere radius. Positive = hard wall, negative = attractive potential strength. |
| `CORE` | float | none | Repulsive core radius (inner excluded region). |
| `ELLIPSOID` | 2 floats | none | Confining ellipsoid: radius and aspect ratio parameter (0=sphere, larger=more variation). |
| `PANCAKE` | 2 floats | none | Oblate (flat) ellipsoid: radius and flattening (0=flat, 1=sphere). |
| `CIGAR` | 2 floats | none | Prolate (elongated) ellipsoid: radius and elongation (0=needle, 1=sphere). |
| `CYLINDER` | float | -1 | Confining cylinder radius. Positive = hard wall, negative = attractive line potential. |
| `WIDTH` | float | -999.9 | Slab spacer width in Ångstroms. |
| `SHIFT` | float | 0.0 | Slab position shift from center. |
| `SLACK` | float | 0.0 | Fractional slack factor for bonding constraints (0-1). |
| `AUTOSLACK` | float | none | Initial slack value that auto-increments on retry. |

**Examples:**

```
#CLUSTER                              ! Non-periodic cluster
#SPHERE=5.0                           ! 5 Å confining sphere
#VACUUM=10.0                          ! Add vacuum around cluster
#SLAB                                 ! Slab geometry
#WIDTH=8.0                            ! 8 Å slab thickness
#CYLINDER=3.0                         ! Confining cylinder
#SLACK=0.1                            ! 10% slack on separations
```

### Chemical Constraints

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `PERMUTE` | flag or string | false | Enable species permutation. Optionally specify species: `#PERMUTE=Fe,Co`. |
| `PERMFRAC` | float | 1.0 | Fraction of atoms to permute (0-1). |
| `MOLECULES` | flag | false | Treat atomic groups as rigid molecular units. |
| `FLIP` | flag | false | Randomly mirror (reflect) structural units. |
| `REMOVE` | flag | false | Remove atoms that end up at the same position. |
| `OCTET` | flag | false | Check if valence electrons are multiple of 8. |
| `HOLE` | float | -1 | Create spherical hole of given radius. |
| `HOLEPOS` | 3 floats | random | Position of hole in fractional coordinates. |
| `VACANCIES` | int or int@species | 0 | Number of vacancies to introduce. Use `@species` to specify which atom type. |

**Examples:**

```
#PERMUTE                              ! Permute all atoms
#PERMUTE=Fe,Co                        ! Only permute Fe and Co
#PERMFRAC=0.5                         ! Permute 50% of atoms
#MOLECULES                            ! Treat as molecules
#VACANCIES=2                          ! Remove 2 random atoms
#VACANCIES=1@O                        ! Remove 1 oxygen
#HOLE=2.0                             ! 2 Å hole
#HOLEPOS=0.5 0.5 0.5                  ! Hole at cell center
```

### Control Parameters

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `MAXTIME` | float | 1.0 | Maximum time in seconds for structure generation attempts. |
| `THREE` | float | -999.9 | Three-body hard-sphere potential parameter. (Currently deactivated) |

---

## Per-Atom Tags

Per-atom tags are specified after the `%` symbol in the POSITIONS block. They override global settings for individual atoms or atomic groups.

### Syntax

```
Element x y z # label % TAG1=value TAG2 TAG3=min-max
```

### Position and Amplitude Tags

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `NUM` | int or range | 1 | Number of atoms of this type. |
| `POSAMP` | float | global | Position randomization amplitude. |
| `MINAMP` | float | global | Minimum position amplitude. |
| `XAMP` | float | global | X-direction amplitude. |
| `YAMP` | float | global | Y-direction amplitude. |
| `ZAMP` | float | global | Z-direction amplitude. |
| `ANGAMP` | float | global | Angular amplitude for rotation (degrees). |

### Property Tags

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `VOL` | float | -1 | Per-atom volume (ų). If all atoms have VOL set, used for total volume. |
| `RAD` | float | 0.0 | Per-atom hard-sphere radius. |
| `OCC` | float or fraction | 1.0 | Occupation factor. Supports fractions like `1/3`. |
| `MULT` | int | -1 | Wyckoff multiplicity. Sets OCC = MULT/num_symm. |
| `SPIN` | float | 0.0 | Magnetic spin moment. |
| `COORD` | int or range | global | Per-atom coordination constraint. |

### Constraint Flags

| Tag | Format | Default | Description |
|-----|--------|---------|-------------|
| `FIX` | flag | false | Fix position completely (during generation AND DFT). |
| `NOMOVE` | flag | false | Fix during buildcell only (free in DFT). |
| `PERM` | flag | false | Allow this atom to be permuted. |
| `ADATOM` | flag | false | Add this atom after supercell creation. |
| `ATHOLE` | flag | false | Place at hole position (if HOLE is set). |

### Neighbor Constraints

| Tag | Format | Description |
|-----|--------|-------------|
| `NN+species` | string | Require this species as nearest neighbor. |
| `NN-species` | string | Forbid this species as nearest neighbor. |

### Examples

```
%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0 # Si1 % NUM=4 COORD=4
O  0.5 0.5 0.5 # O1  % NUM=8 POSAMP=1.0
Fe 0.0 0.0 0.0 # Fe1 % NUM=2 SPIN=2.5 FIX
Al 0.25 0.25 0.25 # Al1 % OCC=1/2 MULT=4
C  0.0 0.0 0.0 # mol % NUM=1 ANGAMP=30
H  0.1 0.0 0.0 # mol % NUM=3
%ENDBLOCK POSITIONS_FRAC
```

---

## Complete Tag Reference (Glossary)

| **Tag** | **Type** | **Default** | **Description** |
|:--------|:---------|:------------|:----------------|
| **ABFIX** | flag | false | Fix a- and b-axes. Place in LATTICE block. |
| **ACONS** | float | 0.5 | Angle constraint; rejects flat cells. |
| **ADJGEN** | int/range | 0 | Adjust general positions; 0 = max general, higher = allow special. |
| **ANGAMP** | float | -1 | Angular amplitude (degrees); -1 = full rotation. |
| **AUTOSLACK** | float | none | Initial slack that auto-increments on retry. |
| **BREAKAMP** | float | -1 | Symmetry breaking displacement amplitude. |
| **CELLADAPT** | flag | false | Allow cell shape change during optimization. |
| **CELLAMP** | float | -1 | Cell variation amplitude; -1 = generate from scratch. |
| **CELLCON** | 6 floats | none | Cell constraints (a, b, c, α, β, γ). |
| **CFIX** | flag | false | Fix c-axis. Place in LATTICE block. |
| **CHIRAL** | flag | false | Restrict to chiral (Sohncke) space groups. |
| **CIGAR** | 2 floats | none | Prolate ellipsoid: radius and elongation. |
| **CLUSTER** | flag | false | Build non-periodic cluster. |
| **COMPACT** | flag | auto | Force Niggli cell reduction. |
| **CONS** | float | 0.4 | Cell shape constraint (0 = free, 1 = cubic). |
| **COORD** | int/range | -1 | Coordination constraint; -1 = none. |
| **CORE** | float | none | Repulsive core radius. |
| **CYLINDER** | float | -1 | Confining cylinder radius or attractive potential. |
| **ELLIPSOID** | 2 floats | none | Ellipsoid: radius and aspect variation. |
| **FIX** | flag | false | Fix unit cell. Place in LATTICE block. |
| **FLIP** | flag | false | Randomly mirror structural units. |
| **FOCUS** | int | 0 | Focus on n-ary compositions. |
| **FORMULA** | string | from atoms | Chemical formula. |
| **HOLE** | float | -1 | Spherical hole radius. |
| **HOLEPOS** | 3 floats | random | Hole position (fractional). |
| **MAXBANGLE** | float | 180.0 | Maximum bond angle (degrees). |
| **MAXTIME** | float | 1.0 | Max generation time (seconds). |
| **MINAMP** | float | 0.0 | Minimum position amplitude. |
| **MINBANGLE** | float | 0.0 | Minimum bond angle (degrees). |
| **MINSEP** | float/pairs | auto | Minimum separations. Supports AUTO. |
| **MOLECULES** | flag | false | Treat groups as molecular units. |
| **NATOM** | int/range | from seed | Total number of atoms. |
| **NFAILS** | int | 0 | Constraint failures to tolerate. |
| **NFORM** | int/range | -1 | Number of formula units. -1 = auto. |
| **NOCOMPACT** | flag | false | Disable Niggli reduction. |
| **NOPUSH** | flag | false | Disable PUSH optimization. |
| **OCTET** | flag | false | Check octet rule. |
| **OVERLAP** | float | -999.9 | Maximum allowed overlap. |
| **PACKING** | float | 0.3401 | Packing fraction for volume estimation. |
| **PANCAKE** | 2 floats | none | Oblate ellipsoid: radius and flattening. |
| **PERMFRAC** | float | 1.0 | Fraction of atoms to permute. |
| **PERMUTE** | flag/string | false | Enable permutation or specify species. |
| **POSAMP** | float | -1 | Position randomization amplitude. |
| **PUSHMAX** | int | 100 | Maximum PUSH iterations. |
| **PUSHSTEP** | float | 0.25 | PUSH step size. |
| **RAD** | float | 0.0 | Global atomic radius. |
| **RASH** | flag | false | Enable RASH algorithm. |
| **RASH_ANGAMP** | float | 30.0 | RASH angular amplitude. |
| **RASH_POSAMP** | float | 1.0 | RASH position amplitude. |
| **REMOVE** | flag | false | Remove duplicate atoms. |
| **SEED** | string | time | Random seed. |
| **SGRANK** | int | 230 | Maximum space group rank. |
| **SHIFT** | float | 0.0 | Slab position shift. |
| **SLAB** | flag | false | Build slab geometry. |
| **SLACK** | float | 0.0 | Bonding slack factor. |
| **SPECIES** | string | from atoms | Species list with modifiers. |
| **SPHERE** | float | -1 | Confining sphere radius. |
| **SPIN** | 2 floats | 0 0 | Total spin and modulation. |
| **SUPERCELL** | int(s) | identity | Supercell transformation. |
| **SURFACE** | flag | false | Build surface model. |
| **SYMM** | string | P1 | Space group by name. |
| **SYMMNO** | int/range/list | 1 | Space group number (1-230). |
| **SYMMOPS** | int/range | 1 | Number of symmetry operations. |
| **SYMMORPHIC** | flag | false | Restrict to symmorphic groups. |
| **SYSTEM** | string | none | Crystal system. |
| **TARGVOL** | float/range | auto | Target volume per formula unit. |
| **THREE** | float | -999.9 | Three-body potential (deactivated). |
| **TIGHT** | flag | false | Tight packing mode. |
| **VACANCIES** | int | 0 | Number of vacancies. |
| **VACUUM** | float | 0.0 | Vacuum padding (Å). |
| **VOL** | float/range | auto | Target cell volume (ų). |
| **VARVOL** | float | from TARGVOL | Variable target volume. |
| **WIDTH** | float | -999.9 | Slab spacer width. |
| **XAMP** | float | -1 | X-direction amplitude. |
| **YAMP** | float | -1 | Y-direction amplitude. |
| **ZAMP** | float | -1 | Z-direction amplitude. |

---

## Usage Examples

### Simple Crystal Search

```
%BLOCK LATTICE_CART
2 0 0
0 2 0
0 0 2
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0 # Si1 % NUM=8
%ENDBLOCK POSITIONS_FRAC

#MINSEP=2.0
#TARGVOL=20
#SYMMOPS=2-4
```

### Binary System with Constraints

```
%BLOCK LATTICE_CART
3 0 0
0 3 0
0 0 3
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0 # Si1 % NUM=2 COORD=4
O  0.5 0.5 0.5 # O1  % NUM=4 COORD=2
%ENDBLOCK POSITIONS_FRAC

#MINSEP=1.4 Si-O=1.6 O-O=2.2 Si-Si=2.4
#TARGVOL=35
#SYMMOPS=2-8
#NFORM=1
```

### High-Pressure Binary Search

```
#SPECIES=Fe,Bi
#NATOM=3-6
#FOCUS=2
#SYMMOPS=2-4
#NFORM=1
#ADJGEN=0-1
#SLACK=0.25
#OVERLAP=0.1
#MINSEP=1-3 AUTO
#COMPACT
```

### Cluster Generation

```
%BLOCK POSITIONS_FRAC
Au 0.0 0.0 0.0 # Au1 % NUM=13
%ENDBLOCK POSITIONS_FRAC

#CLUSTER
#SPHERE=4.0
#VACUUM=10.0
#MINSEP=2.5
```

### Surface/Slab Model

```
%BLOCK LATTICE_CART
4 0 0
0 4 0
0 0 15
#CFIX
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Pt 0.0 0.0 0.0 # Pt_bulk % NUM=4 ZAMP=0.1
Pt 0.0 0.0 0.5 # Pt_surf % NUM=2 ZAMP=1.0
%ENDBLOCK POSITIONS_FRAC

#SLAB
#WIDTH=8.0
#VACUUM=7.0
#MINSEP=2.4
```

### Molecular Unit

```
%BLOCK POSITIONS_FRAC
C 0.0 0.0 0.0 # CH4 % ANGAMP=360
H 0.1 0.1 0.1 # CH4
H -0.1 -0.1 0.1 # CH4
H 0.1 -0.1 -0.1 # CH4
H -0.1 0.1 -0.1 # CH4
%ENDBLOCK POSITIONS_FRAC

#MOLECULES
#MINSEP=2.5
#TARGVOL=40
```

### High-Symmetry Search

```
%BLOCK POSITIONS_FRAC
C 0.0 0.0 0.0 # C1 % NUM=8
%ENDBLOCK POSITIONS_FRAC

#SYSTEM=Cubi
#SYMMOPS=24-48
#MINSEP=1.4
#COMPACT
```
