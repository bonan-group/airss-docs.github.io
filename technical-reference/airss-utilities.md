---
title: "AIRSS Utilities"
layout: single
classes: wide
lang: en
lang-ref: airss-utilities
sidebar:
  nav: "docs"
toc: true
toc_sticky: true
---

This page provides details of the various utilities included in the AIRSS package. They are grouped by function: core workflow scripts, structure generation, structure analysis, relaxation wrappers, file management, job management, and format conversion.

## Core Workflow

### airss.pl

The main AIRSS script. This Perl program performs ab initio random structure searching by repeatedly generating random structures with `buildcell` and relaxing them with a chosen code. The default relaxation engine is CASTEP, but many alternatives are supported.

```console
$ airss.pl -seed Carbon -max 100 -press 10
```

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-seed s` | (required) | Seed name; expects `<seed>.cell` (and `<seed>.param` for CASTEP) |
| `-max n` | 1000000 | Maximum number of structures to generate |
| `-pressure f` | 0.0 | External pressure in GPa |
| `-prand` | false | Randomise pressure between `-pmin` and `-pressure` |
| `-pmin f` | 0.0 | Minimum pressure when randomising |
| `-mpinp n` | 1 | Number of MPI cores per calculation |
| `-steps n` | 400 | Maximum geometry optimisation steps |
| `-build` | false | Build structures only (no relaxation) |
| `-nbest n` | 1 | Best of n random structures kept per trial |
| `-best` | false | Only keep the best structure for each composition |
| `-track` | false | Keep track of good structures during relax-and-shake |
| `-keep` | false | Keep all intermediate files |
| `-pack` | false | Concatenate results into a `.res.tar` archive |
| `-harvest` | false | Collect intermediate structures as `.res` files |
| `-sim f` | 0.0 | Threshold for structure similarity filtering (0 = off) |
| `-symm f` | 0.0 | Symmetrise structures on-the-fly (0 = off) |
| `-nosymm` | false | Disable symmetry detection in output |
| `-num n` | 0 | Number of shake trials from best structure (RASH mode) |
| `-amp f` | -1.5 | Amplitude of ionic displacement for shaking |
| `-camp f` | 0.0 | Amplitude of cell move for shaking |
| `-mode` | false | Push structures along low-lying vibrational modes |
| `-minmode n` | 4 | Lowest vibrational mode to use |
| `-dos` | false | Calculate density of states at the Fermi level |
| `-workdir s` | `.` | Working directory for calculations |
| `-cluster` | false | Use cluster settings for the symmetry finder |
| `-slab` | false | Use slab settings |
| `-exec s` | (auto) | Custom executable name |
| `-launch s` | `mpirun -np` | Custom MPI launcher command |

**Supported relaxation codes** (select with flag):

| Flag | Code |
|------|------|
| (default) | CASTEP |
| `-pp3` | pp3 (built-in pair potentials, 3D periodic) |
| `-pp0` | pp0 (built-in pair potentials, 0D) |
| `-gosh` | gosh (built-in pair potentials, 0D) |
| `-gulp` | GULP |
| `-vasp` | VASP |
| `-qe` | Quantum Espresso |
| `-lammps` | LAMMPS |
| `-gap` | GAP (via QUIP/QUIPPY/ASE) |
| `-psi4` | Psi4 |
| `-python` | Custom Python script (`relax.py`) |
| `-repose` | Repose (data-derived potentials) |
| `-ramble` | Ramble (data-derived potentials, dynamics) |
| `-matsim` | MatterSim (machine learning potential) |
| `-castep` | Use CASTEP in addition to the chosen alternative |

The script creates a `stop` file check at the beginning of each iteration. To stop a running search gracefully, create a file called `stop` in the working directory:

```console
$ touch stop
```

See [AIRSS Examples](../tutorials) for detailed usage instructions.

### crud.pl

The CASTEP Run Daemon (CRUD). A Perl script for high-throughput batch relaxation of pre-existing structures. Structures (as `.res` files) are placed in a `hopper/` subdirectory. The script picks them up, relaxes them, and sorts the results.

```console
$ crud.pl -mpinp 4
```

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-exec s` | `castep` | Executable to use |
| `-launch s` | `mpirun -np` | MPI launcher command |
| `-mpinp n` | 0 | Number of MPI cores (0 = serial) |
| `-repose` | false | Use repose for relaxation |
| `-ramble` | false | Use ramble for dynamics |
| `-keep` | false | Keep all output files (default: clean up) |
| `-nostop` | false | Keep running as a daemon after hopper is empty |
| `-cycle` | false | Retry failed runs (move back to hopper) |
| `-pack` | false | Work with packed `.res.tar` archives |
| `-num n` | 1000 | Maximum number of structures to unpack at once |
| `-workdir s` | `.` | Working directory |

**Workflow:**

1. Place `.res` files in `./hopper/`
2. Place `<seed>.param` (and optionally `<seed>.cell` with non-structural settings) in the current directory
3. Run `crud.pl`
4. Successful results go to `./good_castep/`
5. Failed calculations go to `./bad_castep/`
6. Create a file named `STOP_CRUD` to halt the daemon

### run.pl

A Perl script for running batch CASTEP calculations on existing structure files. Useful for "polishing" results, _i.e._ re-running low-energy structures at higher calculation accuracy for publication.

```console
$ run.pl -mpinp 4
```

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-mpinp n` | 0 | Number of MPI cores (0 = serial) |
| `-conventional` | false | Use conventional cell for CIF inputs (no reduction) |
| `-keep` | false | Keep all output files |
| `-check` | false | Check crystal structure files for validity (CIF only) |

The script processes all `*-*.res` and `*-*.cif` files in the current directory. It reads `<root>.param` and `<root>.cell` for computational settings (where `<root>` is the part of the filename before the first hyphen). A `jobs.txt` file is created to track processed structures and prevent duplicate runs. Failed calculations are moved to `./bad_castep/`.

---

## Structure Generation

### buildcell

This Fortran program reads annotated CASTEP `.cell` files and generates random "sensible" structures from them. It is the core structure generation engine of AIRSS. The type of randomness introduced can be controlled through hash-tagged directives (which are treated as comments and ignored by CASTEP).

```console
$ buildcell < seed.cell > random.cell 2>/dev/null
```

`buildcell` reads from standard input and writes the generated structure to standard output. Diagnostic information is sent to standard error.

For full documentation of all buildcell tags and options, see the [Buildcell Manual](../buildcell-manual/).

### gencell

A bash script that generates a recommended set of CASTEP `.cell` and `.param` files from a supplied unit cell volume and atomic composition. It is strongly recommended as a starting point for most projects.

```console
$ gencell <volume> <units> n x [<species> <number>] ...
```

**Arguments:**

- `volume` -- Target volume per formula unit in cubic Angstroms
- `units` -- Number of formula units
- `species` / `number` -- Pairs of element symbol and atom count

**Example:**

```console
$ gencell 40 1 Si 2 O 4
```

This creates `Si2O4.cell` and `Si2O4.param` with sensible defaults: PBEsol functional, QC5 on-the-fly pseudopotentials, a k-point spacing of 0.07, and buildcell tags for a moderate structure search. The generated `.cell` file includes several commented-out (`##`) alternative options that can be enabled by removing one `#`. A `.par` file for ramble/repose hot-AIRSS is also generated.

### genqe

Similar to `gencell`, but generates input files for Quantum Espresso instead of CASTEP. Creates both a `<seed>.cell` file (for `buildcell`) and a `<seed>.qe` file (for `qe_relax`).

```console
$ genqe <volume> <units> n x [<species> <number>] ...
```

After running, you must manually edit the `.qe` file to specify pseudopotential filenames, atomic masses, and appropriate wavefunction/charge-density cutoffs (`ecutwfc` and `ecutrho`).

---

## Structure Analysis

### cryan

A general-purpose Fortran program for analysing large collections of structure data. Structures are read from standard input in SHELX `.res` format.

```console
$ cat *.res | cryan -s
$ gunzip -c lots.res.gz | cryan -f H2O -r
$ find . -name "*.res" | xargs cat | cryan -m
$ cat H2O-P21c.res | cryan -g
```

Experience suggests that `cryan` is suitable for the analysis of up to about 100,000 structures. Other techniques are required for larger data sets. Type `cryan` with no arguments to see the full list of options.

### ca

A convenient bash wrapper for the `cryan` tool that automatically finds and concatenates all `.res` files in the current directory (including compressed `.res.xz` and archived `.res.tar` files) and pipes them to `cryan`. Uses the same command-line options as `cryan`.

```console
$ ca -r                # Rank structures
$ ca -s                # Summary statistics
$ ca -u 0.01           # Unique structures within tolerance
$ ca -R -r             # Recursively search subdirectories
```

The `-R` flag enables recursive searching through subdirectories, rather than the default of only looking in the current directory.

### symm

Finds the space group of a structure. Works with `.res` files and uses `cabal` for cell standardisation internally.

```console
$ symm test            # Crystal symmetry of test.res
$ symm -cl test        # Cluster (molecular) symmetry of test.res
```

The crystal symmetry mode uses `cabal` to determine the space group. The cluster mode (`-cl`) requires the external `symmol` program and determines the Schoenflies point group.

### csymm

Determines the Schoenflies point group symmetry of a cluster (molecule) from an XYZ file. Requires the external `symmol` program.

```console
$ csymm structure.xyz
```

### comp2minsep

Looks up recommended minimum separations and target volumes for a given composition, based on known structures in the current database. This is used internally by `buildcell` when `#MINSEP=AUTO` is specified.

```console
$ comp2minsep SiO2
```

Outputs `#MINSEP` and `#TARGVOL` tags suitable for use in a `.cell` file.

### fm

A structure filtering tool that removes near-duplicate structures based on pressure, volume, and enthalpy values. Reads `cryan -r` formatted output from standard input.

```console
$ ca -r | fm 0.1              # Filter with threshold 0.1
$ ca -r | fm 0.01 --delete    # Filter and delete duplicate .res files
```

**Arguments:**

- `threshold` (optional, default 0.1) -- Tolerance for matching pressure, volume, and enthalpy
- `--delete` -- Actually remove the duplicate `.res` files from disk (prompts for confirmation)

### mc

A small utility to create a reference convergence file (`ref.conv`). Useful for setting a known energy baseline when monitoring convergence.

```console
$ mc <energy_per_fu> <n_fu>
```

Creates `ref.conv` containing the total energy (energy/fu times number of formula units) at 0 and 100 iterations.

---

## Structure Conversion

### cabal

A versatile structure format conversion tool. Also performs Niggli reduction when input and output formats are the same.

```console
$ cabal in out < seed.in > seed.out
```

**Supported formats:**

| Format | Code | Read | Write |
|--------|------|------|-------|
| CASTEP input | `castep` | Yes | No |
| CASTEP cell | `cell` | Yes | Yes |
| SHELX | `shx`/`res` | Yes | Yes |
| GULP | `gulp` | No | Yes |
| CIF | `cif` | No | Yes |
| Psi4 | `psi4` | No | Yes |
| XTL | `xtl` | Yes | Yes |
| XYZ | `xyz` | Yes | Yes |
| Extended XYZ | `xyze` | Yes | Yes |
| POSCAR | `poscar` | Yes | Yes |
| LAMMPS conf | `conf` | No | Yes |
| Quantum Espresso | `qe` | No | Yes |

**Cell standardisation:** When the input and output formats are the same, `cabal` performs cell reduction. An optional third argument controls the tolerance:

```console
$ cabal cell cell < input.cell > output.cell       # Niggli reduction (tolerance 0)
$ cabal res res 0.1 < input.res > output.res        # Conventional cell (tolerance 0.1)
$ cabal res res -0.1 < input.res > output.res       # Primitive cell (tolerance -0.1)
$ cabal res res -0.01 < input.res > output.res      # Symmetrised cell (tolerance -0.01)
```

See also: `cif2res` in the [External Utilities](../external-utilities/) page.

### castep2res

Converts CASTEP output to a SHELX `.res` file, extracting key results (pressure, volume, enthalpy, spin) into the TITL line and computational parameters into REM lines.

```console
$ castep2res seed > seed.res
$ castep2res -cl seed > seed.res    # Cluster mode (molecular symmetry)
```

This script collects information from `<seed>.castep`, `<seed>-out.cell`, `<seed>.cell`, and optionally `<seed>.magres` and `<seed>.odo` files. The TITL line contains: name, pressure, volume, enthalpy, spin, |spin|, delta, number of atoms, and space group. REM lines record the run date, code version, functional, cutoff, k-point grid, pseudopotentials, runtime, and the original AIRSS command line.

### casteps2res

Batch-converts all `.castep` files in the current directory to `.res` files using `castep2res`. Uses GNU `parallel` if available for faster processing.

```console
$ casteps2res
```

### castepsplit2res

Splits a multi-structure `.castep` file (e.g., from a geometry optimisation trajectory) into individual `.res` files. Reads from standard input.

```console
$ castepsplit2res < seed.castep > harvest.res
```

This is used by `airss.pl -harvest` to extract intermediate structures from geometry optimisation runs.

---

## Cell Transformations

### conv

Converts all `.res` files in the current directory to conventional cells. Uses `cabal res res 0.1` internally. Supports GNU `parallel` for faster processing.

```console
$ conv
```

See also: `niggli` and `prim`.

### niggli

Performs Niggli reduction on all `.res` files in the current directory. Uses `cabal res res 0` internally. Supports GNU `parallel` for faster processing.

```console
$ niggli
```

See also: `conv` and `prim`.

### prim

Converts all `.res` files in the current directory to primitive cells. Uses `cabal res res -0.1` internally. Supports GNU `parallel` for faster processing.

```console
$ prim
```

See also: `conv` and `niggli`.

---

## Relaxation Wrappers

Each relaxation wrapper script interfaces between `airss.pl` and a specific simulation code. They all follow a common pattern: convert the AIRSS `.cell` file to the code's native format, run the relaxation, construct a fake `.castep` file with standardised output, and return pressure, enthalpy, and volume on standard output.

### castep_relax

Performs a self-consistent geometry optimisation using CASTEP. This is the default relaxation engine for `airss.pl`.

```console
$ castep_relax <max_iter> <executable> <similarity> <symmetry> <seed>
```

The script performs the optimisation in stages: first three short rough runs (3 steps each), then a full optimisation until convergence. If `similarity` is non-zero, previously seen structures are detected and skipped. If `symmetry` is non-zero, structures are symmetrised on-the-fly after each stage. A negative `max_iter` performs a single-point energy calculation.

### gulp_relax

Performs a geometry optimisation using GULP.

```console
$ gulp_relax <executable> <cluster_flag> <pressure> <seed>
```

Where `cluster_flag` is 1 for cluster (constant volume) or 0 for periodic (constant pressure) calculations. The GULP library file is expected as `<root>.lib`.

### vasp_relax

Performs a geometry optimisation using VASP.

```console
$ vasp_relax <executable> <seed>
```

Requires `<seed>.POTCAR` and `<seed>.INCAR` files. The script runs 4 consecutive VASP optimisations, copying CONTCAR to POSCAR between each run.

### qe_relax

Performs a self-consistent geometry optimisation using Quantum Espresso.

```console
$ qe_relax <max_iter> <executable> <pressure_GPa> <similarity> <symmetry> <seed>
```

Requires a `<seed>.qe` input file with QE settings. The script follows a similar multi-stage approach to `castep_relax`: three short rough runs followed by a full optimisation. Handles unit conversions between QE's native Rydberg/kbar units and AIRSS's eV/GPa conventions.

### pp3_relax

Performs a geometry optimisation using `pp3`, a built-in pair potential code for 3D periodic systems.

```console
$ pp3_relax <executable> <seed>
```

Reads potential parameters from `<seed>.par` (command-line options) and `<seed>.pp`/`<seed>.ppp`/`<seed>.ddp` (potential data files).

### pp0_relax

Performs a geometry optimisation using `pp0` or `gosh`, built-in pair potential codes for 0D (cluster/molecular) systems.

```console
$ pp0_relax <executable> <seed>
```

Converts the cell to XYZ format, runs the optimisation, and converts back.

### lammps_relax

Performs a geometry optimisation using LAMMPS.

```console
$ lammps_relax <executable> <pressure> <seed>
```

Generates a LAMMPS input script with multiple minimisation stages of decreasing tolerance. Requires a `<seed>.pp` file defining the interatomic potentials in LAMMPS format. Pressure is specified in GPa and converted to atmospheres internally.

> **Note:** This relaxation wrapper is not currently recommended due to issues with structural optimisation.

### gap_relax

Performs a geometry optimisation using GAP (Gaussian Approximation Potentials) through QUIP/QUIPPY/ASE.

```console
$ gap_relax <python_script> <seed>
```

The executable argument should be an ASE relaxation script (e.g., `ase_relax_cell.py`) that reads a `.cell` file from standard input and writes the relaxed structure.

### psi4_relax

Performs a geometry optimisation using Psi4 (quantum chemistry).

```console
$ psi4_relax <executable> <seed>
```

Converts the cell to Psi4 input format, runs a PBE/cc-pvdz optimisation, and converts back. Suitable for molecular (cluster) calculations.

> **Note:** This relaxation wrapper is not currently recommended due to issues with structural optimisation.

### python_relax

Performs a geometry optimisation using a user-provided Python script.

```console
$ python_relax <command> <pressure> <seed>
```

Converts the cell to extended XYZ format, runs the specified Python command with `<pressure>` and `<seed>` as arguments, and expects the script to produce `<seed>-out.xyz` and `<seed>.castep` output files. This provides a flexible interface for custom relaxation workflows.

### repose_relax

Performs a geometry optimisation using `repose`, a data-derived interatomic potential code.

```console
$ repose_relax <executable> <n_omp_threads> <seed>
```

Reads potential data from `<seed>.ddp` or `<seed>.eddp` files and command-line parameters from `<seed>.par`.

### matsim_relax

Performs a geometry optimisation using MatterSim, a machine-learning interatomic potential.

```console
$ matsim_relax <executable> <n_omp_threads> <pressure> <seed>
```

The executable is `msrelax`, a Python script that uses the MatterSim universal potential (`MatterSim-v1.0.0-5M.pth`) with ASE's FIRE optimiser. Supports both fixed-cell and variable-cell relaxation (controlled by `FIX_ALL_CELL` in the `.cell` file).

---

## Res File Management

### rescat

Extracts and displays a specific structure from `.res` files, compressed `.res.xz` files, or `.res.tar` archives.

```console
$ rescat seed-name-123         # Display a specific structure
$ rescat seed1 seed2 seed3     # Display multiple structures
```

Searches through individual `.res` files first, then falls back to packed/compressed archives.

### resdel

Marks structures as removed by replacing `TITL` with `REMOVED` in their `.res` files. The structures are not physically deleted but will be ignored by analysis tools.

```console
$ resdel seed-name-123         # Mark one structure as removed
$ resdel seed1 seed2 seed3     # Mark multiple structures
```

Works on both individual `.res` files and structures within concatenated files.

### reshake

Generates "shaken" variants of existing structures for re-relaxation. Takes low-energy structures, applies random perturbations to atomic positions and cell parameters, and optionally creates supercells.

```console
$ reshake <posamp> <cellamp> <maxatoms> <nstructures> <seed>
```

**Arguments:**

- `posamp` -- Maximum amplitude of position perturbation (Angstroms)
- `cellamp` -- Maximum amplitude of cell perturbation
- `maxatoms` -- Maximum atoms per structure (controls supercell size)
- `nstructures` -- Number of shaken variants per input structure
- `seed` -- Seed name (processes all `<seed>-*.res` files)

Output structures are placed in a `./shook/` subdirectory. The first structure for each input is always the unshaken original (possibly as a supercell).

### resname

Renames all `.res` files in the current directory using a descriptive naming scheme based on their composition, space group, and energy ranking. Requires GNU `parallel`.

```console
$ resname
```

The new names follow the pattern `<root><formula>-<spacegroup>-<energy>-<random>.res`.

### name

Similar to `resname`, but copies files to a `./named/` subdirectory with names based on composition information and a user-specified prefix.

```console
$ name <prefix>
```

### respack

Packs all files matching `<seed>-*.*` into a tar archive and removes the originals.

```console
$ respack <seed>
```

Creates `<seed>.res.tar`.

### resunpack

Unpacks a `.res.tar` archive.

```console
$ resunpack <seed>
```

Extracts files from `<seed>.res.tar`.

### resplit

Splits a concatenated multi-structure `.res` stream (from standard input) into individual `.res` files, organised into directories by seed name.

```console
$ cat packed.res | resplit
$ ca -r | resplit
```

### resview

Displays a structure in standardised form (conventional cell with Niggli reduction). Searches for the structure in `.res` files, `.castep` files, and `.cell` files.

```console
$ resview seed-name-123
```

### ress2xyz

Batch-converts all `.res` files in the current directory to XYZ format.

```console
$ ress2xyz
```

### pack2tar

Splits a concatenated `.res` stream from standard input into individual `.res.tar` archives, grouped by seed name.

```console
$ cat many_structures.res | pack2tar
```

---

## Job Management

### spawn

Submits multiple parallel jobs to remote machines via SSH. Reads the machine list from `~/.spawn`.

```console
$ spawn airss.pl -seed Carbon
```

**The `~/.spawn` file format:**

```
node1 slots=8 root=
node2 slots=8 root=
node3 slots=12 root=
```

- `slots` -- Number of available cores on the node
- `root` -- Path prefix (if remote filesystem mount differs from local)

The script divides the available slots by the `-mpinp` value (if present in the command) to determine the number of jobs per node. Jobs are launched via SSH and run in the background. PID files (`.spawnpids.*`) are created for later use by `despawn`.

> **Note:** Password-free SSH access to the remote nodes is required. The alternative to `spawn` is to use the queuing system of a multiuser cluster.

### spawn-slow

Identical to `spawn`, but launches jobs sequentially (one node at a time) rather than in parallel. This prevents race conditions when launching `run.pl` or `crud.pl`, which can try to grab the same file if started simultaneously.

```console
$ spawn-slow crud.pl -mpinp 4
```

### despawn

Halts all remotely spawned jobs in a controlled manner by reading the `.spawnpids.*` files created by `spawn`/`spawn-slow` and sending kill signals to the corresponding process groups.

```console
$ despawn
```

This is the recommended way to stop spawned jobs.

### stopairss

An aggressive script to kill all spawned jobs. Connects to every node in `~/.spawn` and kills all `buildcell`, `perl`, and `castep` processes owned by you.

```console
$ stopairss
```

> **Warning:** This kills _all_ your jobs on the remote nodes. Use `despawn` in preference.

### remote

A more advanced job submission tool for running AIRSS searches on remote servers. Uses `~/.remote` for server configuration and handles file synchronisation via `rsync`.

```console
$ remote spawn airss.pl -seed Carbon
```

**The `~/.remote` file format:**

```
server1 workspace=/scratch/user
server2 workspace=/tmp/airss
```

The script creates a temporary working directory on each remote server, syncs input files, launches the specified command, and creates helper scripts for monitoring:
- `./refresh` -- Sync results back from remote
- `./status` -- Check the status of remote jobs
- `./finish` -- Collect results and stop remote jobs
- `./delete` -- Remove remote working directories
- `./sendstop` -- Send a stop signal to remote jobs

A `lock` file prevents multiple simultaneous submissions.

---

## Pressure Management

### press

Sets up a series of calculations at different pressures. Creates input files and a hopper directory for use with `crud.pl`.

```console
$ press <p_min> <p_step> <p_num> <seed>
```

**Example:**

```console
$ press 0 10 5 Carbon     # Creates files for 0, 10, 20, 30, 40 GPa
```

For each pressure, creates `<seed>_<pressure>.cell` (with the appropriate `EXTERNAL_PRESSURE` block), copies parameter files, and places existing `.res` files into `./hopper/` with pressure-labelled names.

---

## Format Conversion

### cifsplit

Splits a concatenated multi-structure CIF file (from standard input) into individual `.cif` files, named by their CSD database codes.

```console
$ cifsplit < many_structures.cif
```

### xyzs2res

Batch-converts all `.xyz` files in the current directory to `.res` format. Requires a vacuum padding argument.

```console
$ xyzs2res 10.0          # 10 Angstrom padding for periodic box
```

### xyzesplit2res

Splits a concatenated extended XYZ file (from standard input) into individual `.res` structures.

```console
$ xyzesplit2res < trajectory.xyze
```

### pdbs2res

Batch-converts all `.pdb` files in the current directory to `.res` format. Requires OpenBabel (`obabel`) for the PDB-to-XYZ conversion and a vacuum padding argument.

```console
$ pdbs2res 10.0
```

### lammps2cell

A Perl script that converts LAMMPS output files (`.conf`, `.lammps`, `.lammpstrj`) to a CASTEP `.cell` file. Used internally by `lammps_relax`.

```console
$ lammps2cell <seed> > output.cell
```

---

## Housekeeping

### tidy.pl

Cleans up the output of an AIRSS run. Removes incomplete calculations, compresses archives, and consolidates `.res` files and `.res.tar` archives.

```console
$ tidy.pl
```

The script prompts for confirmation before deleting files. It will refuse to run if a `jobs.txt` file is present (indicating `run.pl` is active). Specifically, it:

1. Removes empty or incomplete `.res` files (and their associated files)
2. Compresses individual `.res.tar` files with `xz`
3. Concatenates per-process `.res.tar` archives into a single `<root>.res.tar`
4. Packs individual `.res` files into the tar archive
5. Removes temporary files, trash directories, and spawn PID files

### check_airss

A diagnostic script that checks whether all required and optional AIRSS dependencies are installed, and runs a set of basic tests.

```console
$ check_airss
```

Checks for:
- **Essential:** `airss.pl`, `run.pl`, `crud.pl`, `castep2res`, `buildcell`, `cryan`, `pp3`, `cabal`, `symmol`, `md5sum`, `bc`
- **Recommended:** `castep`, `optados`, `qhull`/`qconvex`, `xmgrace`, `Rscript`, `obabel`
- **Optional:** `pw.x` (QE), `vasp`, `gulp`, `cif2cell`
- **Very optional:** `lammps`, `hull`, `off_util`

Also checks for CASTEP pseudopotentials (`$PSPOT_DIR`) and the `~/.spawn` file.

### qb

A quick-build utility. Generates one random structure from a seed and opens it for viewing.

```console
$ qb <seed>
```

Runs `buildcell` on `<seed>.cell`, converts the output to `.res` format, and opens it with the system viewer (macOS `open` command).

### airss_version

Prints the current AIRSS version string.

```console
$ airss_version
```
