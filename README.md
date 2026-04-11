# Craton

Craton is a **computational chemistry toolkit** for molecular simulation and drug design. It covers the full workflow from quantum chemistry (QM), molecular mechanics (MM), force-field parameterization, and molecular dynamics (MD) to free-energy perturbation (FEP), with support for MongoDB, Gaussian, GROMACS, and other external tools.

## Overview

### Chemical structure and molecule handling
- **Read/write and format conversion**: SMILES, MOL/SDF, MOL2, PDB, etc., via `molecule_create` and `format_convert`
- **Topology and structure analysis**: rings, aromaticity, chirality, rotatable bonds, functional-group recognition, and fragmentation (MolScalpel)
- **Biomolecules**: protein structure preparation, mutations, sequence alignment, and atom mapping
- **Conformations**: conformer search, dihedral/angle scans, RMSD, local minima, and removal of similar conformers

### Quantum chemistry (QM)
- Multi-stage QM pipeline (Q0–Q10): geometry optimization, frequencies, single-point, charges, dihedral/bond-angle scans, etc.
- Integration with Gaussian; methods and basis sets configurable via `Q*QMSetting` in `configure/configure.yaml`

### Molecular mechanics and force fields
- **Atom typing**: PLFF, GAFF2, empirical force fields, and custom typing files
- **Force-field assignment**: bonds, angles, dihedrals, impropers, nonbonded (including charge methods)
- **Force-field fitting**: parameter fitting and validation against energy, forces, Hessian, torsion penalties, etc. (e.g. OpenMM)

### Molecular dynamics and free energy
- **MD types**: vacuum, solution, liquid, protein, complex; standard MD, pull, etc.
- **FEP types**:
  - **RBFE**: relative binding free energy (r_group, core, charge, mutation, mem-rbfe, pep-rbfe, etc.)
  - **ABFE / HFE**: absolute binding free energy, hydration free energy (ahfe, rhfe, alogp, rlogp, rlogs, etc.)
- **Engine**: primarily GROMACS (gmx); HPC and job settings can be configured

### Other
- **System building**: box setup for solution/liquid/complex/membrane (density, concentration, counterions, etc.)
- **Docking**: ligand/protein preparation, pocket analysis, integration with Vina and others
- **Data and database**: MongoDB read/write (compounds, QM data, etc.), PubChem/PDB/UniProt queries
- **Properties**: thermodynamic properties, IC50–free-energy conversion, chemical-space analysis, etc.

## Usage

### Python API (MolXpert)
Use the unified `molxpert` interface in code or Jupyter:

```python
from craton.craton import molxpert, CRATON_CONFIGURE

# Create molecules from SMILES and run structure analysis
mols = molxpert.molecule_create(["CCO", "c1ccccc1"])
mols = molxpert.molecule_structure(mols)

# Force field and MM energy
mols = molxpert.get_force_field(mols)
energies = molxpert.energy(mols, prop="energy")
```

Override defaults via `molxpert.update_configure(usr_config)` or by editing `configure/configure.yaml`.

### Command line
After installation, use the `craton` command and its subcommands:

- **simulation**: submit or configure MD or QM jobs (e.g. `rbfe`, `abfe`, `solution`, `Q0`–`Q10`, `yaml`)
- **analyze**: FEP result analysis, property statistics, chemical space, train/test split, etc.
- **force_field (ff)**: atom typing, force-field assignment, AM1BCC charges, force-field fitting
- **mm**: single-point energy, multipoles, volume/surface area, moment of inertia, etc.
- **stru**: topology analysis, structure parameters, RMSD
- **prepare**: ligand/protein preparation, protonation, PDB/UniProt queries
- **data**: MongoDB import/export
- **tool**: torsion-scan analysis, peptide generation, thermodynamic properties, format conversion, etc.

Example:

```bash
# Run RBFE workflow with a YAML file (provide your own YAML and inputs)
craton simulation yaml -f your_rbfe.yaml

# Force-field fitting
craton ff fitting <fit_type> <input_files> -o <output_directory> ...
```

## Install

```bash
git clone https://github.com/Gewu-Intelligence/craton.git
cd craton
conda env create -f env.yaml
conda activate craton
pip install -e .
git config --local core.autocrlf false
git config --local core.eol lf
```

- If the default Conda channel is unreachable, copy `.condarc` to your home directory to use the Tsinghua mirror.
