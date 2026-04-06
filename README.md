# Borosilicate Glass + ReaxFF Ion Exchange Workflow

This repository contains a minimal workflow to:
1. Generate an initial sodium-borosilicate configuration.
2. Equilibrate it into a glass using a melt-quench ReaxFF simulation.
3. Run a Na↔K ion-exchange simulation from the equilibrated structure.

## Files

- `generate_borosilicate_glass.py` – creates `borosilicate_init.data`.
- `in.equilibrate_borosilicate_reaxff` – melt-quench equilibration.
- `in.ion_exchange_reaxff` – ion exchange on equilibrated glass.
- `borosilicate_init.data` – generated starting structure.

## Requirements

- LAMMPS build with ReaxFF support (`pair_style reaxff` and `fix qeq/reaxff`).
- A ReaxFF parameter file containing **O/Si/B/Na/K** in the order used below,
  saved as `ffield.reax.borosilicate`.

## Run

```bash
python3 generate_borosilicate_glass.py --atoms 3000 --box 42 --out borosilicate_init.data
lmp -in in.equilibrate_borosilicate_reaxff
lmp -in in.ion_exchange_reaxff
```

## Notes

- The generated data file is a random-packed precursor. The structure becomes
  physically meaningful **after** melt-quench equilibration.
- Ion exchange is represented by near-surface Na→K type conversion followed by
  finite-temperature diffusion under ReaxFF. For a more explicit bath model,
  add a salt region (e.g., nitrate melt) with a compatible parameterization.
