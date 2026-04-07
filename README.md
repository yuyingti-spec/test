# ReaxFF Tensile-Failure Workflow for a Hollow Borosilicate Optical Fiber

This repository provides a **LAMMPS ReaxFF workflow** to:
1. define a **real-size hollow optical-fiber geometry**,
2. create an atomistic shell model with borosilicate composition,
3. equilibrate the model, and
4. run a **tensile stretch-to-failure** test.

> Important: a truly real-size fiber (e.g., 125 µm OD, 8–20 µm ID, 0.5–5 mm gauge length) contains an astronomically large number of atoms for fully atomistic ReaxFF. The provided workflow still builds the exact physical dimensions and estimates atom count, then supports either:
> - exporting a smaller tractable gauge section, or
> - submitting ultra-large runs on extreme HPC.

## Files

- `scripts/generate_hollow_fiber.py`: Builds a hollow cylindrical shell geometry and writes a LAMMPS `data` file.
- `inputs/in.reaxff_tensile_failure.lmp`: ReaxFF relaxation + tensile loading deck.
- `inputs/species.map`: Species mapping for ReaxFF analysis.
- `data/fiber_hollow.data`: Pre-generated initial hollow-fiber data file (50,000 atoms, isotropically scaled from real-size target).

## Typical use

```bash
python3 scripts/generate_hollow_fiber.py \
  --outer-diameter-um 125 \
  --inner-diameter-um 10 \
  --gauge-length-um 500 \
  --density-g-cc 2.23 \
  --max-atoms 1200000 \
  --output data/fiber_hollow.data
```

Then run:

```bash
lmp -in inputs/in.reaxff_tensile_failure.lmp \
    -var data_file data/fiber_hollow.data \
    -var ffield_file ffield.reax.borosilicate
```

## Notes

- The generator always computes and reports the real-size atom estimate.
- If the requested atom count exceeds `--max-atoms`, it automatically downscales all dimensions isotropically while preserving OD/ID ratio and composition.
- To force writing the full real-size model (usually impractical), pass `--allow-huge`.
