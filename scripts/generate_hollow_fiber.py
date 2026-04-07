#!/usr/bin/env python3
"""Generate a hollow borosilicate optical fiber shell for LAMMPS/ReaxFF."""

from __future__ import annotations

import argparse
import math
import random
from dataclasses import dataclass
from pathlib import Path

AVOGADRO = 6.02214076e23
ANGSTROM_PER_UM = 1.0e4


@dataclass
class FiberSpec:
    outer_diameter_um: float
    inner_diameter_um: float
    gauge_length_um: float
    density_g_cc: float


# Approximate borosilicate atom fractions (sum to 1.0).
ATOM_FRACTIONS = {
    "Si": 0.208,
    "B": 0.052,
    "O": 0.664,
    "Na": 0.040,
    "Al": 0.036,
}

ATOM_TYPES = {"Si": 1, "B": 2, "O": 3, "Na": 4, "Al": 5}

ATOMIC_MASS = {
    "Si": 28.085,
    "B": 10.81,
    "O": 15.999,
    "Na": 22.989769,
    "Al": 26.981538,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outer-diameter-um", type=float, default=125.0)
    p.add_argument("--inner-diameter-um", type=float, default=10.0)
    p.add_argument("--gauge-length-um", type=float, default=500.0)
    p.add_argument("--density-g-cc", type=float, default=2.23)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--max-atoms", type=int, default=1_200_000)
    p.add_argument("--allow-huge", action="store_true")
    p.add_argument("--output", type=Path, default=Path("data/fiber_hollow.data"))
    return p.parse_args()


def average_atomic_mass() -> float:
    return sum(ATOM_FRACTIONS[s] * ATOMIC_MASS[s] for s in ATOM_FRACTIONS)


def shell_volume_cc(spec: FiberSpec) -> float:
    ro_cm = (spec.outer_diameter_um / 2.0) * 1e-4
    ri_cm = (spec.inner_diameter_um / 2.0) * 1e-4
    length_cm = spec.gauge_length_um * 1e-4
    return math.pi * (ro_cm**2 - ri_cm**2) * length_cm


def estimate_atom_count(spec: FiberSpec) -> int:
    v_cc = shell_volume_cc(spec)
    mass_g = v_cc * spec.density_g_cc
    m_avg = average_atomic_mass()  # g/mol per atom-average basis
    atoms = (mass_g / m_avg) * AVOGADRO
    return int(atoms)


def isotropic_scaled_spec(spec: FiberSpec, target_atoms: int) -> FiberSpec:
    est_atoms = estimate_atom_count(spec)
    if est_atoms <= target_atoms:
        return spec
    scale = (target_atoms / est_atoms) ** (1.0 / 3.0)
    return FiberSpec(
        outer_diameter_um=spec.outer_diameter_um * scale,
        inner_diameter_um=spec.inner_diameter_um * scale,
        gauge_length_um=spec.gauge_length_um * scale,
        density_g_cc=spec.density_g_cc,
    )


def choose_species(n: int, rng: random.Random) -> list[str]:
    species = list(ATOM_FRACTIONS.keys())
    weights = [ATOM_FRACTIONS[s] for s in species]
    return rng.choices(species, weights=weights, k=n)


def generate_atoms(spec: FiberSpec, natoms: int, seed: int):
    rng = random.Random(seed)
    ro = (spec.outer_diameter_um / 2.0) * ANGSTROM_PER_UM
    ri = (spec.inner_diameter_um / 2.0) * ANGSTROM_PER_UM
    zmax = spec.gauge_length_um * ANGSTROM_PER_UM

    chosen = choose_species(natoms, rng)
    atoms = []
    for idx, sym in enumerate(chosen, start=1):
        # Uniform sampling in annulus via sqrt transform.
        r2 = rng.uniform(ri * ri, ro * ro)
        r = math.sqrt(r2)
        theta = rng.uniform(0.0, 2.0 * math.pi)
        x = r * math.cos(theta)
        y = r * math.sin(theta)
        z = rng.uniform(0.0, zmax)
        charge = 0.0
        atoms.append((idx, ATOM_TYPES[sym], charge, x, y, z))
    return atoms, ro, zmax


def write_lammps_data(path: Path, spec: FiberSpec, natoms: int, atoms, ro: float, zmax: float):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("# Hollow borosilicate optical fiber shell (ReaxFF-ready)\n\n")
        f.write(f"{natoms} atoms\n")
        f.write("5 atom types\n\n")
        f.write(f"{-ro:.6f} {ro:.6f} xlo xhi\n")
        f.write(f"{-ro:.6f} {ro:.6f} ylo yhi\n")
        f.write(f"0.000000 {zmax:.6f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for sym, t in ATOM_TYPES.items():
            f.write(f"{t} {ATOMIC_MASS[sym]:.6f} # {sym}\n")
        f.write("\nAtoms # charge\n\n")
        for a in atoms:
            idx, atype, q, x, y, z = a
            f.write(f"{idx} {atype} {q:.6f} {x:.6f} {y:.6f} {z:.6f}\n")


def main() -> None:
    args = parse_args()
    if args.inner_diameter_um >= args.outer_diameter_um:
        raise SystemExit("inner diameter must be smaller than outer diameter")

    real = FiberSpec(
        outer_diameter_um=args.outer_diameter_um,
        inner_diameter_um=args.inner_diameter_um,
        gauge_length_um=args.gauge_length_um,
        density_g_cc=args.density_g_cc,
    )
    est_real_atoms = estimate_atom_count(real)

    model = real
    natoms = est_real_atoms
    scaled = False

    if est_real_atoms > args.max_atoms and not args.allow_huge:
        model = isotropic_scaled_spec(real, args.max_atoms)
        natoms = estimate_atom_count(model)
        scaled = True

    if natoms <= 0:
        raise SystemExit("calculated atom count is zero; increase geometry size or density")

    atoms, ro, zmax = generate_atoms(model, natoms, args.seed)
    write_lammps_data(args.output, model, natoms, atoms, ro, zmax)

    print("=== Hollow borosilicate fiber generation summary ===")
    print(f"Real-size geometry: OD={real.outer_diameter_um} um, ID={real.inner_diameter_um} um, L={real.gauge_length_um} um")
    print(f"Estimated atoms for real-size geometry: {est_real_atoms:,}")
    if scaled:
        print("Scaled dimensions for tractable run (isotropic):")
        print(f"  OD={model.outer_diameter_um:.6e} um, ID={model.inner_diameter_um:.6e} um, L={model.gauge_length_um:.6e} um")
    print(f"Written atoms: {natoms:,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
