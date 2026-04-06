#!/usr/bin/env python3
"""
Generate an initial borosilicate glass configuration for LAMMPS.

This creates a random-packed starting point (not yet equilibrated) and writes
`borosilicate_init.data`.

Recommended workflow:
  1) python3 generate_borosilicate_glass.py
  2) lmp -in in.equilibrate_borosilicate_reaxff
  3) lmp -in in.ion_exchange_reaxff
"""

from __future__ import annotations

import argparse
import math
import random
from dataclasses import dataclass


@dataclass(frozen=True)
class Species:
    type_id: int
    label: str
    mass: float


SPECIES = [
    Species(1, "O", 15.999),
    Species(2, "Si", 28.085),
    Species(3, "B", 10.81),
    Species(4, "Na", 22.989769),
    Species(5, "K", 39.0983),
]


def allocate_counts(total_atoms: int, composition: dict[str, float]) -> dict[str, int]:
    # Largest-remainder allocation that sums exactly to total_atoms.
    raw = {k: total_atoms * v for k, v in composition.items()}
    counts = {k: int(math.floor(v)) for k, v in raw.items()}
    remainder = total_atoms - sum(counts.values())
    order = sorted(raw.keys(), key=lambda k: raw[k] - counts[k], reverse=True)
    for i in range(remainder):
        counts[order[i % len(order)]] += 1
    return counts


def random_positions(n: int, box: float, min_sep: float, seed: int) -> list[tuple[float, float, float]]:
    rng = random.Random(seed)
    pos: list[tuple[float, float, float]] = []
    max_trials = 30000
    for _ in range(n):
        accepted = False
        for _trial in range(max_trials):
            x = rng.random() * box
            y = rng.random() * box
            z = rng.random() * box
            good = True
            for px, py, pz in pos:
                dx = x - px
                dy = y - py
                dz = z - pz
                if dx * dx + dy * dy + dz * dz < min_sep * min_sep:
                    good = False
                    break
            if good:
                pos.append((x, y, z))
                accepted = True
                break
        if not accepted:
            raise RuntimeError(
                "Failed to place atoms without overlap. "
                "Increase box length or lower min separation."
            )
    return pos


def write_lammps_data(
    out_file: str,
    box: float,
    atom_types: list[int],
    positions: list[tuple[float, float, float]],
) -> None:
    masses = {s.type_id: s.mass for s in SPECIES}
    with open(out_file, "w", encoding="utf-8") as f:
        f.write("# Borosilicate glass initial structure (random packed)\n\n")
        f.write(f"{len(atom_types)} atoms\n")
        f.write(f"{len(SPECIES)} atom types\n\n")
        f.write(f"0.0 {box:.6f} xlo xhi\n")
        f.write(f"0.0 {box:.6f} ylo yhi\n")
        f.write(f"0.0 {box:.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        for tid in sorted(masses):
            f.write(f"{tid} {masses[tid]:.8f}\n")
        f.write("\nAtoms # charge\n\n")

        for i, (atype, (x, y, z)) in enumerate(zip(atom_types, positions), start=1):
            f.write(f"{i} {atype} 0.0 {x:.6f} {y:.6f} {z:.6f}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate borosilicate initial data file")
    parser.add_argument("--atoms", type=int, default=3000, help="Total atom count")
    parser.add_argument("--box", type=float, default=42.0, help="Cubic box length [Angstrom]")
    parser.add_argument("--min-sep", type=float, default=1.35, help="Minimum atom separation [Angstrom]")
    parser.add_argument("--seed", type=int, default=20260406, help="RNG seed")
    parser.add_argument("--out", default="borosilicate_init.data", help="Output LAMMPS data filename")
    args = parser.parse_args()

    # Example sodium borosilicate composition (atomic fractions).
    # Adjust as needed for your target glass chemistry.
    composition = {
        "O": 0.620,
        "Si": 0.220,
        "B": 0.110,
        "Na": 0.050,
        # K in initial glass set to zero. K enters during ion-exchange step.
        "K": 0.000,
    }

    if abs(sum(composition.values()) - 1.0) > 1e-9:
        raise ValueError("Composition fractions must sum to 1.0")

    counts = allocate_counts(args.atoms, composition)
    id_by_label = {s.label: s.type_id for s in SPECIES}

    atom_types: list[int] = []
    for label in ["O", "Si", "B", "Na", "K"]:
        atom_types.extend([id_by_label[label]] * counts[label])

    rng = random.Random(args.seed)
    rng.shuffle(atom_types)

    positions = random_positions(len(atom_types), args.box, args.min_sep, args.seed + 17)
    write_lammps_data(args.out, args.box, atom_types, positions)

    print(f"Wrote {args.out}")
    print("Atom counts:")
    for label in ["O", "Si", "B", "Na", "K"]:
        print(f"  {label}: {counts[label]}")


if __name__ == "__main__":
    main()
