import argparse
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?"
RESIDUAL_ROW = re.compile(
    rf"^\s*(\d+)\s+\|\s+(.+?)\s+\d+\s+({NUMBER})\s+({NUMBER})\s*$"
)
ORIGINAL_RESIDUAL = re.compile(
    rf"^\s*(\d+)\s+\|\s+Projection:\s*([^,]+),.*Original res\.:\s*({NUMBER})"
)
FLUID_SECTION = re.compile(r"-{3,}\s*Fluid\s*-{3,}", re.IGNORECASE)

def read_residuals(log_path: Path) -> dict[str, dict[str, tuple[list[int], list[float]]]]:
    """Return solver residuals indexed by simulation and field."""
    values: dict[str, dict[str, list[tuple[int, float]]]] = {
        "Forward": defaultdict(list),
        "Adjoint": defaultdict(list),
    }

    current_simulation = "Forward"
    with log_path.open(encoding="utf-8", errors="replace") as log_file:
        for line in log_file:
            if "Adjoint fluid" in line:
                current_simulation = "Adjoint"
            elif FLUID_SECTION.search(line):
                current_simulation = "Forward"

            original_match = ORIGINAL_RESIDUAL.match(line)
            if original_match is not None:
                step, field, original_residual = original_match.groups()
                residual = float(
                    original_residual.replace("D", "E").replace("d", "e")
                )
                values[current_simulation]["Original res."].append((int(step), residual))
                continue

            match = RESIDUAL_ROW.match(line)
            if match is None:
                continue

            step, field, _start_residual, final_residual = match.groups()
            simulation = "Adjoint" if "adjoint" in field.lower() else current_simulation
            residual = float(final_residual.replace("D", "E").replace("d", "e"))
            values[simulation][field.strip()].append((int(step), residual))

    result = {}
    for simulation, fields in values.items():
        result[simulation] = {
            field: (
                [step for step, _ in samples],
                [residual for _, residual in samples],
            )
            for field, samples in fields.items()
        }
    return result


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("log", type=Path, help="path to the Neko log file")
args = parser.parse_args()

residuals = read_residuals(args.log)
for simulation, fields in residuals.items():
    if not fields:
        parser.error(f"no {simulation.lower()} residuals found in {args.log}")

figure, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
for axis, (simulation, fields) in zip(axes, residuals.items()):
    for field, (steps, values) in fields.items():
        axis.semilogy(steps, values, label=field)

    axis.set_title(f"{simulation} solver residuals")
    axis.set_xlabel("Time step")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend()

axes[0].set_ylabel("Residual")
figure.tight_layout()
plt.show()