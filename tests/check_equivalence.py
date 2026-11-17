"""
Check that two versions of MinGens produce identical outputs.

Usage, from the MinGens folder, with both versions compiled (make):

    python tests/check_equivalence.py <reference_dir> <candidate_dir> <data_dir> "<theta list>"

Example, comparing the old main (checked out as a worktree) with this one,
on the toy data at three scales:

    git worktree add ~/MinGens_old_draft old_draft
    (cd ~/MinGens_old_draft && make)
    python tests/check_equivalence.py ~/MinGens_old_draft . data "0.1 0.2 0.3"

<reference_dir> and <candidate_dir> are MinGens folders (each with its own
compiled src/). <data_dir> holds the csv point clouds (pandas index column
first, then x, y). Theta is a radius.

Each version runs in its own Python process, because both define modules
with the same names (CySPT, DriverN, ComputeMinGensDP) that cannot be
loaded twice in one process. For every theta we run
computeMinGensAndDeathPoints on all csv files of <data_dir>, with plots
off, and require every output file (generators, edge list, death points)
to be present in both versions and identical, value for value.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

# Code run in a separate process for one version: we put that version's
# src/ first on sys.path, so its own modules are the ones imported.
RUNNER = r'''
import sys
src_dir, thresh, data_dir, out_dir = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
sys.path.insert(0, src_dir)
import CySPT
from ComputeMinGensDP import computeMinGensAndDeathPoints
print("CySPT implementation:", getattr(CySPT, "CYSPT_IMPLEMENTATION", "original (no marker)"))
computeMinGensAndDeathPoints(Thresh=thresh, data_folder=data_dir, output_folder=out_dir,
                             plot=False, picture_folder=None, deathpoints=True,
                             verbose=False, square=False)
'''


def run_version(mingens_dir: str, thresh: float, data_dir: str, out_dir: str) -> float:
    """
    Run one MinGens version in a subprocess; return the wall-clock seconds.

    Parameters
    ----------
    mingens_dir : str
        MinGens folder of the version to run (its src/ must be compiled).
    thresh : float
        Scale theta (radius).
    data_dir : str
        Folder with the csv point clouds.
    out_dir : str
        Output folder, created empty.
    """
    os.makedirs(out_dir)
    src_dir = os.path.join(os.path.abspath(mingens_dir), 'src')
    t0 = time.perf_counter()
    completed = subprocess.run(
        [sys.executable, '-c', RUNNER, src_dir, str(thresh), os.path.abspath(data_dir), out_dir],
        capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if completed.returncode != 0:
        print(completed.stdout)
        print(completed.stderr)
        raise RuntimeError(f'MinGens run failed in {mingens_dir}')
    first_line = completed.stdout.splitlines()[0] if completed.stdout else ''
    print(f'    {mingens_dir}: {first_line}, {elapsed:.2f} s')
    return elapsed


def compare_folders(ref_out: str, new_out: str) -> list:
    """
    Compare all .npy outputs of two runs.

    Returns a list of problems (empty if every file exists in both runs and
    holds identical arrays).
    """
    problems = []
    ref_files = sorted(f for f in os.listdir(ref_out) if f.endswith('.npy'))
    new_files = sorted(f for f in os.listdir(new_out) if f.endswith('.npy'))

    # Files written by only one version
    only_ref = sorted(set(ref_files) - set(new_files))
    only_new = sorted(set(new_files) - set(ref_files))
    for f in only_ref:
        problems.append(f'{f}: only in reference')
    for f in only_new:
        problems.append(f'{f}: only in candidate')

    # Files written by both: shapes and values must match exactly
    for f in sorted(set(ref_files) & set(new_files)):
        ref = np.asarray(np.load(os.path.join(ref_out, f), allow_pickle=True))
        new = np.asarray(np.load(os.path.join(new_out, f), allow_pickle=True))
        if ref.shape != new.shape:
            problems.append(f'{f}: shape {ref.shape} vs {new.shape}')
        elif not np.array_equal(ref, new):
            problems.append(f'{f}: values differ')
    return problems


def main() -> None:
    ref_dir, new_dir, data_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    thetas = [float(x) for x in sys.argv[4].split()]

    scratch = tempfile.mkdtemp(prefix='mingens_equiv_')
    all_ok = True
    try:
        for thresh in thetas:
            print(f'theta = {thresh}')
            ref_out = os.path.join(scratch, f'ref_{thresh}')
            new_out = os.path.join(scratch, f'new_{thresh}')
            run_version(ref_dir, thresh, data_dir, ref_out)
            run_version(new_dir, thresh, data_dir, new_out)

            problems = compare_folders(ref_out, new_out)
            n_files = len([f for f in os.listdir(ref_out) if f.endswith('.npy')])
            if problems:
                all_ok = False
                for p in problems:
                    print(f'    DIFFERENCE {p}')
            else:
                print(f'    {n_files} output files identical')
    finally:
        shutil.rmtree(scratch)

    print('\nALL IDENTICAL' if all_ok else '\nDIFFERENCES FOUND')
    sys.exit(0 if all_ok else 1)


if __name__ == '__main__':
    main()
