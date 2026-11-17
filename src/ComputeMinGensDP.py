## MINIMAL GENERATORS AND DEATH TIMES
# This script defines a function to compute minimal generators for H_1 of an Alpha filtration
# plus the death point of the persistent homology classes they belong to
#
# VERSION 2026-11 (parallel slices + incremental death points). Same
# outputs as the original (kept as Old_ComputeMinGensDP.py).
#
# Parallel slices:
#   The work on one csv file is the module-level function processFile;
#   computeMinGensAndDeathPoints runs it on every file, one after the other
#   (n_jobs=1, default) or in parallel over slices (n_jobs > 1). Slices are
#   independent, so the outputs do not depend on n_jobs. A slice with no
#   cycles saves an empty death-point file and moves on; figures are closed
#   after saving; files are processed in sorted order.
#
# Incremental death points:
#   The death of a minimal generator c is the first triangle (in filtration
#   order, after theta) at which c lies in the span of the boundaries of
#   all triangles so far. The original answered this at every step and for
#   every alive cycle by copying the dense boundary matrix and reducing c
#   from scratch. We instead keep:
#     - pivots: lowest edge -> reduced boundary column, for the columns
#       added so far (standard left-to-right reduction over Z2);
#     - for each cycle, its residual: the cycle reduced against the pivots.
#   A residual whose lowest edge is not a pivot is not a combination of the
#   reduced columns, so the cycle is not yet a boundary. When a new triangle
#   adds a pivot p, only residuals with lowest edge p can change; we reduce
#   those. A cycle dies at the first step where its residual becomes zero,
#   which is the same test as the original, at the same step, so the death
#   value is the same entry of AllTriScales.
#   Columns are stored as Python ints used as bit vectors over the edges:
#   bit k is edge k of FullEL; adding two columns mod 2 is XOR; the lowest
#   (largest-index) edge of a column is bit_length() - 1.

import os
import numpy as np
import scipy as sp
import scipy.sparse
import networkx as nx
import pandas as pd
import scipy.spatial
import datetime
import concurrent.futures

import matplotlib.pyplot as plt

from Plotting import draw_2d_simplicial_complex

import DriverN as DriverN
from acomplex import AlphaComplex


def computeDeathPoints(DataSample, cycles, eL, Thresh, verbose=False):
    '''Death scale of each minimal generator, by incremental reduction.

    INPUT:
    DataSample (N x 2 array) the point cloud
    cycles (list of 0-1 lists) minimal generators, one list per cycle,
        indexed by the edges of eL
    eL (#edges x 2 array) edge list of the complex at scale Thresh
    Thresh (float) scale theta (radius)

    OUTPUT:
    (#cycles,) array of death scales (radius); np.inf if a cycle never dies.
    '''
    nCycles = len(cycles)
    deaths = np.full(nCycles, np.inf)

    # Full Alpha filtration: triangles in filtration order, with scales
    AllTriList, FullEL, AllTriScales, _ = AlphaComplex(DataSample, np.inf)

    # Full edge list, unique and sorted lexicographically (as the original)
    FullEL = [tuple(x.tolist()) for x in FullEL]
    FullEL = list(set(FullEL))
    FullEL = sorted(FullEL)

    # Edge -> row index, instead of FullEL.index (a linear search)
    edgeRow = {}
    for k, e in enumerate(FullEL):
        edgeRow[(int(e[0]), int(e[1]))] = k

    # Boundary of each triangle as a bit vector over the full edge list
    triCols = []
    for t in AllTriList:
        T = sorted(int(v) for v in t)
        e0 = edgeRow[(T[0], T[1])]
        e1 = edgeRow[(T[1], T[2])]
        e2 = edgeRow[(T[0], T[2])]
        triCols.append((1 << e0) ^ (1 << e1) ^ (1 << e2))

    if verbose:
        print('Full boundary columns computed ', flush=True)

    # Each cycle as a bit vector over the full edge list
    mapEdges = [edgeRow[(int(t[0]), int(t[1]))] for t in eL]
    cycleCols = []
    for c in cycles:
        col = 0
        for jj, entry in enumerate(c):
            if entry > 0:
                col ^= (1 << mapEdges[jj])
        cycleCols.append(col)

    # Pivots of the reduced boundary columns: lowest edge -> column
    pivots = {}

    def reduce(col):
        '''Reduce a column against the current pivots; return the residual.'''
        while col:
            lowEdge = col.bit_length() - 1
            pivotCol = pivots.get(lowEdge)
            if pivotCol is None:
                return col
            col ^= pivotCol
        return 0

    # Reduce the boundary columns of the triangles with scale <= theta
    firstCols = len([x for x in AllTriScales if x <= Thresh])
    for h in range(firstCols):
        col = reduce(triCols[h])
        if col:
            pivots[col.bit_length() - 1] = col

    if verbose:
        print('Reduced initial boundary (' + str(firstCols) + ' columns)', flush=True)

    # Residual of each cycle against the boundaries up to theta
    residuals = [reduce(c) for c in cycleCols]

    # Walk the triangles after theta, as the original: step stepInd adds
    # column firstCols + stepInd and reports scale steps[stepInd]
    steps = [x for x in AllTriScales if x >= Thresh]
    alive = list(range(nCycles))

    for stepInd, s in enumerate(steps):

        colIdx = firstCols + stepInd
        if colIdx >= len(triCols):
            # The original would raise IndexError here; this can only happen
            # if a scale equals theta exactly and some cycle never dies.
            break

        # Reduce the new triangle; a nonzero result adds a pivot
        newCol = reduce(triCols[colIdx])
        newLow = -1
        if newCol:
            newLow = newCol.bit_length() - 1
            pivots[newLow] = newCol

        # Only residuals whose lowest edge is the new pivot can change
        stillAlive = []
        for i in alive:
            r = residuals[i]
            if r and newLow >= 0 and (r.bit_length() - 1) == newLow:
                r = reduce(r)
                residuals[i] = r
            if r == 0:
                deaths[i] = s
                if verbose:
                    print('Death of cycle', i, ' found at step ' + str(stepInd) + ' and scale ' + str(s), flush=True)
            else:
                stillAlive.append(i)
        alive = stillAlive

        if len(alive) == 0:
            if verbose:
                print('Reached the end of cycles to kill')
            break

    return deaths


def processFile(filename, Thresh, data_folder, output_folder, plot, picture_folder,
                deathpoints, verbose, square):
    '''Compute minimal generators and death points of one csv file.

    Writes, in output_folder:
        <name>.npy          minimal generators, (#edges x #cycles) 0-1 matrix
        <name>EdgeList.npy  edge list, (#edges x 2) point indices
        <name>DP.npy        death scale of each generator, (#cycles,)
    and, if plot, two figures in picture_folder.

    Returns (filename, seconds) for the timing summary.
    '''
    t_start = datetime.datetime.now()

    if verbose:
        print('Processing: ', filename, flush=True)

    file_path = os.path.join(data_folder, filename)
    name, _ = os.path.splitext(filename)

    # Read
    Points = pd.read_csv(file_path)
    Points = Points.to_numpy()[:, 1:3]  # Turn to numpy

    if plot:

        # generate limits to plot
        x_min = np.min(Points[:, 0])
        x_max = np.max(Points[:, 0])
        y_min = np.min(Points[:, 1])
        y_max = np.max(Points[:, 1])

        x_padding = (x_max - x_min)/10
        y_padding = (y_max - y_min)/10

        limits = [[x_min - x_padding, x_max + x_padding],
                  [y_min - y_padding, y_max + y_padding]]

    DataSample = Points.copy()
    if verbose:
        print('Number of points: ' + str(DataSample.shape), flush=True)

    # generate complex
    TriList, EL, TriScales, EdgeScales = AlphaComplex(DataSample, Thresh)

    # adapt lists of edges and simplices
    EL = [x.tolist() for x in EL]
    simplices = EL + TriList

    # Generate a list of Thresh parameters to use
    epsList = [Thresh]

    if plot:
        fig = draw_2d_simplicial_complex(simplices, limits, pos=DataSample, markedEdges=None, fig=None)
        fig.savefig(os.path.join(picture_folder, name+'Complex.png'))
        plt.close(fig)

    start = datetime.datetime.now()

    ## COMPUTE THE GENS
    Filtr, cycles, eL = DriverN.getFiltrBasisN(DataSample, epsList, method='Alpha', shortEdges=False)

    end = datetime.datetime.now()

    if verbose:
        print('Duration - ' + str(end-start), flush=True)

    # Adjust output
    cycles = cycles.as_list()
    cycles = [x[0].tolist() for x in cycles]
    eL = np.array(eL)

    if verbose:
        print('Complex has ' + str(len(EL)) + ' edges', flush=True)
        print('Basis has ' + str(len(cycles)) + ' cycles', flush=True)

    if len(cycles) != 0:
        SHB = np.matrix(cycles).T
    else:
        SHB = np.zeros((len(eL), 1))

    Sum = np.sum(SHB, axis=1)
    mEdges = [e for j, e in enumerate(eL) if Sum[j] > 0]

    SHB_filename = os.path.join(output_folder, name+'.npy')
    eL_filename = os.path.join(output_folder, name+'EdgeList.npy')

    np.save(SHB_filename, SHB)
    np.save(eL_filename, eL)

    if plot:
        fig = draw_2d_simplicial_complex(simplices, limits, pos=DataSample, markedEdges=mEdges, fig=None)
        fig.savefig(os.path.join(picture_folder, name+'.png'))
        plt.close(fig)

    if deathpoints:

        if verbose:
            print('Computing death point of basis cycles', flush=True)

        if len(cycles) == 0:
            # We save an empty death-point file so every slice has one
            np.save(os.path.join(output_folder, name+'DP.npy'), np.zeros(0))
            seconds = (datetime.datetime.now() - t_start).total_seconds()
            return (filename, seconds)

        DeathPoints = computeDeathPoints(DataSample, cycles, eL, Thresh, verbose=verbose)

        if square:
            DeathPoints = DeathPoints**2

        DP_filename = os.path.join(output_folder, name+'DP.npy')
        np.save(DP_filename, DeathPoints)

    seconds = (datetime.datetime.now() - t_start).total_seconds()
    return (filename, seconds)


def _processFileTask(task):
    '''Unpack one task tuple for the process pool (the pool passes a single argument).'''
    return processFile(*task)


def computeMinGensAndDeathPoints(Thresh, data_folder, output_folder, plot=False, picture_folder=None,
                                 deathpoints=True, verbose=False, square=False, n_jobs=1):
    '''Function to compute minimal generators and death points

    Analyzes all csv files in data_folder, for each generates the Alpha complex at scale Thresh, then
    computes the minimal generators for H_1 of that complex. If deathpoints, computes the Alpha filtration
    scale of the death of the persistent homology class of each generator.
    Outputs to the output folder, based on filenames.

    INPUT:
    Thresh (float) scale of the Alpha complex
    data_folder (path) path to the input image (csv file)
    output_folder (path) path to the folder to store the output
    plot (bool) whether to produce images
    picture_folder (path) where to store images
    deathpoints (bool) whether to compute the death scale of the classes
    verbose (bool) print as it executes
    square (bool) if True, return the SQUARED death scale of cycles. Otherwise, the regular one.
    n_jobs (int) number of slices processed in parallel; 1 = one after the other.
        Outputs are the same for any n_jobs.
    '''

    # Check types are correct

    if not isinstance(Thresh, (float, int)):
        raise TypeError(f"Thresh must be a float, but got {type(Thresh).__name__}")
    if not isinstance(data_folder, (str, bytes, os.PathLike)):
        raise TypeError(f"data_folder must be a valid path (str, bytes, or os.PathLike), but got {type(data_folder).__name__}")
    if not isinstance(output_folder, (str, bytes, os.PathLike)):
        raise TypeError(f"output_folder must be a valid path (str, bytes, or os.PathLike), but got {type(output_folder).__name__}")
    if not isinstance(plot, bool):
        raise TypeError(f"plot must be a bool, but got {type(plot).__name__}")
    if picture_folder is not None and not isinstance(picture_folder, (str, bytes, os.PathLike)):
        raise TypeError(f"picture_folder must be a valid path (str, bytes, or os.PathLike), but got {type(picture_folder).__name__}")
    if not isinstance(deathpoints, bool):
        raise TypeError(f"deathpoints must be a bool, but got {type(deathpoints).__name__}")
    if not isinstance(verbose, bool):
        raise TypeError(f"verbose must be a bool, but got {type(verbose).__name__}")
    if not isinstance(square, bool):
        raise TypeError(f"square must be a bool, but got {type(square).__name__}")
    if not isinstance(n_jobs, int) or n_jobs < 1:
        raise ValueError(f"n_jobs must be a positive int, but got {n_jobs}")

    # Check valid inputs
    if Thresh <= 0:
        raise ValueError("Thresh must be a positive number")
    if plot and picture_folder is None:
        raise ValueError("Must specify a location to save pictures if plot is True")
    if not os.path.exists(data_folder):
        raise FileNotFoundError(f"The data_folder path '{data_folder}' does not exist")
    if not os.path.isdir(data_folder):
        raise NotADirectoryError(f"The data_folder path '{data_folder}' is not a directory")
    if not os.path.exists(output_folder):
        raise FileNotFoundError(f"The output_folder path '{output_folder}' does not exist")
    if not os.path.isdir(output_folder):
        raise NotADirectoryError(f"The output_folder path '{output_folder}' is not a directory")
    if picture_folder is not None and not os.path.exists(picture_folder):
        raise FileNotFoundError(f"The picture_folder path '{picture_folder}' does not exist")
    if picture_folder is not None and not os.path.isdir(picture_folder):
        raise NotADirectoryError(f"The picture_folder path '{picture_folder}' is not a directory")

    # One task per csv file, in sorted order
    csv_files = sorted(f for f in os.listdir(data_folder) if f.endswith('.csv'))
    tasks = [(f, Thresh, data_folder, output_folder, plot, picture_folder, deathpoints, verbose, square)
             for f in csv_files]

    t_total = datetime.datetime.now()
    timings = []

    if n_jobs == 1:
        for task in tasks:
            timings.append(processFile(*task))
    else:
        # Each worker handles whole slices; an error in any slice is raised
        # here by f.result(), stopping the run.
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = [executor.submit(_processFileTask, task) for task in tasks]
            for f in concurrent.futures.as_completed(futures):
                timings.append(f.result())

    total_seconds = (datetime.datetime.now() - t_total).total_seconds()
    print(f'MinGens: {len(tasks)} slices in {total_seconds:.1f} s with n_jobs={n_jobs}', flush=True)
    for fname, secs in sorted(timings):
        print(f'  {fname}: {secs:.1f} s', flush=True)
