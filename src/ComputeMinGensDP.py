## MINIMAL GENERATORS AND DEATH TIMES
# This script defines a function to compute minimal generators for H_1 of an Alpha filtration
# plus the death point of the persistent homology classes they belong to

import os
import numpy as np
import scipy as sp
import scipy.sparse
import networkx as nx
import pandas as pd
import scipy.spatial
import datetime

from plotting import draw_2d_simplicial_complex

import DriverN as DriverN
from acomplex import AlphaComplex

# import the cythonized low function from CyUtils 
from CyUtils import low as low
    
        
### DEFINE

def computeMinGensAndDeathPoints(Thresh, data_folder, output_folder, plot=False, picture_folder=None, 
                                 deathpoints=True, verbose=False, square=False):
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


    
    # Process each file in the data folder
    for filename in os.listdir(data_folder):

        if filename.endswith('.csv'):
            if verbose:
                print('Processing: ',filename, flush=True)

            file_path = os.path.join(data_folder, filename)        
            name, _ = os.path.splitext(filename)
            
            if verbose:
                print('File = '+name, flush=True)

            # if it's not already been processed
            #if not os.path.exists(os.path.join(picture_folder,name+'.png')):
            if verbose:
                print('Process', flush=True)

            # Read
            Points = pd.read_csv(file_path)
            Points = Points.to_numpy()[:,1:3] # Turn to numpy
            
            if plot:
                
                # generate limits to plot
                x_min = np.min(Points[:,0])
                x_max = np.max(Points[:,0])
                y_min = np.min(Points[:,1])
                y_max = np.max(Points[:,1])
                
                x_padding = (x_max - x_min)/10
                y_padding = (y_max - y_min)/10
                
                limits = [ [x_min - x_padding, x_max + x_padding] ,
                          [y_min - y_padding  , y_max + y_padding] ]

            DataSample = Points.copy()
            if verbose:
                print('Number of points: '+ str(DataSample.shape), flush=True)
                print('Build complex', flush=True)
            # Process
            
            # Legacy
            isItShort = False
            
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

            if verbose:
                print('Begin reading- ' + str(datetime.datetime.now()), flush=True)

            start = datetime.datetime.now()

            ## COMPUTE THE GENS
            Filtr, cycles, eL = DriverN.getFiltrBasisN(DataSample, epsList, method='Alpha', shortEdges = isItShort)


            end = datetime.datetime.now()
            
            if verbose:
                print('Duration - ' + str(end-start), flush=True)

            # Adjust output
         
            cycles = cycles.as_list()
            cycles = [ x[0].tolist() for x in cycles ]
            eL = np.array(eL)
            
            if verbose:
                if len(cycles) != 0:
                    print('The shape of cycles is ',str(len(cycles[0])) , " x ",str(len(cycles)))
                else:
                    print('Cycles is empty')
                    
                print('Shape of eL ', eL.shape)

                if isItShort:
                    print('Complex has ' + str(len(TriEdges)) + ' edges', flush=True)

                else:
                    print('Complex has ' + str(len(EL)) + ' edges', flush=True)

                print('Basis has ' + str(len(cycles)) + ' cycles', flush=True)

            
            if len(cycles) != 0:
                SHB = np.matrix( cycles ).T
            else:
                SHB = np.zeros( (len(eL),1) )


            Sum = np.sum(SHB, axis=1)
            mEdges = [ e for j,e in enumerate(eL) if Sum[j] > 0  ]


            SHB_filename = os.path.join(output_folder, name+'.npy')
            eL_filename = os.path.join(output_folder, name+'EdgeList.npy')

            np.save(SHB_filename, SHB)
            np.save(eL_filename, eL)

            #move a bit

            if plot:
                fig = draw_2d_simplicial_complex(simplices, limits, pos=DataSample, markedEdges=mEdges, fig=None)
                fig.savefig(os.path.join(picture_folder, name+'.png'))
                if verbose:
                    print('Saved pics', flush=True)

            
            if deathpoints:

                ## PERSISTENCE PAIRS 
                if verbose:
                    print('Computing death point of basis cycles', flush=True)
                    
                if len(cycles) == 0:
                    
                    if verbose:
                        print('No cycles to compute', flush=True)
                    return


                PersIntervals = np.zeros((len(cycles), 2))

                for i in range(len(cycles)):
                    PersIntervals[i,1] = np.infty


                AllTriList, FullEL , AllTriScales, _ = AlphaComplex(DataSample, np.infty)

                #print('How many scales? ', len(AllTriScales), flush=True)

                FullEL = [tuple(x.tolist()) for x in FullEL]
                FullEL = list(set(FullEL))
                FullEL = [list(x) for x in FullEL]

                FullEL = list(sorted(FullEL))
                
                #print('Len of the full EdgeList = ', len(FullEL), flush=True)

                MapEdges = [FullEL.index(t.tolist()) for t in eL]

                #print('MapEdges = ',len(MapEdges), flush=True)

                FullBMatrix = np.zeros( ( len(FullEL), len(AllTriList) ) , dtype=int)

                # Generate the full boundary matrix of triangles
                for i,t in enumerate(AllTriList):

                    T = sorted(t)

                    #this does not require mapping to the larger edgelist
                    # as it is already in the large edgelist
                    e0 = FullEL.index( [T[0], T[1]] )
                    e1 = FullEL.index( [T[1], T[2]] )
                    e2 = FullEL.index( [T[0], T[2]] )

                    FullBMatrix[ e0 , i ] = 1
                    FullBMatrix[ e1 , i ] = 1
                    FullBMatrix[ e2 , i ] = 1
                    
                if verbose:
                    print('Full boundary matrix computed ', flush=True)


                def sumMod2( v , w):

                    return np.mod(v+w,2)

                def reduceCol(ii,mat):

                    l = low(mat[:,ii])
                    j = 0

                    while True:

                        if j == ii:

                            break

                        if low( mat[:,j] ) == l and l != -1:

                            mat[:,ii] = sumMod2( mat[:,ii], mat[:,j]) 
                            l = low(mat[:,ii])
                            j = 0

                        else:
                            j += 1

                    return l == -1


                #Translate cycles in the larger edgelist
                Cycles = np.zeros((len(FullEL), len(cycles)), dtype=int)

                if verbose:
                    print('Translating cycles ', flush=True)
                    
                for i,c in enumerate(cycles):

                    # TRANSLATE CYCLES IN THE LARGER EDGELIST

                    for jj, entry in enumerate(c):

                        if entry > 0:

                            # Map to the full edgelist
                            LargerIndex = MapEdges[jj]
                            Cycles[LargerIndex,i] = 1
                            

                # Turn to list
                Cycles = [ Cycles[:,i] for i in range(Cycles.shape[1]) ]


                # The death time is the smallest scale at which sufficiently many triangles appear to fill a cycle
                # This will require several calls of a linear system

                # Initialize up until Thresh

                firstCols = len( [x for x in AllTriScales if x <= Thresh] )

                BMatrix = np.zeros((len(FullEL), firstCols), dtype=int)

                # Since it is slicing, it MAKES A COPY and there is no referencing problem, i.e. 
                #changes to BMatrix do not affect FullBMatrix
                
                BMatrix[:, :firstCols] = FullBMatrix[:, :firstCols]

                if verbose:
                    print('Initial boundary matrix computed ', flush=True)
                #print(BMatrix)

                # reduce the initial boundary matrix 
                if verbose:
                    print('Start reducing initial matrix (' + str(firstCols)+' columns) ' , flush=True)
                    
                for h in range(firstCols):

                    reduceCol(h, BMatrix)
                    if verbose and h % 20 == 0:
                        print('Done ', str(h), ' out of ', str(firstCols), flush=True)
                        
                if verbose:
                    print('Finished reducing initial matrix', flush=True)

                ## SAVE this matrix for later! We need it to compute birth times
                import copy
                ReducedTriBMatrix = copy.deepcopy(BMatrix)

                steps = [x for x in AllTriScales if x >= Thresh]

                LenSteps = len(steps)

                CyclesToKill = [x for x in range(len(Cycles))]

                for stepInd, s in enumerate(steps):
                    
                    if verbose:
                        print('Working step ', stepInd, ' of ', LenSteps, flush=True)

                    BMatrix = np.c_[BMatrix, FullBMatrix[:,firstCols+stepInd]] # add new column

                    reduceCol( firstCols+stepInd , BMatrix ) # reduce the last column added

                    # now for every cycle in Cycles (the ones in the LARGE edgelist!)
                    for i, cycle in enumerate(Cycles):

                        if i in CyclesToKill:

                            system = np.c_[BMatrix, cycle]

                            #print('System matrix is ')
                            #print(system)

                            lastCol = system.shape[1] - 1

                            #print('Last column is ', lastCol)

                            if reduceCol( lastCol , system ): # if it gets all zero it's linearly dependent

                                # FOUND the death step
                                if verbose:
                                    print('Death of cycle',i,' found at step ' +str(stepInd)+' and scale ' + str(s), flush=True)
                                PersIntervals[i,1] = s

                                CyclesToKill.remove(i)

                            if len(CyclesToKill) == 0:
                                break

                    if len(CyclesToKill) == 0:
                        if verbose:
                            print('Reached the end of cycles to kill')
                        break


                #print('Persistence Pairs: ', PersIntervals)
                
                DeathPoints = PersIntervals[:,1]
                
                if square:
                    DeathPoints = DeathPoints**2
                    
                    if verbose:
                        print('Squaring death filtration values')

                DP_filename = os.path.join(output_folder, name+'DP.npy')
                np.save(DP_filename, DeathPoints)

