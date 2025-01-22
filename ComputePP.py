## Computing persistence pairs of the minimal representatives found

import math
import numpy as np
import scipy as sp
import scipy.sparse
from random import random
import networkx as nx
import pickle as pk
import pandas as pd
import scipy.spatial
import datetime

import sys
import os 
os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"
import matplotlib.pyplot as plt

from acomplex import AlphaComplex

import DriverN
import GeometryN
import ScaffoldN


# import the cythonized low function from CyUtils 
from CyUtils import low as low

## define stuff

# def getEL(A):
#     edgesList = []
#     for i in range(len(A)):
#         for j in range(i,len(A)):
#             if A[i,j] > 0:
#                 edgesList.append([i,j])
#     return edgesList

from scipy.spatial import Delaunay
import scipy as sp

# def AlphaComplex(points, thresh, short):
#     ''' Takes a point cloud (Npoints x Ndims) and returns the Alpha complex at threshold thresh '''
    
#     NPoints = points.shape[0]
    
#     #build the circumsphere
#     def circumsphereRad(P):
#         '''
#         TAKEN FROM YOHAI REHANI
#         :param P: array Nxd of N points of dimension d
#         :return: center and squared radius of the circum-sphere of P
#         '''
#         p1 = P[0, :]
#         A = P[1:, :] - p1
#         Pnorm = np.sum(np.power(P, 2), 1)
#         b = 0.5*(Pnorm[1:] - Pnorm[0])
#         invA = np.linalg.pinv(A)
#         c0 = invA.dot(b)
#         F = sp.linalg.null_space(A)
#         if F.size != 0:
#             z = np.transpose(F).dot(p1-c0)
#             c = c0 + F.dot(z)
#         else:
#             c = c0
#         R = np.sum(np.power(p1-c, 2))
#         return R

#     Del = Delaunay(points)
#     #triangles = points[Del.simplices]
    
#     TriList = [ t for t in Del.simplices if circumsphereRad(points[t]) <= thresh**2 ]
#     TriList = [ np.array(sorted(t)) for t in TriList ]
    
#     # THIS IS THE FULL ALPHA COMPLEX
#     EdgeList = []
#     CandidateEdges = []
    
#     for t in Del.simplices:
            
#         e1 = tuple(sorted((t[0] ,t[1])))
#         e2 = tuple(sorted((t[0] ,t[2])))
#         e3 = tuple(sorted((t[1] ,t[2])))  
        
#         CandidateEdges.extend( [e1,e2,e3] )
    
#     CandidateEdges = list(set(CandidateEdges))
    
#     for e in CandidateEdges:
        
#         if circumsphereRad(points[[e[0],e[1]],:]) <= thresh**2 :
            
#             EdgeList.append(  (e[0],e[1]) )
            
#     EdgeList = [ np.array(e) for e in EdgeList ]
    
#     # THIS GIVES THE SHORT ALPHA COMPLEX
#     TriEdges = []
    
#     for t in TriList:
#         TriEdges.append( (t[0] ,t[1]) )
#         TriEdges.append( (t[0] ,t[2]) )
#         TriEdges.append( (t[1] ,t[2]) )
        
#     TriEdges = list(set(TriEdges))
#     TriEdges = [ np.array(e) for e in TriEdges ]
        
#     if not short:
#         return TriList, EdgeList  # This gives the true Alpha complex
#     else:
#         return TriList, TriEdges  # This gives the short Alpha complex
    
import itertools

def draw_2d_simplicial_complex(simplices, pos=None, return_pos=False, fig = None, markedEdges=None):
    """
    Draw a simplicial complex up to dimension 2 from a list of simplices, as in [1].
        
        Args
        ----
        simplices: list of lists of integerss
            List of simplices to draw. Sub-simplices are not needed (only maximal).
            For example, the 2-simplex [1,2,3] will automatically generate the three
            1-simplices [1,2],[2,3],[1,3] and the three 0-simplices [1],[2],[3].
            When a higher order simplex is entered only its sub-simplices
            up to D=2 will be drawn.
        
        pos: dict (default=None)
            If passed, this dictionary of positions d:(x,y) is used for placing the 0-simplices.
            The standard nx spring layour is used otherwise.
           
        ax: matplotlib.pyplot.axes (default=None)
        
        return_pos: dict (default=False)
            If True returns the dictionary of positions for the 0-simplices.
            
        References
        ----------    
        .. [1] I. Iacopini, G. Petri, A. Barrat & V. Latora (2018)
               "Simplicial Models of Social Contagion".
               arXiv preprint arXiv:1810.07031..
    """

    
    #List of 0-simplices
    nodes =list(set(itertools.chain(*simplices)))
    
    #List of 1-simplices
    edges = list(set(itertools.chain(*[[tuple(sorted((i, j))) for i, j in itertools.combinations(simplex, 2)] for simplex in simplices])))

    #List of 2-simplices
    triangles = list(set(itertools.chain(*[[tuple(sorted((i, j, k))) for i, j, k in itertools.combinations(simplex, 3)] for simplex in simplices])))
    
    fig, ax = plt.subplots(figsize=(7,7))
    ax.set_xlim([0, 270])      
    ax.set_ylim([0, 270])
    # ax.set_xlim([-3, 3])      
    # ax.set_ylim([-3, 3])
    ax.get_xaxis().set_ticks([])  
    ax.get_yaxis().set_ticks([])
    ax.axis('on')
       
    if pos is None:
        # Creating a networkx Graph from the edgelist
        G = nx.Graph()
        G.add_edges_from(edges)
        # Creating a dictionary for the position of the nodes
        pos = nx.spring_layout(G)
        
    # Drawing the edges
    for i, j in edges:
        (x0, y0) = pos[i]
        (x1, y1) = pos[j]
        line = plt.Line2D([ x0, x1 ], [y0, y1 ],color = 'blue', zorder = 1, lw=1.5)
        ax.add_line(line);
    
    # Filling in the triangles
    for i, j, k in triangles:
        (x0, y0) = pos[i]
        (x1, y1) = pos[j]
        (x2, y2) = pos[k]
        tri = plt.Polygon([ [ x0, y0 ], [ x1, y1 ], [ x2, y2 ] ],
                          edgecolor = 'white', facecolor = plt.cm.Blues(0.6),
                          zorder = 2, alpha=0.4, lw=0.5)
        ax.add_patch(tri);
        
    # AGGIUNTA MIA 
    # HIGHLIGHTED EDGES
    if markedEdges is not None:
            for i,j in markedEdges:
                (x0 , y0) = pos[i]
                (x1 , y1) = pos[j]
                line = plt.Line2D([ x0, x1 ], [y0, y1 ],color = u'#ff7f0e', zorder = 3, lw=2)
                ax.add_line(line);

    # Drawing the nodes 
    #for i in nodes:
    for i in range(len(pos)):
        (x, y) = pos[i]
        #  radius was 0.1
        circ = plt.Circle([ x, y ], radius = 0.03, zorder = 4, lw=0.5,
                          edgecolor = 'Black', facecolor = u'#ff7f0e')
        ax.add_patch(circ);
        
#     for i in nodes:
        
#         (x, y) = posLabels[i]
#         string = ' $x_{' + str(i) + '}$'
        
#         ax.text(x,y,string, color='black')
        

    return fig

def oldLow(v):
    '''Low function '''

    ind = np.where(v == 1)
    try:
        low = np.max(ind)
        return low
    except Exception as e:
        return -1
    
    
def circumsphereRad(P):
            '''
            TAKEN FROM YOHAI REHANI
            :param P: array Nxd of N points of dimension d
            :return: center and squared radius of the circum-sphere of P
            '''
            p1 = P[0, :]
            A = P[1:, :] - p1
            Pnorm = np.sum(np.power(P, 2), 1)
            b = 0.5*(Pnorm[1:] - Pnorm[0])
            invA = np.linalg.pinv(A)
            c0 = invA.dot(b)
            F = sp.linalg.null_space(A)
            if F.size != 0:
                z = np.transpose(F).dot(p1-c0)
                c = c0 + F.dot(z)
            else:
                c = c0
            R = np.sum(np.power(p1-c, 2))
            return R
        
### COMPUTE STUFF

data_folder = os.path.join(os.getcwd(), 'data', 'compScans')
subsamples_folder = os.path.join(os.getcwd(), 'saved', 'scans', 'Alpha', '1800', 'Subsamples')
output_folder = os.path.join(os.getcwd(), 'saved', 'scans', 'Alpha', '1800', 'Scaffolds')
picture_folder = os.path.join(os.getcwd(), 'saved', 'scans', 'Alpha', '1800', 'Pictures')

# data_folder = os.path.join(os.getcwd(), 'data', 'SmallTest')
# subsamples_folder = os.path.join(os.getcwd(), 'saved', 'SmallTest', 'Subsamples')
# output_folder = os.path.join(os.getcwd(), 'saved', 'SmallTest', 'Scaffolds')
# picture_folder = os.path.join(os.getcwd(), 'saved', 'SmallTest', 'Pictures')

# subsamples_folder = os.path.join(os.getcwd(), 'saved', 'scans','Short','Test1800', 'Subsamples')
# output_folder = os.path.join(os.getcwd(), 'saved', 'scans','Short','Test1800', 'Scaffolds')
# picture_folder = os.path.join(os.getcwd(), 'saved', 'scans','Short','Test1800', 'Pictures')

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)
os.makedirs(subsamples_folder, exist_ok=True)
os.makedirs(picture_folder, exist_ok=True)

# Process each file in the "Subsamples" folder
for filename in os.listdir(subsamples_folder):
    # if filename.endswith('.csv') and filename.startswith('0101'):
    if filename.endswith('.csv'):
        
        print('Processing: ',filename)
        
        file_path = os.path.join(subsamples_folder, filename)        
        name, _ = os.path.splitext(filename)
        
        print('File = '+name, flush=True)
        
        # if it's not already been processed
        #if not os.path.exists(os.path.join(picture_folder,name+'.png')):
            
        print('Process', flush=True)

        # Read
        Points = pd.read_csv(file_path)
        Points = Points.to_numpy()[:,1:3] # Turn to numpy

        # print('Subsample', flush=True)
        # #Sample
        # NPoints = 1800
        # sample = np.random.randint(Points.shape[0], size=NPoints)
        # sample = np.unique(sample)

        # It was already subsampled
        DataSample = Points.copy()
        
        print('Number of points: '+ str(DataSample.shape), flush=True)

        # # Save subsampled point cloud
        # pd.DataFrame(DataSample).to_csv(os.path.join(subsamples_folder, name + '-1800.csv'))
        # # save indices of subsampling 
        # np.save(os.path.join(subsamples_folder, name + '-Indices'),sample)


        print('Build complex', flush=True)
        # Process
        
        
#         # THIS IS WRONG!! CORRECT TO CONSIDER ONLY EDGES IN THE DELAUNAY COMPLEX
        D = scipy.spatial.distance_matrix(DataSample, DataSample, p=2)
        
#         # THIS IS WRONG!! CORRECT TO CONSIDER ONLY EDGES IN THE DELAUNAY COMPLEX
#         D = np.zeros(D.shape)
#         for i in range(DataSample.shape[0]):
#             for j in range(DataSample.shape[0]):
                
#                 dist = np.sqrt(circumsphereRad(DataSample[[i,j]]))
#                 D[i,j] = dist
#                 D[j,i] = dist
        
        Thresh = 3.7
        #Thresh = 0.3
        #Thresh = .28
        #Thresh = .6103
        #Thresh = 0.52
        
        isItShort = False

#         fD = Hopes.filterMatrix(D, T=Thresh)
#         EL = getEL(fD)
        
#         EL = getEL(D)

        #TriList, TriEdges = AlphaComplex(DataSample, Thresh, True) # give the short edges
        
        #_ , EL = AlphaComplex(DataSample, Thresh, False) # give the full Alpha complex edges
        
        TriList, EL, TriScales, EdgeScales = AlphaComplex(DataSample, Thresh)
        
        EL = [x.tolist() for x in EL]

        simplices = EL + TriList
        # ShortSimplices = TriEdges + TriList
        # if isItShort:
        #     simplices = ShortSimplices
        # else:
        #     simplices = AlphaSimplices


        epsList = [Thresh]
        
        # remove "-1800" from name
        name = name[:-5]
        
        fig = draw_2d_simplicial_complex(simplices, pos=DataSample, markedEdges=None, fig=None)
        fig.savefig(os.path.join(picture_folder, name+'Complex.png'))
        

        print('Begin reading- ' + str(datetime.datetime.now()), flush=True)

        start = datetime.datetime.now()

        ## IF YOU COMPUTE THE GENS
        Filtr, cycles, eL = DriverN.getFiltrBasisN(DataSample, epsList, method='Alpha', shortEdges = isItShort)
        
        # IF YOU READ THE MIN GENS
        # cycles = np.load(os.path.join( output_folder, name + '.npy' ))
        # eL = np.load(os.path.join( output_folder, name + 'EdgeList.npy' ))

        end = datetime.datetime.now()

        print('Duration - ' + str(end-start), flush=True)

        ## IF YOU READ THE MIN GENS
        ##cycles = cycles.tolist()
        #cycles = [ cycles[:,i].tolist() for i in range(cycles.shape[1]) ]
        
        # IF YOU COMPUTE THE MIN GENS
        cycles = cycles.as_list()
        cycles = [ x[0].tolist() for x in cycles ]
        eL = np.array(eL)
        
        print('The shape of cycles is ',str(len(cycles[0])) , " x ",str(len(cycles)))
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
        
        fig = draw_2d_simplicial_complex(simplices, pos=DataSample, markedEdges=mEdges, fig=None)
        fig.savefig(os.path.join(picture_folder, name+'.png'))
        print('Saved pics', flush=True)
        
        #print(eL)

        ## PERSISTENCE PAIRS 
        
        print('Computing persistence intervals of basis cycles', flush=True)
        

        PersIntervals = np.zeros((len(cycles), 2))
        
        for i in range(len(cycles)):
            PersIntervals[i,1] = np.infty
            
            
        AllTriList, FullEL , AllTriScales, _ = AlphaComplex(DataSample, np.infty)
        
        
        
#         Del = Delaunay(DataSample)
    
#         AllSimpList = sorted( Del.simplices , key = lambda t : circumsphereRad(DataSample[t]) )
        
#         AllTriList = [ list(sorted(x)) for x in AllSimpList if len(x)==3] # All triangles
#         AllTriScales = [ np.sqrt(circumsphereRad(DataSample[t])) for t in AllTriList ] # with their values
        
        # StepEdgeList = [ list(sorted(x)) for x in AllSimpList if len(x)==2] # All edges
        # StepEdgeScales = [ np.sqrt(circumsphereRad(DataSample[t])) for t in StepEdgeList ] # with their values
        # StepEdgeScales = [x for x in StepEdgeScales if x <= Thresh] # ONLY KEEP UNTIL THRESH
        # StepEdgeList = StepEdgeList[ 0 : len(StepEdgeScales) ]
        
        print('How many scales? ', len(AllTriScales), flush=True)
        
        #StepTriList = [ x for x in AllTriList if circumsphereRad(DataSample[x]) <= Thresh**2 ]
        
        #FullEL = []
        
        if isItShort: # if short alpha complex, the edges are just the boundaries of AllTriList

            for t in AllTriList:
                FullEL.append( (t[0] ,t[1]) )
                FullEL.append( (t[0] ,t[2]) )
                FullEL.append( (t[1] ,t[2]) )
                
            # Also, it could be this way 
            # _, FullEl = AlphaComplex(DataSample, AllTriScales[-1]+1.0, short=True)
                
        else: # if instead it's the Alpha complex, you need to consider all possible simplices
            # _, FullEL = AlphaComplex(DataSample, AllTriScales[-1]+1.0, short=False)
            # FullEL = [tuple(x.tolist()) for x in FullEL]
            pass

        FullEL = [tuple(x.tolist()) for x in FullEL]
        FullEL = list(set(FullEL))
        FullEL = [list(x) for x in FullEL]
        
        FullEL = list(sorted(FullEL))
        
        print('Len of the full EdgeList = ', len(FullEL), flush=True)
        
        MapEdges = [FullEL.index(t.tolist()) for t in eL]
        
        print('MapEdges = ',len(MapEdges), flush=True)
        
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
            
        print('Full boundary matrix computed ', flush=True)
        #print(FullBMatrix)
        
        def sumMod2( v , w):

            return np.mod(v+w,2)

        def reduceCol(ii,mat):

            #print('Reducing column ',ii)

            l = low(mat[:,ii])
            j = 0

            #print('  Low of column ',ii, ' is ',l)
            while True:

                if j == ii:
                    #print('Reached last column, reduction finished')
                    break

                if low( mat[:,j] ) == l and l != -1:

                    #print(' DEPENDENCE, j = ',j, 'low = ', l)
                    mat[:,ii] = sumMod2( mat[:,ii], mat[:,j]) 

                    #print('New column is')
                    #print(mat[:,ii])

                    l = low(mat[:,ii])
                    #print('New low is ',l)

                    j = 0

                else:
                    j += 1
                    #print('Next ',j)


            #print(not np.any(mat[:,i]))
            #return not np.any(mat[:,ii])
            return l == -1
                          

        #Translate cycles in the larger edgelist
        Cycles = np.zeros((len(FullEL), len(cycles)), dtype=int)
        
        print('Translating cycles ', flush=True)
        for i,c in enumerate(cycles):
            
            # TRANSLATE CYCLES IN THE LARGER EDGELIST
            
            for jj, entry in enumerate(c):
                
                if entry > 0:
                    
                    # Map to the full edgelist
                    LargerIndex = MapEdges[jj]
                    Cycles[LargerIndex,i] = 1
                    
                

            #edges = [ FullEL[j] for j,v in enumerate(Cycles[:,i].tolist()) if v!=0  ]
            #birthVal = max([ D[e[0] , e[1]] for e in edges  ])
            
            #birthVal = max([np.sqrt(circumsphereRad(DataSample[[e[0], e[1]]])) for e in edges])

            #PersIntervals[i,0] = birthVal

        # Turn to list
        Cycles = [ Cycles[:,i] for i in range(Cycles.shape[1]) ]
        
        #print(AllTriScales)
        # The death time is the smallest scale at which sufficiently many triangles appear to fill a cycle
        # This will require several calls of a linear system
        
        # Initialize up until Thresh

        firstCols = len( [x for x in AllTriScales if x <= Thresh] )

        BMatrix = np.zeros((len(FullEL), firstCols), dtype=int)

        # Since it is slicing, it MAKES A COPY and there is no referencing problem, i.e. changes to BMatrix do not affect FullBMatrix
        BMatrix[:, :firstCols] = FullBMatrix[:, :firstCols]

        print('Initial boundary matrix computed ', flush=True)
        #print(BMatrix)

        # reduce the initial boundary matrix 
        print('Start reducing initial matrix (' + str(firstCols)+' columns) ' , flush=True)
        for h in range(firstCols):

            reduceCol(h, BMatrix)
            if h % 20 == 0:
                print('Done ', str(h), ' out of ', str(firstCols), flush=True)

        print('Finished reducing initial matrix', flush=True)
        
        ## SAVE this matrix for later! We need it to compute birth times
        import copy
        ReducedTriBMatrix = copy.deepcopy(BMatrix)

        steps = [x for x in AllTriScales if x >= Thresh]

        LenSteps = len(steps)
        
        CyclesToKill = [x for x in range(len(Cycles))]

        for stepInd, s in enumerate(steps):

            print('Working step ', stepInd, ' of ', LenSteps, flush=True)
            #print('Step ',stepInd, ', reducing')
            BMatrix = np.c_[BMatrix, FullBMatrix[:,firstCols+stepInd]] # add new column

            #print('Bmatrix is now ')
            #print(BMatrix)

            reduceCol( firstCols+stepInd , BMatrix ) # reduce the last column added

            #print('After reduction, Bmatrix is now ')
            #print(BMatrix)

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
                        print('Death of cycle',i,' found at step ' +str(stepInd)+' and scale ' + str(s), flush=True)
                        PersIntervals[i,1] = s

                        CyclesToKill.remove(i)

                    if len(CyclesToKill) == 0:
                        break
                
            if len(CyclesToKill) == 0:
                print('Reached the end of cycles to kill')
                break
                
        ## NOW BIRTH TIME
        # Birth time initialization
        
        print('Computing birth times of basis cycles', flush=True)
        
        # Compute edges up to step Thresh
#         _ , StepEdgeList, _ , StepEdgeScales = AlphaComplex(DataSample, Thresh)
        
#         # StepEdgeList = list(sorted(StepEdgeList, key = lambda e : circumsphereRad(DataSample[e]) ))
#         # StepEdgeScales = [ np.sqrt( circumsphereRad(DataSample[e]) ) for e in StepEdgeList ]
        
        BirthIntervals = np.zeros((len(cycles), 2))
        for i in range(len(cycles)):
            BirthIntervals[i, 0] = -1  # Initialize birth time to -1 (not yet found)

#         print(len(StepEdgeList) , len(StepEdgeScales))
        
#         # BUILD THE MATRIX OF EDGES WRT THE LARGER EDGELIST!
        
#         # Initialize to shape Number of edges in the larger edgelist x N of edges up to step ThisStep 
#         FullEdgeMatrix = np.zeros((len(FullEL), len(StepEdgeList)), dtype=int)
        
#         for i,e in enumerate(StepEdgeList):
            
#             E = sorted(e) # probably unnecessary
            
#             # not very pythonic
#             EIndex = FullEL.index( [E[0], E[1]] ) # find the index of E in the FullEL
            
#             FullEdgeMatrix[ EIndex , i ] = 1
                           
#         # Initialize the boundary matrix with the reduced triangle boundary matrix
#         # is this already a deep copy? Yes
#         BMatrixBirth = ReducedTriBMatrix
        
#         TrianglesLenAtThisStep = BMatrixBirth.shape[1]
#         NCyclesToBirth = len(Cycles)

#         # Iterate through all scales up to Thresh
#         for stepInd, s in enumerate(StepEdgeScales):
                
#             # ADD to BMatrixBirth the right column from FullEdgeMatrix
            
#             BMatrixBirth = np.c_[ BMatrixBirth , FullEdgeMatrix[ : , stepInd ] ]
            
#             # REDUCE the new column
            
#             reduceCol( TrianglesLenAtThisStep + stepInd , BMatrixBirth) # reduce the last column

#             # Now, check each cycle for dependence
#             for i, cycle in enumerate(Cycles):
#                 if BirthIntervals[i, 0] == -1:  # Only check if birth time hasn't been found
#                     # Create the system matrix with triangles and edges
#                     system = np.c_[BMatrixBirth, cycle]

#                     lastCol = system.shape[1] - 1

#                     if reduceCol(lastCol, system):  # If it goes to zero, it's dependent and then it is born
#                         # FOUND the birth step
#                         BirthIntervals[i, 0] = s  # Record the birth scale
#                         NCyclesToBirth -= 1 # There is one less to compute
#                         print('Birth of cycle', i, 'found at step ' + str(stepInd) + ' and scale ' + str(s), flush=True)
#             if NCyclesToBirth == 0:
#                 break
                
        PersIntervals[:,0] = BirthIntervals[:,0]

        # Output final birth intervals
        print('Birth intervals computed', flush=True)

        #print(len(AllTriList))
        print('Persistence Pairs: ', PersIntervals)
        
        PP_filename = os.path.join(output_folder, name+'PP.npy')
        np.save(PP_filename, PersIntervals)
        
        print('** Square pairs')
        print(PersIntervals**2)
        