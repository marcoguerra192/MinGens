# Plotting function 

import os 
os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"
import matplotlib.pyplot as plt

import itertools

def draw_2d_simplicial_complex(simplices, pos=None, return_pos=False, fig = None, markedEdges=None):
    """
    Draw a simplicial complex up to dimension 2 from a list of simplices
        
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
