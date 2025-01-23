# BUILD THE ALPHA COMPLEX VIA GUDHI

import numpy as np
from gudhi import AlphaComplex as gudhiAlphaComplex

def AlphaComplex( points, thresh ):
    
    # Assuming points is a numpy array
    points = points.tolist()
    
    alpha = gudhiAlphaComplex(points = points)
    stree = alpha.create_simplex_tree()
    
    simplices = list(stree.get_filtration())
    
    TriList = [ np.array(x[0]) for x in simplices if len(x[0])==3 and x[1] <= thresh**2 ]
    EdgeList = [ np.array(x[0]) for x in simplices if len(x[0])==2 and x[1] <= thresh**2 ]
    
    TriScales = [ np.sqrt(x[1]) for x in simplices if len(x[0])==3 and x[1] <= thresh**2 ]
    EdgeScales = [ np.sqrt(x[1]) for x in simplices if len(x[0])==2 and x[1] <= thresh**2 ]
    
    return TriList, EdgeList, TriScales, EdgeScales
    

