#cython: language_level=3
#cython : profile = True

# -*- coding: utf-8 -*-
# TESI: CLASSE per il problema SHORTEST PATH TREE
#
# FAST VERSION (2026-11). Same algorithm and same output as the original
# CySPT.pyx (kept as Old_CySPT.pyx). Changes, all implementation-level:
#   1. Edge lookup by dictionary instead of list.index (getHeatMap,
#      closeCycle): O(1) per edge instead of O(#edges).
#   2. closeCycle builds the cycle from the list of its edge indices in one
#      step, instead of summing one dense heat map per edge.
#   3. The Z2 label of each sentinel edge with respect to a support vector,
#      <Sup, ann(e)> mod 2, is computed once for all sentinel edges with one
#      matrix product, instead of once per use with a NumPy dot product.
#   4. computeLabels uses a deque for its breadth-first traversal and works
#      on a Python list per root, instead of popping from the front of a
#      list and writing NumPy scalars one by one.
#   5. The orthogonality test of findShortestNonOrtho is evaluated for all
#      (root, non-tree edge) candidates at once with NumPy, in the same
#      order as before, so the chosen cycle is the same.
# fixNTE is deliberately unchanged: it fixes the order of the candidates,
# and that order decides ties between equally short cycles.
#
# fast-v3 (2026-11). One correctness fix and further speed-ups, same
# outputs otherwise (see CHANGES.md):
#   6. inv_z2 inverts over Z2 exactly, by Gauss-Jordan elimination. The
#      previous version used A^-1 * det(A) from floating point, which loses
#      the integers once det(A) passes 2^53: wrong inverses from ~50 x 50,
#      i.e. H1 of dimension ~100 at one scale.
#   7. UpdateMehlhorn recomputes only the rows of B = L C that change.
#   8. mulz2 multiplies in float64 (BLAS); the entries are small integers,
#      exact in float64.
#   9. getAnnotation reduces rows stored as bit vectors (Python ints). It
#      accepts the dense d2 as before, or the rows as bit vectors directly.
#  10. Candidate cycle lengths come from tree distances, computed once per
#      graph (d(r,u) + d(r,v) + l(uv) - 2 d(r, lca(u,v))); a cycle vector is
#      built only for the shortest candidates.
#  11. computeLabels fills all trees at once, level by level;
#      cycleAnnotation is vectorised; with draws, identical tied cycles
#      (the same cycle reached from different roots) are compared once.
# The sympy import is removed: it was no longer used.

import numpy as np
from math import floor, ceil
import scipy.sparse
from scipy.sparse import csgraph
import functools
import GeometryN as GeometryN
import multiprocessing as mlt
from multiprocessing import Pipe
from collections import deque
import os

# Marker to check which implementation is loaded (the compiled module is
# used in preference to the .pyx, so a rebuild is needed after editing)
CYSPT_IMPLEMENTATION = 'fast-v3'



def gf2_inverse(matrix):
    """ Inverse of a square matrix over Z2, by Gauss-Jordan elimination.

        Exact at any size: rows are 0/1 and adding rows is XOR, so nothing
        grows and nothing is rounded. Raises np.linalg.LinAlgError if the
        matrix is singular over Z2. Returns an int64 array.
    """
    n = matrix.shape[0]
    augmented = np.zeros((n, 2 * n), dtype=np.uint8)
    augmented[:, :n] = np.asarray(matrix, dtype=np.int64) % 2
    augmented[:, n:] = np.eye(n, dtype=np.uint8)
    for column in range(n):
        candidates = np.nonzero(augmented[column:, column])[0]
        if candidates.size == 0:
            raise np.linalg.LinAlgError("matrix is singular over Z2")
        pivot = column + candidates[0]
        if pivot != column:
            augmented[[column, pivot]] = augmented[[pivot, column]]
        rows_to_clear = np.nonzero(augmented[:, column])[0]
        rows_to_clear = rows_to_clear[rows_to_clear != column]
        augmented[rows_to_clear] ^= augmented[column]
    return augmented[:, n:].astype(np.int64)


def rows_to_bits(rows):
    """ Rows of a 0/1 matrix as Python ints, bit k of row r = rows[r, k].

        If `rows` is already a sequence of Python ints (bit vectors), it is
        returned as a list unchanged.
    """
    if len(rows) == 0:
        return []
    if isinstance(rows[0], (int, np.integer)):
        return [int(r) for r in rows]
    dense = np.asarray(rows, dtype=np.uint8)
    packed = np.packbits(dense, axis=1, bitorder='little')
    return [int.from_bytes(r.tobytes(), 'little') for r in packed]


class SimplexGraph:
    """ Classe che implementa un grafo come matrice di adiacenza pesata, e fornisce metodi per ottenere lo shortest path tree
        a partire da un nodo arbitrario, più l'elenco dei non-tree edges.
        Input costruttore: Il numero di vertici, o il numero di vertici e la matrice dei pesi
    """

    def __init__(self, numVertices, weights, maximal=None, connection=None, eps=None):
        self.NVert = numVertices
        self.Weights = scipy.sparse.coo_matrix(weights)
        # coo_matrix può essere istanziata a partire da una matrice densa di tipo np.ndarray
        self.Weights = self.Weights.tocsr()
        # al momento ci serve che sia csr per accedere ai pesi

        self.dist = None
        self.pred = None
        self.Followers = []
        self.Sentinel = None
        self.SentinelSet = None
        self.CCC = None
        self.CCV = None
        self.Origins = None
        self.Ann = None
        self.dimH1 = None
        self.dimB1 = None
        self.dimZ1 = None
        self.Labels = -1 * np.ones((self.NVert, self.NVert), dtype=int)
        self.Support = []  # vettori di supporto
        self.SHB = None  # base minima di omologia
        self.NTE = None  # array di liste di NTE
        self.closedCycles = {}  # dizionario di cicli già chiusi
        # la chiave è la terna (source, v1, v2)
        self.allDraws = True
        self.cycleBasis = []

        # multiprocessing pipe object to return partial results
        self.pipe = connection
        self.eps = eps

        # matrici per il calcolo con la versione ricorsiva
        self.L = None  # matrice dei label m(S_i , a(e_j))
        self.C = None  # matrice m(S_i, C_j)
        self.B = None  # matrice m(S_i, C_j)
        self.S = None  # matrice S_i

        # Structures of the fast version, filled in later:
        # AnnMat (dimH1 x #sentinel) annotation matrix as a plain int array,
        # SentinelIdx: sentinel edge (v1, v2) -> its column in AnnMat,
        # NTE_root, NTE_u, NTE_v, NTE_sent: the candidates as parallel arrays.
        self.AnnMat = None
        self.SentinelIdx = None
        self.NTE_root = None
        self.NTE_u = None
        self.NTE_v = None
        self.NTE_sent = None
        # Cache of fixNTE: the list of all directed edges (as returned by
        # scipy.sparse.find) and, per connected component, the set of those
        # edges whose first vertex lies in the component.
        self._edgeW = None
        self._compEdgeSet = None

        # diversifichiamo il caso in cui abbiamo la filtrazione dal caso
        # singolo step
        if (maximal == None):
            rw, cl, _ = scipy.sparse.find(self.Weights)  # righe e colonne sparse
            self.Edges = list(zip(list(rw), list(cl)))  # lista di coppie di vertici
            # scelgo solo quelli da un vertice più basso ad uno più alto
            self.Edges = [x for x in self.Edges if x[0] < x[1]]
            # funziona perchè non ci sono self loops

            def orderEdges(edgeList):
                """ Funzione che ordina la lista di edge secondo qualche criterio. Per ora il criterio è prima rispetto al primo vertice,
                poi rispetto al secondo.
                """
                return sorted(edgeList, key=lambda t: (t[0], t[1]))
            # ora è una lista ordinata di edges, senza ripetizioni!
            self.Edges = orderEdges(self.Edges)
            # può essere usata per tenere traccia dell'ordine degli edges
            self.NEdges = len(self.Edges)  # numero di edges IN SENSO INDIRETTO
            self.WEdges = []  # lista delle lunghezze (pesi) degli edges
            for e in self.Edges:
                self.WEdges.append(self.Weights[e[0], e[1]])

            self.WEdges = np.array(self.WEdges)  # trasformo in array
            self.NonZ = int(self.Weights.count_nonzero() / 2)
        else:  # se ho definito la matrice più fine della filtrazione
            # maximal è una tupla = (edgelist, weightlist)
            # la edgelist deve già essere ordinata
            self.Edges = maximal[0]
            self.NEdges = len(self.Edges)
            self.WEdges = maximal[1]
            self.NonZ = int(self.Weights.count_nonzero() / 2)

        # We index the edges once: (v1, v2) with v1 < v2, as Python ints,
        # maps to the position of the edge in self.Edges.
        self.EdgeIndex = {}
        for k, e in enumerate(self.Edges):
            self.EdgeIndex[(int(e[0]), int(e[1]))] = k

    def getEdgeList(Weights):
        """ Metodo Statico per ottenere l'edgelist da una matrice massimale
            In input solo la matrice dei pesi
        """
        rw, cl, _ = scipy.sparse.find(Weights)  # righe e colonne sparse
        edgelist = list(zip(list(rw), list(cl)))  # lista di coppie di vertici
        # scelgo solo quelli da un vertice più basso ad uno più alto
        edgelist = [x for x in edgelist if x[0] < x[1]]
        # funziona perchè non ci sono self loops

        def orderEdges(edgeList):
            """ Funzione che ordina la lista di edge secondo l'ordine lessicografico
            """
            return sorted(edgeList, key=lambda t: (t[0], t[1]))
        edgelist = orderEdges(edgelist)  # ora è una lista ordinata di edges, senza ripetizioni!
        WEdges = []  # lista delle lunghezze (pesi) degli edges
        for e in edgelist:
            WEdges.append(Weights[e[0], e[1]])
        WEdges = np.array(WEdges)  # trasformo in array

        return (edgelist, WEdges)

    def getEdgeVector(self):  # ordinamento degli archi
        return self.Edges

    def edgeIndex(self, e):
        """ Position of edge e = (a, b) in self.Edges, in either orientation.
            Raises ValueError if the edge does not exist (as getHeatMap did).
        """
        a = int(e[0])
        b = int(e[1])
        if a > b:
            a, b = b, a
        try:
            return self.EdgeIndex[(a, b)]
        except KeyError:
            print('Wrong Edge!', e)
            errstr = "L'edge ", e, " non esiste!"
            raise ValueError(errstr)

    def getHeatMap(self, edgeList):
        """ Prende in input una lista di coppie di vertici (v1, v2), li ordina in modo che v1 < v2, e restituisce un vettore 0-1 lungo             quanto la dimensione dello spazio vettoriale delle 1-catene, che descrive la 1-catena corrispondente.
            Se l'input contiene un edge inesistente o un self-loop, solleva ValueError
        """
        heatMap = np.zeros(self.NEdges, dtype=int)  # vettore di zeri
        # Duplicates set the same entry to 1 twice, as the original's set did
        for val in edgeList:
            heatMap[self.edgeIndex(val)] = 1
        return heatMap

    def set_weights_by_coo(self, data, rows, col):
        """ Istanzia la matrice dei pesi nel formato sparso attraverso tre array data, rows e col"""
        self.Weights = scipy.sparse.coo_matrix((data, (rows, col)), shape=(self.NVert, self.NVert))

    def computeShortestPath(self, withDistances=False):
        """ Restituisce la matrice dei predecessori nello Shortest Path Tree. Ogni riga si riferisce ad un nodo
            di origine dello SPT. Se l'argomento opzionale è vero ritorna anche la matrice delle distanze
        """
        # va convertito in un formato su cui può fare i calcoli
        (self.dist, self.pred) = csgraph.shortest_path(
            self.Weights.tocsr(), directed=False, return_predecessors=True)

        # Aggiungiamo il calcolo dei successori. Utile da fare una volta invece di n
        # Viene utile una struttura di lista^3
        for source in range(self.NVert):  # per ogni sorgente
            link = self.pred[source, :]  # predecessori
            forward = [[] for i in range(self.NVert)]  # lista di liste vuote lunga NVert
            for (foll, prev) in enumerate(link):
                if (prev != -9999):  # se è connesso e non è source
                    forward[prev].append(foll)
            self.Followers.append(forward)

    def shortestPathTree(self, source, withD=False, justEdges=False):
        """ Restituisce la matrice dello shortest path tree centrata nel vertice source. Se withD è vero
            la matrice è pesata, altrimenti è solo 0-1. Se source non è un vertice valido solleva un'eccezione
            I VERTICI VANNO DA 0 A N-1!!
            Se justEdges == True restituisce solo i non tree edges
            Input: L'indice del nodo sorgente, opzionale un bool per le distanze
            Output: La matrice sparsa del grafo SPT
        """
        if (type(source) is not int or source < 0 or source >= self.NVert):  # da qui in poi source è un valore affidabile
            raise Exception('Indice del vertice non valido!')

        if withD:
            if (self.dist is None or self.pred is None):
                _ = self.computeShortestPath(withD)
        else:
            if self.pred is None:
                _ = self.computeShortestPath(withD)

        # estraiamo dalla matrice dei precedenti la riga che ci serve
        link = np.array(self.pred[source, :])
        # predisponiamo i vettori data, rows e col per creare la matrice del SPT
        data = []
        rows = []
        col = []  # man mano faremo append dei dati in queste liste
        # scorriamo il vettore e creiamo i link
        for i in range(len(link)):
            prev = link[i]
            if (prev == -9999):  # non appartiene alla stessa componente connessa => non ha predecessori
                pass
            else:
                rows.append(prev)
                col.append(i)  # aggiungi un edge dal predecessore del nodo i-esimo al nodo i-esimo
                rows.append(i)  # devo aggiungere anche il simmetrico!
                col.append(prev)
                if withD:  # se devo aggiungere le distanze le cerco in self.dist
                    data.append(self.dist[prev, i])
                    # aggiungo anche il simmetrico (stessa distanza)
                    data.append(self.dist[prev, i])
                else:  # se no inserisco solo 1 per indicare l'edge fra prev e i
                    data.append(1)
                    data.append(1)

        if withD:
            data = np.array(data, dtype=float)
        else:
            data = np.array(data, dtype=int)
        rows = np.array(rows)
        col = np.array(col)
        SPT = scipy.sparse.coo_matrix((data, (rows, col)), shape=(self.NVert, self.NVert))
        # crea matrice del grafo
        # Per trovare i non-tree edge: trovare gli elementi di Weights che NON SONO in SPT
        # restituisce vettori di righe, colonne e valori non nulli
        rw, cl, _ = scipy.sparse.find(self.Weights)
        # ignoro il vettore di dati, che non ci interessa
        edgeSPT = list(zip(list(rows), list(col)))  # creo una lista di edge per SPT
        edgeW = list(zip(list(rw), list(cl)))  # stessa cosa per il grafo originale
        # MA CONTIENE ANCHE TUTTI GLI EDGES NELLE COMPONENTI NON CONNESSE!
        # Affinchè il metodo funzioni, è necessario che il nodo source non abbia -9999 per distinguerlo da quelli davvero disconnessi
        link[source] = source  # creiamo un "finto" self-loop
        edgeW = [x for x in edgeW if link[x[0]] != -9999]
        # limito a solo gli edge fra nodi che sono raggiungibili da source
        NTE = list(set(edgeW) - set(edgeSPT))  # differenza insiemistica tra i set di edges
        # contiene ancora entrambe le versioni dell'edge (v1,v2 e v2,v1). Imponiamo solo quella (v1,v2) con v1<v2
        NTE = [x for x in NTE if (x[0] < x[1])]

        if justEdges:
            return NTE
        else:
            return (SPT, NTE)

    def spanningTree(self):
        """ Calcola uno spanning tree a partire dalla stessa matrice dei predecesori di SPT. Questa volta, oltre a lavorare su tutte le               componenti connesse, restituisce tutti i sentinel edges.
        """
        if self.pred is None:
            _ = self.computeShortestPath(True)
        # prepariamo le liste per la matrice del grafo
        rows = []
        col = []
        data = []
        self.CCV = list(np.zeros(self.NVert, dtype=int))  # lista indice delle componenti connesse
        self.CCC = 0  # contatore delle componenti connesse
        self.Origins = []  # lista lunga come il numero di CC. In posizione k contiene l'origine della componente connessa k+1-esima
        # che sarà necessario utilizzare per calcolare i sentinel edges!
        while (0 in self.CCV):  # finchè non ho assegnato ogni vertice ad una componente
            self.CCC += 1
            vert = self.CCV.index(0)  # trova il primo vertice non assegnato ad alcuna componente
            # Questo vertice diventerà l'origine della CC!! Sarà necessario usare lui per calcolare i sentinel edges!
            self.Origins.append(vert)  # questa CC origina da vert
            self.CCV[vert] = self.CCC  # assegna questo vertice alla nuova componente
            link = np.array(self.pred[vert, :])  # estrai i predecessori connessi a vert
            for i in range(len(link)):
                prev = link[i]
                if (prev == -9999):  # se è disconnesso fai nulla
                    pass
                else:
                    self.CCV[i] = self.CCC  # se è connesso assegna lo stesso indice di componente
                    rows.append(prev)
                    # aggiungi un edge dal predecessore del nodo i-esimo al nodo i-esimo
                    col.append(i)
                    rows.append(i)  # devo aggiungere anche il simmetrico!
                    col.append(prev)
                    data.append(1)
                    data.append(1)
        data = np.array(data, dtype=int)
        rows = np.array(rows)
        col = np.array(col)
        ST = scipy.sparse.coo_matrix((data, (rows, col)), shape=(
            self.NVert, self.NVert))    # matrice del grafo Spanning Tree
        # per trovare i sentinel edges. Questo spanning tree è già generato a partire dall'origine di ogni CC!
        rw, cl, _ = scipy.sparse.find(self.Weights)  # tutti gli edge del grafo originale
        edgeW = list(zip(list(rw), list(cl)))  # lista edge grafo originale
        edgeST = list(zip(list(rows), list(col)))  # lista edge spanning tree
        self.Sentinel = list(set(edgeW) - set(edgeST))  # differenza insiemistica
        # contiene ancora entrambe le versioni dell'edge (v1,v2 e v2,v1). Imponiamo solo quella (v1,v2) con v1<v2
        self.Sentinel = [x for x in self.Sentinel if (x[0] < x[1])]
        self.Sentinel = sorted(self.Sentinel, key=lambda x: (x[0], x[1]))

        # We also keep the sentinel edges as a set of int pairs, for the
        # membership check in closeCycle (a list search before)
        self.SentinelSet = set((int(x[0]), int(x[1])) for x in self.Sentinel)
        return (ST, self.Sentinel, self.CCV)

    def closeCycle(self, E, force=None):
        """
        Funzione che prende in input un sentinel edge E, e sfruttando la matrice
        self.pred e l'ordinamento self.Edges, restituisce un vettore dello spazio
        delle 1-catene che descrive il ciclo identificato da E nello Shortest Path Tree.
        Lo SPT è quello radicato nel nodo ORIGINE DELLA CC DI E!
        Se l'input non è quello atteso, solleva ValueError.
        Se è definito force, forza il nodo sorgente ad essere quello. Utile per
        generare il candidate set.
        Sfrutta la struttura self.closedCycles per verificare se una certa chiusura
        è già stata calcolata

        Fast version: we collect the indices of the edges of the cycle (E,
        then the two tree paths back to the source) and build the 0-1 vector
        once, taking each count mod 2 so that edges shared by the two paths
        cancel, as the original sum of heat maps did.
        """

        if (force is not None):  # se è fissata l'origine
            key = (force, E[0], E[1])  # la chiave è (source, v1,v2)
            try:  # trova se è già stato calcolato
                cycle = self.closedCycles[key]
                return cycle
            except:
                pass

        try:
            # index of the closing edge E; raises if E does not exist
            edgeIdx = [self.edgeIndex(E)]
        except Exception:
            raise Exception

        # controllo che sia un sentinel edge (solo se force non è definito)
        if (self.SentinelSet is not None and force is None):
            a = int(E[0])
            b = int(E[1])
            if a > b:
                a, b = b, a
            if (a, b) not in self.SentinelSet:
                raise ValueError("E non è un sentinel edge!")

        # fissiamo i bordi del sentinel edge
        v1 = E[0]
        v2 = E[1]

        if (force is not None):
            if (force < 0 or force >= self.NVert):
                raise ValueError("Parametro force non valido!")
            else:
                source = force
        else:
            # troviamo il NODO ORIGINE della CC di E
            source = self.Origins[self.CCV[v1] - 1]
            key = (source, E[0], E[1])  # calcoliamo la chiave del ciclo da chiudere
        # Calcoliamo i pred, se non è già stato fatto
        if (self.pred is None):
            self.pred = self.computeShortestPath(False)

        link = self.pred[source, :]  # vettore dei precedenti rispetto a source
        pr = link[v2]
        if (pr == -9999 and v2 != source):  # significa che v2 è irraggiungibile da source, ci deve essere stato un errore
            raise ValueError(
                'Critical Error! Source e v2 dovrebbero essere per costruzione nella stessa CC')
        follow = v2
        while (pr != -9999):  # segue all'indietro i precedenti verso source (che ha link -9999)
            edgeIdx.append(self.edgeIndex((pr, follow)))
            follow = pr
            pr = link[pr]

        pr = link[v1]
        if (pr == -9999 and source != v1):  # significa che v1 è irraggiungibile da source, sono in due componenti connesse diverse
            raise ValueError(
                'Critical Error! Source e v1 dovrebbero essere per costruzione nella stessa CC')
        follow = v1
        while (pr != -9999):  # finchè non torno a source
            edgeIdx.append(self.edgeIndex((pr, follow)))
            follow = pr
            pr = link[pr]

        # Count how many times each edge appears, then keep it iff odd
        cycle = np.zeros(self.NEdges, dtype=int)
        np.add.at(cycle, edgeIdx, 1)
        cycle = cycle % 2

        # aggiorniamo il dizionario dei cicli chiusi
        self.closedCycles[key] = cycle

        return cycle

    def getCycleBase(self):
        """ Calcola una base dei cicli a partire dai sentinel edges. Per ottenere l'annotazione
        """
        # se ancora non l'ha fatto calcola i sentinel edges
        if (self.Sentinel is None):
            _ = self.spanningTree()

        for e in self.Sentinel:
            self.cycleBasis.append(self.closeCycle(e))

    def getAnnotation(self, d2):
        """
        d2: border matrix of the set of faces of the simplicial complex, borders are written as rows.
            Either a 0/1 matrix (one row per triangle, one column per edge
            of self.Edges), or a list of the rows as bit vectors (Python
            ints, bit k = edge k), which avoids the dense matrix.
        Z: basis of 1-dim-cycles. Each row is a cycle obtained from a spanning tree and a sentinel edge.
        returns an annotation of sentinel edges. Non sentinel edges have annotation equal to 0.

        fast-v3: same reduction, in the same order, so the same annotation:
        boundary rows first, then the fundamental cycles, each reduced by
        the rows whose lowest entry (largest edge index) it shares. Rows are
        bit vectors and XOR adds them over Z2, instead of dense rows with a
        deletion per zero row.
        """
        # pivots maps a lowest edge index to the reduced row that has it,
        # and, for reduced cycle rows, to the index of that cycle.
        pivots = {}

        # The boundaries. A row that reduces to zero is dependent and dropped.
        dimB1 = 0
        for row in rows_to_bits(d2):
            low = row.bit_length() - 1
            while row and low in pivots:
                row ^= pivots[low][0]
                low = row.bit_length() - 1
            if row:
                pivots[low] = (row, None)
                dimB1 += 1

        # The fundamental cycles. For each we record which earlier reduced
        # cycle rows were added to it (boundaries are zero in H1).
        cycleRows = rows_to_bits(np.asarray(self.cycleBasis, dtype=np.uint8)) if len(self.cycleBasis) else []
        used = []
        pivotCycles = []
        for index, row in enumerate(cycleRows):
            combo = 0
            low = row.bit_length() - 1
            while row and low in pivots:
                pivotRow, pivotCycle = pivots[low]
                row ^= pivotRow
                if pivotCycle is not None:
                    combo ^= 1 << pivotCycle
                low = row.bit_length() - 1
            if row:
                # an independent class: the cycle is a basis element of H1
                pivots[low] = (row, index)
                pivotCycles.append(index)
                combo ^= 1 << index
            used.append(combo)

        # A[c, j]: coordinate of fundamental cycle c on the j-th independent cycle
        A = np.zeros((len(cycleRows), len(pivotCycles)), dtype=np.int64)
        for j, p in enumerate(pivotCycles):
            for c, combo in enumerate(used):
                A[c, j] = (combo >> p) & 1

        """observation: the number of rows of A is the dimension of the 1-dim-cycle group;
           the number of columns is the dimension of the 1st homology group
        """
        self.dimB1 = dimB1
        self.dimZ1 = np.shape(A)[0]
        self.dimH1 = np.shape(A)[1]
        self.Ann = np.matrix(A).transpose()

        self.AnnDict = {}
        for i, e in enumerate(self.Sentinel):
            ann = self.Ann[:, i]
            self.AnnDict[e] = ann

        # Fast version: the same annotations as a plain (dimH1 x #sentinel)
        # int array, with the column of each sentinel edge in SentinelIdx.
        # Column j of AnnMat is the annotation of self.Sentinel[j].
        self.AnnMat = np.asarray(self.Ann, dtype=np.int64)
        self.SentinelIdx = {}
        for j, e in enumerate(self.Sentinel):
            self.SentinelIdx[(int(e[0]), int(e[1]))] = j

        # Column of each edge of self.Edges in AnnMat (-1: not a sentinel)
        self.edgeSentinelColumn = np.array(
            [self.SentinelIdx.get((int(a), int(b)), -1) for a, b in self.Edges], dtype=np.int64)

        return (self.Ann, self.dimB1, self.dimZ1, self.dimH1)


    def innerProd(self, S, C):
        """ Scalar Product over Z2 of support vector S and the ANNOTATION of cycle C
        """
        return np.mod(np.dot(np.array(S).transpose(), np.array(C)), 2)

    def sentinelBits(self, Sup):
        """ Z2 label of every sentinel edge with respect to support vector Sup.

            Entry j is <Sup, annotation of self.Sentinel[j]> mod 2, i.e. the
            value innerProd(Sup, AnnDict[Sentinel[j]]) computed for all j at
            once. Returns a 1-D int array of length #sentinel edges.
        """
        sup = np.asarray(Sup, dtype=np.int64).ravel()
        bits = sup.dot(self.AnnMat) % 2
        return np.asarray(bits, dtype=np.int64).ravel()

    def computeLabels(self, Sup):
        """ Calcola i labels rispetto al support vector S_i Sup

            fast-v3: the label of a node is the label of its parent plus the
            bit of the edge between them (0 for non-sentinel edges), mod 2;
            the root has label 0. All trees are filled at once, level by
            level. Same values as before; entries never visited keep their
            previous value.
        """
        if getattr(self, 'treeLevels', None) is None:
            self._prepareTrees()
        bits = self.sentinelBits(Sup)
        bitsExt = np.append(bits, 0)
        n = self.NVert
        self.Labels[np.arange(n), np.arange(n)] = 0
        for rows, nodes, parents, columns in self.treeLevels:
            self.Labels[rows, nodes] = (self.Labels[rows, parents] + bitsExt[columns]) % 2

    def _prepareTrees(self):
        """ Depths, levels and ancestor tables of all shortest path trees
            (row r of each (n, n) array refers to the tree rooted at r).
            Used by computeLabels and by the candidate lengths.
        """
        n = self.NVert
        parent = np.asarray(self.pred, dtype=np.int64)
        reachable = parent >= 0
        rows, nodes = np.nonzero(reachable)

        # depth of each node in each tree, level by level from the roots
        depth = np.full((n, n), -1, dtype=np.int64)
        depth[np.arange(n), np.arange(n)] = 0
        levels = []
        level = 0
        while True:
            parentDepth = depth[rows, parent[rows, nodes]]
            atLevel = (parentDepth == level) & (depth[rows, nodes] < 0)
            if not np.any(atLevel):
                break
            levelRows = rows[atLevel]
            levelNodes = nodes[atLevel]
            depth[levelRows, levelNodes] = level + 1
            levels.append((levelRows, levelNodes, parent[levelRows, levelNodes]))
            level += 1

        # column in AnnMat of each sentinel edge, as a symmetric matrix
        sentinelColumn = np.full((n, n), -1, dtype=np.int64)
        for (a, b), column in self.SentinelIdx.items():
            sentinelColumn[a, b] = column
            sentinelColumn[b, a] = column
        self.treeLevels = [(r, x, p, sentinelColumn[p, x]) for r, x, p in levels]

        # ancestors[j][r, x]: ancestor of x 2^j levels up in tree r
        # (the root, and nodes outside its component, are their own parent)
        first = np.where(reachable, parent, np.arange(n)[np.newaxis, :])
        ancestors = [first]
        for _ in range(max(1, int(np.ceil(np.log2(max(level, 1)))) + 1)):
            previous = ancestors[-1]
            ancestors.append(np.take_along_axis(previous, previous, axis=1))
        self.treeDepth = depth
        self.treeAncestors = ancestors

    def _lowestCommonAncestor(self, roots, u, v):
        """ Lowest common ancestor of u and v in the tree rooted at each root. """
        a = u.copy()
        b = v.copy()
        depthA = self.treeDepth[roots, a]
        depthB = self.treeDepth[roots, b]
        swap = depthA < depthB
        a[swap], b[swap] = b[swap], a[swap]
        difference = np.abs(depthA - depthB)
        for j, table in enumerate(self.treeAncestors):
            move = ((difference >> j) & 1).astype(bool)
            a[move] = table[roots[move], a[move]]
        for table in reversed(self.treeAncestors):
            differ = table[roots, a] != table[roots, b]
            a[differ] = table[roots[differ], a[differ]]
            b[differ] = table[roots[differ], b[differ]]
        return np.where(a == b, a, self.treeAncestors[0][roots, a])

    def _candidateLengths(self):
        """ Length of every (root r, non-tree edge uv) candidate cycle.

            The candidate is the tree path r -> u, the edge uv and the tree
            path v -> r. The two paths share the stretch r -> lca(u, v),
            which cancels mod 2, so the length is
            d(r, u) + d(r, v) + l(uv) - 2 d(r, lca(u, v)):
            the same value as WEdges . (cycle vector).
        """
        roots = self.NTE_root
        if roots.size == 0:
            return np.zeros(0)
        u = self.NTE_u
        v = self.NTE_v
        lca = self._lowestCommonAncestor(roots, u, v)
        dist = np.asarray(self.dist)
        weights = self.Weights
        if scipy.sparse.issparse(weights):
            weights = weights.tocsr()
        edgeLengths = np.asarray(weights[u, v], dtype=np.float64).ravel()
        return dist[roots, u] + dist[roots, v] + edgeLengths - 2.0 * dist[roots, lca]

    def nonTreeEdges(self, p):
        """ Non-tree edges of the shortest path tree rooted at p, as a list
            of (v1, v2) with v1 < v2, IN THE SAME ORDER as
            shortestPathTree(p, justEdges=True).

            The original rebuilds, for every root, the set of all edges of
            the root's connected component and iterates over
            (that set) - (tree edges). Here we build the component edge set
            once per component, from the identical sequence of edges, and
            reuse it. A set built by the same insertions in the same order
            has the same internal layout, and the difference has the same
            sizes on both sides, so iterating it gives the same order.
        """
        if self.CCV is None:
            # Components unknown: fall back to the original method
            return self.shortestPathTree(p, justEdges=True)

        if self._edgeW is None:
            rw, cl, _ = scipy.sparse.find(self.Weights)
            self._edgeW = list(zip(list(rw), list(cl)))
            self._compEdgeSet = {}

        # Edge set of p's component, built from the same filtered sequence
        # the original uses (edges whose first vertex is reachable from p)
        comp = self.CCV[p]
        if comp not in self._compEdgeSet:
            link = np.array(self.pred[p, :])
            link[p] = p
            filtered = [x for x in self._edgeW if link[x[0]] != -9999]
            self._compEdgeSet[comp] = set(filtered)
        edgeSet = self._compEdgeSet[comp]

        # Tree edges of the shortest path tree of p, in both orientations:
        # every reachable vertex other than p is joined to its predecessor
        link = self.pred[p, :]
        nodes = np.nonzero(link != -9999)[0]
        prevs = link[nodes]
        treeEdges = set(zip(prevs.tolist(), nodes.tolist()))
        treeEdges.update(zip(nodes.tolist(), prevs.tolist()))

        NTE = list(edgeSet - treeEdges)
        NTE = [x for x in NTE if (x[0] < x[1])]
        return NTE

    def fixNTE(self):
        """ Calcola un array di NTE per ogni sorgente p, e la salva in self.NTE
            Fast version: same lists, same order, via nonTreeEdges.
        """
        self.NTE = []
        for p in range(self.NVert):
            nte = self.nonTreeEdges(p)
            for e in nte:
                self.NTE.append((p, e))

    def buildCandidateArrays(self):
        """ Store the (root, non-tree edge) candidates as parallel arrays, in
            the order of self.NTE: root, the two endpoints, and the column of
            the edge in AnnMat (-1 if the edge is not a sentinel edge).
        """
        nCand = len(self.NTE)
        self.NTE_root = np.empty(nCand, dtype=np.int64)
        self.NTE_u = np.empty(nCand, dtype=np.int64)
        self.NTE_v = np.empty(nCand, dtype=np.int64)
        self.NTE_sent = np.empty(nCand, dtype=np.int64)
        for k, (p, e) in enumerate(self.NTE):
            u = int(e[0])
            v = int(e[1])
            self.NTE_root[k] = p
            self.NTE_u[k] = u
            self.NTE_v[k] = v
            self.NTE_sent[k] = self.SentinelIdx.get((u, v), -1)

    def resetLabels(self):
        """ Resetta i labels per poterli calcolare rispetto al nuovo vettore di supporto
        """
        self.Labels = -1 * np.ones((self.NVert, self.NVert), dtype=int)

    def cycleAnnotation(self, C):
        """ Calcola l'annotazione del ciclo C (somma delle ann degli edges)

            fast-v3: vectorised, same values; returns a (dimH1, 1) int array.
        """
        onCycle = np.asarray(C).ravel() > 0
        columns = self.edgeSentinelColumn[onCycle]
        columns = columns[columns >= 0]
        return (self.AnnMat[:, columns].sum(axis=1) % 2).reshape(-1, 1)


    def findShortestNonOrtho(self, i):
        """ Funzione che genera una lista di cicli non-ortogonali (m = 1) a Sup
            che è l' i-esimo vettore di supporto.
            I NTE seguono la convenzione lessicografica
            If allDraws ritorna una lista di tutti i cicli minimi a parità!

            The orthogonality test is evaluated for all candidates at once.
            A candidate (root s, non-tree edge (u, v)) is non-orthogonal to
            Sup iff Labels[s, u] + Labels[s, v] + bit(u, v) = 1 (mod 2), with
            bit(u, v) the sentinel bit of the edge (0 if not sentinel). The
            selected candidates keep the order of self.NTE, so ties are
            resolved exactly as before.

            fast-v3: candidate lengths are computed once per graph (initSup),
            and cycle vectors are built only for the shortest candidates.
            With draws, identical tied cycles are compared once: a copy of a
            kept cycle is never kept again, and a copy of a rejected one is
            rejected for the same reason, so the result is unchanged.
        """
        Sup = self.S[:, i-1]
        _ = self.computeLabels(Sup)

        # Sentinel bits, with an extra 0 at the end for index -1
        bits = self.sentinelBits(Sup)
        bitsExt = np.append(bits, 0)
        edgeBit = bitsExt[self.NTE_sent]

        labelU = self.Labels[self.NTE_root, self.NTE_u]
        labelV = self.Labels[self.NTE_root, self.NTE_v]
        orthoTest = (labelU + labelV + edgeBit) % 2
        selected = np.nonzero(orthoTest == 1)[0]

        if getattr(self, 'candidateLength', None) is None:
            self.candidateLength = self._candidateLengths()
        lengths = self.candidateLength[selected]

        # MULTIPROCESSING!
        if self.pipe is not None:
            self.pipe.send([self.eps, mlt.current_process().pid, i, self.dimH1, self.NonZ])

        # The tree-distance lengths add the same weights as
        # WEdges . (cycle vector), in another order, so they may differ in
        # the last bits for non-integer weights. We therefore shortlist the
        # candidates within a tiny tolerance of the minimum, build their
        # cycles, and choose among them with the original lengths and the
        # original rules: candidates outside the shortlist are longer anyway.
        smallest = lengths.min()
        tolerance = 1e-9 * max(abs(smallest), np.finfo(np.float64).tiny)
        shortlist = selected[lengths <= smallest + tolerance]
        cycles = [self.closeCycle(self.NTE[k][1], force=int(self.NTE_root[k])) for k in shortlist]
        exactLengths = [np.dot(self.WEdges, c) for c in cycles]
        shortest = min(exactLengths)
        tied = [c for c, l in zip(cycles, exactLengths) if l == shortest]

        if not self.allDraws:  # trovane solo uno: il primo dei più corti
            return np.array(tied[0], dtype=int)

        # tutti i pareggi omologhi al primo, distinti
        first = tied[0]
        targetAnnotation = self.cycleAnnotation(first)
        kept = [first]
        seen = {np.asarray(first).tobytes()}
        for cycle in tied[1:]:
            key = np.asarray(cycle).tobytes()
            if key in seen:
                continue
            seen.add(key)
            if np.array_equal(self.cycleAnnotation(cycle), targetAnnotation):
                kept.append(cycle)
        return [np.array(c, dtype=int) for c in kept]

    def updateSup(self, newC, index):
        """ Funzione che aggiorna i support vector per ogni nuovo ciclo minimo
            della base di omologia. Fa in pratica un Grahm-Schmidt in Z2
            Deve calcolare l'annotazione del ciclo.
            index deve andare da 0 a g-1
        """
        sup_i = self.Support[index]  # support vector del nuovo ciclo
        AnnNewC = self.cycleAnnotation(newC)  # annotazione del nuovo ciclo

        for j in range(index+1, self.dimH1):  # per S_j da i+1 a g
            sup_j = self.Support[j]
            self.Support[j] = (sup_j + sup_i * self.innerProd(AnnNewC, sup_j)) % 2  # somma ( =sottrae) la proiezione

    def initSup(self, Draws=True):
        """ Function to initialize the Support vectors structures self.Support and
            self.S. Plus it preallocates L,C and B. Then fixes the non-tree edges.
            Also, initializes the Gemoetry.Basis for self.SHB (Edges are computed)
            REQUIRES THE ANNOTATION TO HAVE BEEN COMPUTED!
        """
        self.S = np.eye(self.dimH1, dtype=np.int64)
        for i in range(self.dimH1):
            self.Support.append(self.S[:, i])

        # init the matrices
        self.L = np.zeros((self.dimH1, len(self.Sentinel)), dtype=np.int64)
        self.C = np.zeros((len(self.Sentinel), self.dimH1), dtype=np.int64)
        self.B = np.zeros((self.dimH1, self.dimH1), dtype=np.int64)

        # compute labels as a starter: row i of L holds the sentinel bits of
        # support vector S[:, i] (same values as the original double loop)
        for i in range(self.dimH1):
            self.L[i, :] = self.sentinelBits(self.S[:, i])

        # fix non-tree edges as we need them for the base case
        _ = self.fixNTE()
        self.buildCandidateArrays()

        # FIND A FILTER BETWEEN EDGES AND SENTINEL EDGES
        Filter = np.array([1 if ((int(x[0]), int(x[1])) in self.SentinelSet) else 0 for x in self.Edges])
        self.SentFilter = np.nonzero(Filter)

        self.SHB = GeometryN.Basis(self.NEdges, edgelist=self.Edges)

        self.allDraws = Draws

        # fast-v3: tree structures and candidate lengths, once per graph
        self._prepareTrees()
        self.candidateLength = self._candidateLengths()

    def mulz2(self, a, b):
        """ Product over Z2. In float64, so that NumPy uses BLAS: the entries
            are sums of 0/1 products, integers far below 2^53, so exact.
        """
        product = np.matmul(np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64))
        return np.mod(product, 2)

    def sumz2(self, a, b):
        return np.mod(a+b, 2).astype(np.float64)

    def inv_z2(self, a):
        """ Inverse over Z2, exact (Gauss-Jordan). The previous version,
            mod(round(inv(a) * det(a)), 2), is exact only while det(a) and
            the adjugate stay below 2^53: it failed from about 50 x 50.
        """
        return gf2_inverse(a)

    def UpdateMehlhorn(self, i, k):
        """ Update function optimized for the divide et impera recursive approach
            as per Mehlhorn (2004), Dey (2018).
            Uses instance variables L,C,B and S. Receives recursion parameters
            i and k.
        """
        lam = i + floor(k/2)-1
        mu = i+k-1

        # slice matrix B to obtain X and Y
        X = self.B[i - 1:lam, i-1:lam]
        Y = self.B[lam+1-1:mu, i-1:lam]

        # By definition of A
        A = self.mulz2(Y, self.inv_z2(X))

        # Update Support vectors
        S1 = np.zeros((self.dimH1, mu - lam + 1), dtype=int)
        S1 = self.sumz2(self.S[:, lam+1-1:mu], self.mulz2(self.S[:, i-1:lam], A.T))
        self.S[:, lam+1-1:mu] = S1

        # Mirror in the self.Sup structure
        for col in range(lam, mu):
            self.Support[col] = self.S[:, col]

        # Update Labels
        # Generate W and Z
        m, _ = A.shape
        W = np.concatenate((A, np.eye(m)), axis=1)

        Z = self.L[i-1:mu, :]

        L1 = self.mulz2(W, Z)
        self.L[lam+1-1:mu, :] = L1

        # At last, update matrix B = L C. fast-v3: only the rows lam..mu-1 of
        # L have changed, and C has not changed in this update, so only those
        # rows of B change; recomputing them keeps B equal to L C exactly.
        self.B[lam+1-1:mu, :] = self.mulz2(self.L[lam+1-1:mu, :], self.C)

    def projOnSent(self, C):
        """
            Compute the Sentinel-Edge trace of cycle C. Return a vector that is
            as long as self.Sentinel.
            Uses the structure self.SentFilter which projects a cycle onto the
            space of sentinel edges.
            In case it is a list of draws, pick the first!
        """
        if type(C) is list:  # in case it is a list of draws, pick the first!
            C = C[0]

        C = np.array(C, dtype=int)
        return C[self.SentFilter]

    def ExtendBasis(self, i, k):
        """ Recursive procedure to compute the minimal homology basis. Follows
            Mehlhorn(2004) and Dey (2018).
        """
        if k == 0:  # If homology is trivial it returns the empty basis
            return
        if k == 1:  # BASE CASE. Will search for the shortest non-orthogonal cycle
            newCycle = self.findShortestNonOrtho(i)

            # save the new cycle to the SHB
            self.SHB.addEl(newCycle)
            # compute the Sentinel-edge trace of the new cycle
            self.C[:, i-1] = self.projOnSent(newCycle)
            self.B[:, i-1] = self.mulz2(self.L, self.C[:, i-1])
        else:
            self.ExtendBasis(i, floor(k/2.0))  # extend by floor(k/2) elements
            self.UpdateMehlhorn(i, k)  # updates support vectors
            self.ExtendBasis(i + floor(k/2.0), ceil(k/2.0))

    # THE END
