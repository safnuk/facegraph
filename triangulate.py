#triangulate.py
import math
import numpy as np
import scipy.optimize as opt
#import enthought.mayavi.mlab as mlab

class Vertex:
    """ Class which encodes incidence data at a vertex.

    	Incidence data should be encoded as index data to the positions in lists of
    	vertices, edges and triangles.

    	Note: this class does not record vertex position.
    """
    def __init__(self):
        self.incidentVertices = []
        self.incidentTriangles = []
        self.incidentEdges = []

    def __repr__(self):
        return ('<Vertex: V' + str(self.incidentVertices) + '; E' + 
                str(self.incidentEdges) + '; T' + str(self.incidentTriangles) + '>')

    def __str__(self):
        return ('<Vertex: V' + str(self.incidentVertices) + '; E' + 
                str(self.incidentEdges) + '; T' + str(self.incidentTriangles) + '>')

class Edge:
    """ Edge class for encoding connectivity structure. Vertices should be a 2-element
    	list of indices to the vertex list  (endpoints of the edge). 

    	Incident triangles are recorded as indices to the list of triangles.
    """
    def __init__(self, vertices):
        self.vertices = vertices
        self.vertices.sort()
        self.incidentTriangles = []

    def __repr__(self):
        return '<Edge: V' + str(self.vertices) + '; T' + str(self.incidentTriangles) + '>'

    def __str__(self):		
        return '<Edge: V' + str(self.vertices) + '; T' + str(self.incidentTriangles) + '>'

    def getOtherVertex(self, v):
        i = self.vertices.index(v)
        if i == 0:
            return self.vertices[1]
        else:
            return self.vertices[0]

class Triangle:
    """ Triangle class used to encode connectivty structure. Vertices should be a 
    	counter-clockwise list of indices. Edge indices are also listed counterclockwise,
    	with the edge being opposite to the vertex.

    	i.e. if a triangle has vertices a, b, c with opposite edges A, B, C, then vertices 
    	listed in the order (a, b, c) measn edges should be listed in the order (A, B, C).
    """
    def __init__(self, vertices):
        self.vertices = vertices
        self.edges = [-1, -1, -1]

    def __repr__(self):
        return '<Triangle: V' + str(self.vertices) + '; E' + str(self.edges) + '>'

    def __str__(self):		
        return '<Triangle: V' + str(self.vertices) + '; E' + str(self.edges) + '>'



class Triangulation:
    """ Class which encodes the incidence relations of vertices, edges and faces 
    	of a triangulation initialized with a list of points [[x0, y0, z0], ...] 
    	and a list of triples of indices [[v0, v1, v2], ...] in counterclockwise order. 

    	The indices for the points are implied from the given ordering.  

    	This class also can replace the induced edge-length metric with a circle-
    	packing metric which minimizes the least-squares difference between edge lengths.
    """
    def __init__(self, points, faces, radii=None, edgeInfo=None):

        self.faces = np.array(faces)		# data structure used to display surface with mlab
        self.x = np.empty(len(points))
        self.y = np.empty(len(points))
        self.z = np.empty(len(points))
        self.vertices = []
        self.edges = []
        self.triangles = []

        for index, point in enumerate(points):
            self.vertices.append(Vertex())
            self.x[index] = point[0]
            self.y[index] = point[1]
            self.z[index] = point[2]

        self.edgeDict = {}	# dictionary used to map (v1,v2) pairs to edge index (in edge list)
        edgeIndex = 0
        for triangleIndex, face in enumerate(faces):
            vertexList = [face[0], face[1], face[2]]
            triangle = Triangle(vertexList)
            self.triangles.append(triangle)
            for vertex in vertexList: 			# record triangle incidence in vertices
                self.vertices[vertex].incidentTriangles.append(triangleIndex)

            edgeList = 	[(face[1], face[2]), (face[2], face[0]), (face[0], face[1])]
            for i, e in enumerate(edgeList):
                if e[0] > e[1]:
                    edgeList[i] = (e[1], e[0])

            for i, vertexPair in enumerate(edgeList):
                # if edge not already created, make new one
                if vertexPair not in self.edgeDict:
                    self.edgeDict[vertexPair] = edgeIndex
                    edgeIndex += 1
                    edge = Edge(list(vertexPair))
                    self.edges.append(edge)
                    currentEdgeIndex = edgeIndex - 1
                else:
                    currentEdgeIndex = self.edgeDict[vertexPair]
                    edge = self.edges[currentEdgeIndex]
                # add incidence data to vertices, edge, triangles
                edge.incidentTriangles.append(triangleIndex)
                for vertex in vertexPair:
                    self.vertices[vertex].incidentEdges.append(currentEdgeIndex)
                triangle.edges[i] = currentEdgeIndex
            for v in face:
                self.vertices[v].incidentTriangles.append(triangleIndex)
                for w in face:
                    if w <> v:
                        self.vertices[v].incidentVertices.append(w)
        self.uniquifyIncidenceData()   # remove duplicate edge and triangle references from vertices
        if radii is not None:
            self.radii = np.array(radii)
        if edgeInfo is not None:
            self.circleAngles = np.zeros(len(self.edges))
            self.edgeLengths = np.zeros(len(self.edges))
            for edgeKey, info in edgeInfo.items():
                eIndex = self.edgeDict[edgeKey]
                self.circleAngles[eIndex] = info[0]
                self.edgeLengths[eIndex] = info[1]
        if (radii is None) and (edgeInfo is None):   # should only need to be performed once
            self.improveTriangulation(verbose=True)
        else:
            self.findBoundaries() # normally done when improving triangulation
        self.calcBoundaryCycles()
        self.innerAngles = np.r_[[[0.0, 0.0, 0.0]]*len(self.triangles)]

    def uniquifyIncidenceData(self):
        """ Uniquify incidence lists. """
        for vertex in self.vertices:
            vertex.incidentEdges = uniquify(vertex.incidentEdges)
            vertex.incidentVertices = uniquify(vertex.incidentVertices)
            vertex.incidentTriangles = uniquify(vertex.incidentTriangles)

    def removeIsolatedVertices(self, verbose):
        """ Looks for unreferenced vertices and removes them from the list. Updates all
        	vertex indices to account for the shift.

        	These vertices cause problems with Ricci flow.

        	Note: this function does not modify self.faces, so must call 
        	self.recalcFaces() afterwords.	
        """
        vertexDict = {}
        isolatedVertices = []
        # make dictionary relating old vertex indices to new
        for index, vertex in enumerate(self.vertices):
            if len(vertex.incidentVertices) == 0:  #isolated vertex condition
                isolatedVertices.append(index)
            else:
                vertexDict[index] = index - len(isolatedVertices)
        # remove isolated vertices
        for position, index in enumerate(isolatedVertices):
            self.vertices.pop(index - position)
        self.x = np.delete(self.x, isolatedVertices)
        self.y = np.delete(self.y, isolatedVertices)
        self.z = np.delete(self.z, isolatedVertices)

        # adjust indices of all vertex, edge, triangle data
        for vertex in self.vertices:
            vertex.incidentVertices = [vertexDict[x] for x in vertex.incidentVertices]
        for edge in self.edges:
            edge.vertices[0] = vertexDict[edge.vertices[0]]
            edge.vertices[1] = vertexDict[edge.vertices[1]]
        for triangle in self.triangles:
            for j in range(3):
                triangle.vertices[j] = vertexDict[triangle.vertices[j]]
        if verbose:
            print "Removed " + str(len(isolatedVertices)) + " isolated vertices."

    def improveTriangulation(self, verbose=False):
        self.removeIsolatedVertices(verbose)	
        self.findBoundaries()
        self.adjustBivalentVertices(verbose)
        self.adjustTrivalentVertices(verbose)
        self.adjustBoundaryTrivalentVertices(verbose)
        self.adjustQuadvalentVertices(verbose)

        self.recalcFaces()

    def adjustBivalentVertices(self, verbose):
        counter = 0
        for v, vertex in enumerate(self.vertices):
            if len(vertex.incidentVertices) == 2:		
                counter += 1
                #if verbose:
                #	print "Adjusting bivalent vertex ", v, vertex
                v1, v2, v3, v4, e1, e2, e3, e4 = self.__findNeighborsBivalent(v, vertex)
                triangle1 = Triangle([v, v3, v1])
                triangle2 = Triangle([v, v2, v4])
                edge5 = Edge([v, v3])
                edge6 = Edge([v, v4])
                t1 = len(self.triangles)
                t2 = t1 + 1
                e5 = len(self.edges)
                e6 = e5 + 1

                triangle1.edges = [e3, e1, e5]
                triangle2.edges = [e4, e6, e2]
                self.edges.extend([edge5, edge6])
                self.triangles.extend([triangle1, triangle2])
                self.boundaryEdges = np.concatenate((self.boundaryEdges, np.array([1,1])))
                for x in [e1, e2, e3, e4]:
                    self.boundaryEdges[x] = 0
                for x in [v1, v2]:
                    self.boundaryVertices[x] = 0

                for x in [e1, e3, e5]:
                    self.edges[x].incidentTriangles.append(t1)
                for x in [e2, e4, e6]:
                    self.edges[x].incidentTriangles.append(t2)
                for x in [v, v1, v3]:
                    self.vertices[x].incidentTriangles.append(t1)
                for x in [v, v2, v4]:
                    self.vertices[x].incidentTriangles.append(t2)
                for x in [v3, v4]:
                    self.vertices[x].incidentVertices.append(v)
                for x in [v, v3]:
                    self.vertices[x]. incidentEdges.append(e5)
                for x in [v, v4]:
                    self.vertices[x].incidentEdges.append(e6)
                self.vertices[v].incidentVertices.extend([v3, v4])
        if verbose:
            print "Improved " + str(counter) + " bivalent vertices."

    def adjustTrivalentVertices(self, verbose):
        counter = 0
        for v, vertex in enumerate(self.vertices):
            if len(vertex.incidentVertices) == 3 and self.boundaryVertices[v] == 0:		
                #if verbose:
                #	print "Adjusting trivalent vertex ", v, vertex
                v1, v2, v3 = vertex.incidentVertices
                e1 = self.findEdge(v2, v3)
                e2 = self.findEdge(v1, v3)
                e3 = self.findEdge(v1, v2)
                ranking = [(self.rateTrivalentEdgeSwap(x), x) for x in [e1, e2, e3]]
                ranking.sort()
                ranking.reverse()
                l1, e1 = ranking[0]
                if l1 > 0:
                    self.edgeSwap(e1)
                    counter += 1
                elif verbose:
                    print "Could not fix trivalent vertex ", v, vertex
        if verbose:
            print "Improved " + str(counter) + " trivalent vertices."

    def adjustBoundaryTrivalentVertices(self, verbose):
        counter = 0
        badCount = 0
        for v, vertex in enumerate(self.vertices):
            if len(vertex.incidentVertices) == 3 and self.boundaryVertices[v] <> 0:		
                #if verbose:
                #	print "Adjusting boundary trivalent vertex ", v, vertex
                v1, v2, v3 = vertex.incidentVertices
                edges = [self.findEdge(v1, v2), self.findEdge(v1, v3),
                         self.findEdge(v2, v3)]
                edges.remove(-1)
                ranking = [(self.rateQuadvalentSwap(x), x) for x in edges]
                ranking.sort()
                ranking.reverse()
                l1, e1 = ranking[0]
                if l1 >= 1:
                    self.edgeSwap(e1)
                    counter += 1
                else:
                    badCount += 1
                    #print "Could not improve boundary trivalent vertex ", v, vertex
        if verbose:
            print "Improved " + str(counter) + " trivalent boundary vertices."
            print "Could not improve " + str(badCount) + " vertices."

    def adjustQuadvalentVertices(self, verbose):
        counter = 0
        badCount = 0
        for v, vertex in enumerate(self.vertices):
            if len(vertex.incidentVertices) == 4 and self.boundaryVertices[v] == 0:		
                #if verbose:
                #	print "Adjusting quadvalent vertex ", v, vertex
                incVertices= vertex.incidentVertices[:]
                edges = [self.findEdge(v1, v2) for v1 in incVertices for v2 in incVertices]
                edges = uniquify(edges)
                edges.remove(-1)
                ranking = [(self.rateQuadvalentSwap(x), x) for x in edges]
                ranking.sort()
                ranking.reverse()
                l1, e1 = ranking[0]
                if l1 >= 1:
                    self.edgeSwap(e1)
                    counter += 1
                else:
                    badCount += 1
                    #print "Could not improve boundary trivalent vertex ", v, vertex

        if verbose:
            print "Improved " + str(counter) + " quadvalent vertices."
            print "Could not improve " + str(badCount) + " vertices."

    def findBadVertices(self):
        badList = []
        for v, vertex in enumerate(self.vertices):
            if ((len(vertex.incidentVertices) == 4 and self.boundaryVertices[v] == 0)
                or (len(vertex.incidentVertices) == 3 and self.boundaryVertices[v] <> 0)):		
                incVertices= vertex.incidentVertices[:]
                edges = [self.findEdge(v1, v2) for v1 in incVertices for v2 in incVertices]
                edges = uniquify(edges)
                edges.remove(-1)
                angles = [self.circleAngles[x] for x in edges]
                angles.sort()
                angles.reverse()
                if angles[0] < 1.570797:
                    badList.append(v)
        return badList


    def __findNeighborsBivalent(self, v, vertex):		
        t = vertex.incidentTriangles[0]
        triangle = self.triangles[t]
        i = triangle.vertices.index(v)
        if i == 0:
            e1, e2 = triangle.edges[2], triangle.edges[1]
        elif i == 1:
            e1, e2 = triangle.edges[0], triangle.edges[2]		
        elif i == 2:
            e1, e2 = triangle.edges[1], triangle.edges[0]

        v1 = self.edges[e1].getOtherVertex(v)
        v2 = self.edges[e2].getOtherVertex(v)

        L = [e for e in self.vertices[v1].incidentEdges if self.boundaryEdges[e] <> 0]
        L.remove(e1)
        e3 = L[0]
        L = [e for e in self.vertices[v2].incidentEdges if self.boundaryEdges[e] <> 0]
        L.remove(e2)
        e4 = L[0]

        v3 = self.edges[e3].getOtherVertex(v1)
        v4 = self.edges[e4].getOtherVertex(v2)

        return v1, v2, v3, v4, e1, e2, e3, e4

    def rateTrivalentEdgeSwap(self, e):
        """ Rates the improvement to the triangulation if edge e is swapped.

        	Assumes that swapping e will make it incident to a non-boundary
        	trivalent vertex.
        	Rating:
        		0	e is a boundary edge
        		0	swapping e will create a new trivalent or 
        			boundary 2-valent vertex
        		1	swapping e will create a new quad-valent or boundary 3-valent
        			vertex
        		2	swapping e will not make any bad vertices
        """
        if self.boundaryEdges[e] <> 0:
            return 0
        v1, v2 = self.edges[e].vertices
        valence1 = len(self.vertices[v1].incidentEdges)
        valence2 = len(self.vertices[v2].incidentEdges)
        if self.boundaryVertices[v1] <> 0:
            valence1 += 1
        if self.boundaryVertices[v2] <> 0:
            valence1 += 2
        if min([valence1, valence2]) <= 4:
            return 0
        if min([valence1, valence2]) == 5:
            return 1
        else:
            return 2

    def rateQuadvalentSwap(self, e):
        """ Rates the improvement to the triangulation if edge e is swapped.

        	Assumes that swapping e will make it incident to a boundary
        	trivalent vertex or non-boundary 4-valent vertex.
        	Rating:
        		0	e is a boundary edge
        		0 	swapping e will create a new trivalent or boundary
        			2-valent vertex
        		0	swapping e will replace it with a shorter edge
        		2	swapping e will not make any bad vertices
        		1 + originalLength / New Length if swapping creates
        			a new 4-valent or boundary 3-valent vertex, with new
        			longer edge 
        			(longer edges are preferred so that circle angle > pi/2)
        """
        if self.boundaryEdges[e] <> 0:
            return 0
        vList = self.edges[e].vertices
        valences = [len(self.vertices[x].incidentEdges) for x in vList]
        for i in range(2):
            if self.boundaryVertices[vList[i]] <> 0:
                valences[i] += 1
        M = min(valences)
        if M <= 4:
            return 0
        if M == 5:
            oppositeVertices = self.findOppositeVertices(e)
            originalEdgeLength = self.calcPointDistance(*vList)
            newEdgeLength = self.calcPointDistance(*oppositeVertices)
            if newEdgeLength <= originalEdgeLength:
                return 0
            return (1 + (originalEdgeLength / newEdgeLength))
        else:
            return 2

    def findEdge(self, v1, v2):
        """ Finds the index of the edge joining v1 and v2.

        	Assumes that vertex labels are ordered for each edge.
        	i.e. edge.vertices[0] < edge.vertices[1].

        	Returns -1 if no edge found
        """
        if v1 < v2:
            w1, w2 = v1, v2
        elif v1 > v2:
            w1, w2 = v2, v1
        else:
            return -1
        vertex = self.vertices[w1]
        for e in vertex.incidentEdges:
            edge = self.edges[e]
            if edge.vertices[1] == w2:
                return e
        # failed - no match found
        return -1

    def findOppositeVertices(self, edgeIndex):
        """ Given a non-boundary edge, finds the vertices opposite the edge
        	in each triangle incident to it.
        """
        edge = self.edges[edgeIndex]
        if len(edge.incidentTriangles) < 2:
            raise BoundaryEdgeError("Cannot find opposing vertices to a boundary edge.")
        vList = []
        for t in edge.incidentTriangles:
            triangle = self.triangles[t]
            index = triangle.edges.index(edgeIndex)
            vList.append(triangle.vertices[index])
        return vList

    def removePureBoundaryTriangles(self):
        """ Removes triangles with all three vertices on the boundary. These are not allowed 
        	for Ricci flow to converge.

        	Triangles are recognized by having a vertex incident to a single triangle. In this case,
        	a new vertex is added to the middle, and three new triangles take its place.

        	This function recursively calls itself until no more pure boundary triangles remain.
        	However, there is no depth limit to the recursion.

        	Note: this function does not modify self.faces, so must call 
        	self.recalcFaces() afterwords.
        """
        badTriangles = []
        for vertex in self.vertices:
            if len(vertex.incidentTriangles) == 1:
                badTriangles.extend(vertex.incidentTriangles)
        badTriangles.sort()
        # Note: problems come up here if list is not uniquefied, but this indicates
        # bad triangulation data. Currently, nothing is done to handle this case.

        if len(badTriangles) > 0:
            # Construct dictionary relating old triangle indices to new
            badTriangles.append(len(self.triangles))
            pos = 0
            triangleDict = {}
            isolatedEdges = []
            for index, triangle in enumerate(self.triangles):
                if index < badTriangles[pos]:
                    triangleDict[index] = index - pos	
                elif index == badTriangles[pos]:
                    pos += 1
                    triangle = self.triangles[index]
                    # remove references to bad triangle
                    for vIndex in triangle.vertices:
                        self.vertices[vIndex].incidentTriangles.remove(index)
                    for eIndex in triangle.edges:
                        self.edges[eIndex].incidentTriangles.remove(index)
            badTriangles.remove(len(self.triangles))
            # remove offending triangles
            for position, index in enumerate(badTriangles):
                self.triangles.pop(index - position)
            # adjust triangle indices of vertices and edges
            for vIndex, vertex in enumerate(self.vertices):
                vertex.incidentTriangles = [triangleDict[x] for x in vertex.incidentTriangles]
                if len(vertex.incidentTriangles) == 0:  # now an isolated vertex
                    for j in vertex.incidentVertices:
                        v2 = self.vertices[j]
                        v2.incidentVertices.remove(vIndex)
                    vertex.incidentVertices = []
                    isolatedEdges.extend(vertex.incidentEdges)
            self.__removeIsolatedEdges(isolatedEdges)
            for edge in self.edges:
                edge.incidentTriangles = [triangleDict[x] for x in edge.incidentTriangles]
            self.removePureBoundaryTriangles()  #loop again to remove any newly created bad triangles

    def __removeIsolatedEdges(self, isolatedEdges):
        if len(isolatedEdges) > 0:
            isolatedEdges.sort()
            isolatedEdges.append(len(self.edges))
            pos = 0
            edgeDict = {}
            for index, edge in enumerate(self.edges):
                if index < isolatedEdges[pos]:
                    edgeDict[index] = index - pos
                elif index == isolatedEdges[pos]:
                    pos += 1
                    edge = self.edges[index]
                    for vIndex in edge.vertices:
                        self.vertices[vIndex].incidentEdges.remove(index)

            isolatedEdges.remove(len(self.edges))
            for position, index in enumerate(isolatedEdges):
                self.edges.pop(index - position)
            for vertex in self.vertices:
                vertex.incidentEdges = [edgeDict[x] for x in vertex.incidentEdges]
            for triangle in self.triangles:
                triangle.edges = [edgeDict[x] for x in triangle.edges]


    def findBoundaries(self):
        """ Identifies all edges on the boundary, and identifies which component it lies on
        	Also lists the boundary data for the incident vertices

        	Currently, it does not distinguish boundary components - possible values
        	are 0 (not on boundary) or 1 (on boundary)
        """
        self.boundaryEdges = np.zeros(len(self.edges), dtype='int32')
        self.boundaryVertices = np.zeros(len(self.vertices), dtype='int32')
        for index, edge in enumerate(self.edges):
            if len(edge.incidentTriangles) == 1:  # we've found a boundary
                self.boundaryEdges[index] = 1
                v1, v2 = edge.vertices
                self.boundaryVertices[v1] = 1
                self.boundaryVertices[v2] = 1

    def calcBoundaryCycles(self):
        """ Calculates the oriented cycle of edges running around each boundary component.
        	Compiles a list of boundary cycle, each of which is a list of edges running in 
        	oriented direction around boundary (assuming counterclockwise orientation
        	of triangle vertices, surface will be on left side)
        """
        self.boundaryCycles = []
        bEdges = [int(x) for x in self.boundaryEdges]
        firstEdge = self.getNextBoundaryComponent(bEdges)
        while firstEdge <> -1:
            bEdges[firstEdge] = 0		# mark off edges that have been hit
            boundaryCycle = [firstEdge]
            nextEdge = self.getNextBoundaryEdge(firstEdge)
            while nextEdge <> firstEdge:
                boundaryCycle.append(nextEdge)
                bEdges[nextEdge] = 0
                nextEdge = self.getNextBoundaryEdge(nextEdge)
            self.boundaryCycles.append(boundaryCycle[:])
            firstEdge = self.getNextBoundaryComponent(bEdges)

    def getNextBoundaryComponent(self, bEdges):
        """ Returns the first nonzero entry in list boundaryEdges. If all entries are zero,
        	returns -1.
        """
        for e, status in enumerate(bEdges):
            if status <> 0:
                return e
        return -1

    def getNextBoundaryEdge(self, edgeIndex):
        """ Returns the edge adjacent to given edge along the boundary (using induced orientation).

        	If passed edge is not a boundary, an exception is raised.
        """
        if self.boundaryEdges[edgeIndex] == 0:
            message = "Edge " + str(edgeIndex) + " is not a boundary edge."
            raise BoundaryEdgeError(message)
        # get ordered pair of vertices adjacent to edge
        tIndex = self.edges[edgeIndex].incidentTriangles[0]
        triangle = self.triangles[tIndex]
        i = triangle.edges.index(edgeIndex)
        if i == 0:
            v1, v2 = triangle.vertices[1:3]
        elif i == 2:
            v1, v2 = triangle.vertices[0:2]
        else:
            v1, v2 = triangle.vertices[2], triangle.vertices[0]
        edgeList = [x for x in self.vertices[v2].incidentEdges if self.boundaryEdges[x]]
        edgeList.remove(edgeIndex)
        return edgeList[0]

    def edgeSwap(self, edgeIndex):
        """ Replaces an edge with its complementary edge (Whitehead move)

        	<|>     ---->     <->		
        """
        edge = self.edges[edgeIndex]
        if len(edge.incidentTriangles) < 2:
            raise BoundaryEdgeError("Cannot swap a boundary edge")
        v1, v2 = edge.vertices
        # t1 on right of edge when going v1 -> v2
        t1, t2 = self.determineTriangleSides(edge)

        # gather relevant edges and vertices
        # triangle1 and 2 can be modified because they will be replaced
        triangle1 = self.triangles[t1]
        triangle2 = self.triangles[t2]
        triangle1.vertices.remove(v1)
        triangle1.vertices.remove(v2)
        triangle2.vertices.remove(v1)
        triangle2.vertices.remove(v2)
        i1 = triangle1.edges.index(edgeIndex)
        i2 = triangle2.edges.index(edgeIndex)
        triangle1.edges.remove(edgeIndex)
        triangle2.edges.remove(edgeIndex)
        v3 = triangle1.vertices[0]
        v4 = triangle2.vertices[0]
        if i1 == 1:
            e2, e1 = triangle1.edges
        else:
            e1, e2 = triangle1.edges
        if i2 == 1:
            e3, e4 = triangle2.edges
        else:
            e4, e3 = triangle2.edges

        # adjust incidence data
        vertex1, vertex2, vertex3, vertex4 = [self.vertices[x] for x in [v1, v2, v3, v4]]
        edge1, edge2, edge3, edge4 = [self.edges[x] for x in [e1, e2, e3, e4]]
        triangle1.vertices = [v1, v3, v4]
        triangle1.edges = [edgeIndex, e3, e1]
        triangle2.vertices = [v2, v4, v3]
        triangle2.edges = [edgeIndex, e2, e4]
        if v3 < v4:
            edge.vertices = [v3, v4]
        else:
            edge.vertices = [v4, v3]
        edge2.incidentTriangles.remove(t1)
        edge2.incidentTriangles.append(t2)
        edge3.incidentTriangles.remove(t2)
        edge3.incidentTriangles.append(t1)

        vertex1.incidentTriangles.remove(t2)
        vertex1.incidentEdges.remove(edgeIndex)
        vertex1.incidentVertices.remove(v2)
        vertex2.incidentTriangles.remove(t1)
        vertex2.incidentEdges.remove(edgeIndex)
        vertex2.incidentVertices.remove(v1)
        vertex3.incidentTriangles.append(t2)
        vertex3.incidentEdges.append(edgeIndex)
        vertex3.incidentVertices.append(v4)
        vertex4.incidentTriangles.append(t1)
        vertex4.incidentEdges.append(edgeIndex)
        vertex4.incidentVertices.append(v3)

    def determineTriangleSides(self, edge):
        v1, v2 = edge.vertices
        t1, t2 = edge.incidentTriangles
        triangle1 = self.triangles[t1]
        i1 = triangle1.vertices.index(v1)
        i2 = triangle1.vertices.index(v2)

        if (i1, i2) in [(0,2), (2,1), (1,0)]:
            return t1, t2
        else:
            return t2, t1

    def checkTriangleCycleCondition(self):
        """ Looks for vertices in the middle of a triangle, and removes them.

        	Saves having to first calculate circle angles and then check cocycle condition
        	necessary for convergence of Ricci flow.
        """
        badVertices = []
        for index, vertex in enumerate(self.vertices):
            if (len(vertex.incidentVertices) == 3) and (self.boundaryVertices[index] == 0):
                badVertices.append[index]
        if len(badVertices) > 0:
            pass
            #self.checkCycleCondition()

    def findOrientedEdgeVertices(self, edgeIndex, triangleIndex):
        """ Given a triangle containing referenced edge, returns the ordered list
        	of vertices v0, v1, v2 making up the triangle with v0 --> v1 being the
        	given edge.
        """
        triangle = self.triangles[triangleIndex]
        i = triangle.edges.index(edgeIndex)
        if i == 0:
            return [triangle.vertices[1], triangle.vertices[2], triangle.vertices[0]]
        if i == 1:
            return [triangle.vertices[2], triangle.vertices[0], triangle.vertices[1]]
        if i == 2:
            return [triangle.vertices[0], triangle.vertices[1], triangle.vertices[2]]

    def calcCurvatures(self):
        """ Calculates curvatures at all vertices.

        	K = 2pi - inner angle sum at vertex, unless vertex is on boundary, then
        	K = 2pi - 2(inner angle sum).
        	Assumes that self.innerAngles is calculated (eg by calcInnerAngles())
        """
        self.curvatures = np.r_[[2*math.pi] * len(self.vertices)]
        for index, vertex in enumerate(self.vertices):
            for triangle in vertex.incidentTriangles:
                #identify vertex in triangle 
                vertexOrder = self.triangles[triangle].vertices.index(index)
                self.curvatures[index] -= self.innerAngles[triangle, vertexOrder]
                if self.boundaryVertices[index]:
                    self.curvatures[index] -= self.innerAngles[triangle, vertexOrder]
        """
		# deal with isolated vertices (very cludgey fix!)
		for index, K in enumerate(self.curvatures):
			if K == 2*math.pi:
				self.curvatures[index] = 0
		"""

    def calcOptimalRadii(self):
        """ Attempts to recreate edge lengths of triangulation using a circle packing metric.

        The radius at each vertex is bound between 0 and infinity, while the angles of 
        intersection of the circles is between pi/2 and pi (i.e. the triangle determining
        edge length lActual from radii r1 and r2 has obtuse angle at angle opposite to lActual)
        """
        self.calcInitialEdgeLengths()
        self.calcInitialRadii()
        bounds = [(0, None) for v in self.vertices]
        min = opt.fmin_l_bfgs_b(errorFunc, self.initialRadii, 
                                errorGrad, args=(self, ), bounds=bounds)
        #min = opt.fmin_tnc(errorFunc, self.initialRadii, errorGrad, args=(self, ),
        #	bounds=bounds)
        self.radii = min[0]
        print "Min error = " + str(errorFunc(self.radii, self))
        self.calcCircleAngles()


    def calcInitialEdgeLengths(self):
        """ Calculate actual edge lengths of the triangulation based on distance formula in
        	3-space.
        """
        numberOfEdges = len(self.edges)
        X = np.empty(numberOfEdges)
        Y = np.empty(numberOfEdges)
        Z = np.empty(numberOfEdges)
        x = self.x
        y = self.y
        z = self.z

        for i, edge in enumerate(self.edges):
            j, k = edge.vertices[0],  edge.vertices[1] 
            X[i] = x[j] - x[k]
            Y[i] = y[j] - y[k]
            Z[i] = z[j] - z[k]
        self.initialEdgeLengths = np.sqrt(X**2 + Y**2 + Z**2)

    def calcPointDistance(self, v1, v2):
        """ Calculate the Euclidean (actual) distance between the two points.
        """
        P1 = np.array([self.x[v1], self.y[v1], self.z[v1]])
        P2 = np.array([self.x[v2], self.y[v2], self.z[v2]])
        D = P1 - P2
        return np.linalg.norm(D)

    def calcInitialRadii(self, percent=.70):
        """ Sets initial vertex radii to be 70% of the average edge length

        Assumes caller has first initalize initalEdgeLengths, e.g. by calling
        calcInitialEdgeLengths()
        """		
        self.initialRadii = np.empty(len(self.vertices))
        for index, vertex in enumerate(self.vertices):
            radius = 0.0
            for edge in vertex.incidentEdges:
                radius += self.initialEdgeLengths[edge]
            l = len(vertex.incidentEdges)
            if l <> 0:
                radius = (radius * percent) / l
            self.initialRadii[index] = radius

    def calcCircleAngles(self):
        """ Calculates angles of intersection for circles in circle packing metric.

        	Needs to know self.initialEdgeLengths and self.radii, 
        	e.g. first call findOptimalRadii().
        	This routine will also calculate the induced edge lengths and initializes
        	self.edgeLengths.
        """
        self.circleAngles = np.zeros(len(self.edges))
        self.edgeLengths = np.zeros(len(self.edges))

        for index, edge in enumerate(self.edges):
            v1, v2 = edge.vertices
            r1, r2 = self.radii[v1], self.radii[v2]
            lBar = self.initialEdgeLengths[index]
            if lBar >= r1 + r2:	#tangent circles 
                self.edgeLengths[index] = r1 + r2
                self.circleAngles[index] = math.pi
            elif math.cosh(lBar) <= math.cosh(r1) * math.cosh(r2):  #maximum overlap allowed
                self.edgeLengths[index] = math.acosh(math.cosh(r1) * math.cosh(r2)) 
                self.circleAngles[index] = math.pi/2
            else:	# arbitrary obtuse angle  (no error with actual edge length)
                self.edgeLengths[index] = lBar
                angle = 0.0
                try:  # watch for divide by 0!
                    temp = ((math.cosh(r1) * math.cosh(r2) - math.cosh(lBar)) 
                            / (math.sinh(r1) * math.sinh(r2)))
                    angle = math.acos(temp)
                except:
                    angle = math.pi
                self.circleAngles[index] = angle

    def calcInnerAngles(self):
        """ Calculate inner angles of all triangles.
        	Assumes that self.edgeLengths are calculated
        """
        coshEdgeLengths = np.cosh(self.edgeLengths)
        sinhEdgeLengths = np.sinh(self.edgeLengths)
        for index, triangle in enumerate(self.triangles):
            lengths = np.array([self.edgeLengths[edge] for edge in triangle.edges])
            coshLengths = np.array([coshEdgeLengths[edge] for edge in triangle.edges])
            sinhLengths = np.array([sinhEdgeLengths[edge] for edge in triangle.edges])
            coshOffsets = np.array([coshLengths[1]*coshLengths[2], 
                                    coshLengths[0] * coshLengths[2],
                                    coshLengths[0] * coshLengths[1]])
            sinhOffsets = np.array([sinhLengths[1]*sinhLengths[2], 
                                    sinhLengths[0] * sinhLengths[2],
                                    sinhLengths[0] * sinhLengths[1]])
            self.innerAngles[index] = np.arccos((coshOffsets - coshLengths) / sinhOffsets )


    def calcEdgeLengths(self):
        """ Calculates edge lengths using circle packing metric.

        	Function assumes that self.radii has been calculated (e.g. by ricciflow)
        	and also circleAngles (e.g. by calling calcCircleAngles()). 
        """
        coshRadii = np.cosh(self.radii)
        sinhRadii = np.sinh(self.radii)
        cosAngles = np.cos(self.circleAngles)
        for index, edge in enumerate(self.edges):
            v1, v2 = edge.vertices
            length = math.acosh(coshRadii[v1] * coshRadii[v2] - 
                                sinhRadii[v1] * sinhRadii[v2] * cosAngles[index])
            #if length == 0:
            #	length = 1.01e-10
            self.edgeLengths[index] = length

    def recalcFaces(self):
        """ Reconstruct face triples from triangulation data.

        	Usually called after initializing, which destroys the original face data by
        	removing isolated vertices and subdividing boundary triangles. This recreates
        	the data needed to properly display the surface.
        """
        self.faces = np.empty((len(self.triangles), 3), dtype='int32')
        for index, triangle in enumerate(self.triangles):
            self.faces[index,:] = np.array(triangle.vertices)

    def calcEulerChar(self):
        return len(self.vertices) - len(self.edges) + len(self.triangles)

    def display(self, scalars=None, representation='surface'):
        import enthought.mayavi.mlab as mlab
        if scalars is None:
            return mlab.triangular_mesh(self.x, self.y, self.z, self.faces, representation=representation)
        else:
            return mlab.triangular_mesh(self.x, self.y, self.z, self.faces, scalars=scalars, representation=representation)

def errorFunc(radii, triangulation):
    """ Calculates the sum of squares of errors between initial and
    	current edge lengths. 

    	Does not explicitly calculate edge lengths unless they are outside angle bounds.
    """
    totalError = 0.0
    error = 0.0
    for index, edge in enumerate(triangulation.edges):
        v1, v2 = edge.vertices
        r1, r2 = radii[v1], radii[v2]
        lBar = triangulation.initialEdgeLengths[index]
        if lBar > r1 + r2:
            error = (r1 + r2 - lBar)**2
            totalError += error
        elif math.cosh(lBar) < math.cosh(r1) * math.cosh(r2):
            error = (math.acosh(math.cosh(r1) * math.cosh(r2) ) - lBar)**2
            totalError += error
    return totalError

def errorGrad(radii, triangulation):
    """ Gradient of the error function. """
    grad = np.zeros(len(radii))
    deriv = 0.0
    for vertexIndex, vertex in enumerate(triangulation.vertices):
        for edgeIndex in vertex.incidentEdges:
            v1, v2 = triangulation.edges[edgeIndex].vertices
            if v1 <> vertexIndex:   # swap so that current vertex is v1
                v1, v2 = v2, v1
            r1, r2 = radii[v1], radii[v2]
            lBar = triangulation.initialEdgeLengths[edgeIndex]
            if lBar > r1 + r2:
                deriv = 2*(r1 + r2 - lBar)
                grad[vertexIndex] += deriv
            elif math.cosh(lBar) < math.cosh(r1) * math.cosh(r2):
                lActual = math.acosh(math.cosh(r1) * math.cosh(r2) )
                deriv = (2 * (lActual - lBar) * 
                         math.sinh(r1) * math.cosh(r2) / math.sinh(lActual))
                grad[vertexIndex] += deriv
    return grad



def uniquify(seq, idfun=None):  
    """ Given a list of items, function removes duplicates but preserves the ordering.	
    """
    if idfun is None: 
        def idfun(x): return x 
    seen = {} 
    result = [] 
    for item in seq: 
        marker = idfun(item) 
        # in old Python versions: 
        # if seen.has_key(marker) 
        # but in new ones: 
        if marker in seen: continue 
        seen[marker] = 1 
        result.append(item) 
    return result			

class BoundaryEdgeError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
