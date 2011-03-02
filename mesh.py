# mesh.py

import math
import numpy as np

#import scipy.optimize as opt

class Vertex:
    """ Class which encodes incidence data at a vertex.

    	Incidence data should be encoded as index data to the positions in lists of
    	vertices, edges and triangles.

    	Note: this class does not record vertex position, i.e. (x,y,z)
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
        self.vertices = np.array(vertices)
        self.vertices.sort()
        self.incidentTriangles = np.array([-1,-1])

    def __repr__(self):
        return '<Edge: V' + str(self.vertices) + '; T' + str(self.incidentTriangles) + '>'

    def __str__(self):		
        return '<Edge: V' + str(self.vertices) + '; T' + str(self.incidentTriangles) + '>'

    def getOtherVertex(self, v):
        if v == self.vertices[0]:
            return self.vertices[1]
        if v == self.vertices[1]:
            return self.vertices[0]
        raise VertexIndexError("Vertex is not incident to edge")    

class Triangle:
    """ Triangle class used to encode connectivty structure. Vertices should be a 
    	counter-clockwise list of indices. Edge indices are also listed counterclockwise,
    	with the edge being opposite to the vertex.

    	i.e. if a triangle has vertices a, b, c with opposite edges A, B, C, then vertices 
    	listed in the order (a, b, c) measn edges should be listed in the order (A, B, C).
    """
    def __init__(self, vertices):
        self.vertices = np.array(vertices)
        self.edges = np.array([-1, -1, -1])

    def __repr__(self):
        return '<Triangle: V' + str(self.vertices) + '; E' + str(self.edges) + '>'

    def __str__(self):		
        return '<Triangle: V' + str(self.vertices) + '; E' + str(self.edges) + '>'

class Mesh:
    """ Class which encodes the incidence relations of vertices, edges and faces 
    	of a triangulation initialized with a list of points [[x0, y0, z0], ...] 
    	and a list of triples of indices [[v0, v1, v2], ...] in counterclockwise order. 

    	The indices for the points are implied from the given ordering.  
    """
    def __init__(self, points, faces, radii=None, edgeInfo=None):

        self.faces = np.array(faces)		# data structure used to display surface with mlab

        # initialize data arrays in preparation for moving information into
        # the class
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
class VertexIndexError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class BoundaryEdgeError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
