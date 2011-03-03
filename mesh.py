# mesh.py

import math
import numpy as np

#import scipy.optimize as opt

MAX_DEGREE = 15

class Vertex:
    """ Class which encodes incidence data at a vertex.

    	Incidence data should be encoded as indices to the positions in lists of
    	vertices, edges and triangles.

        degree is the number of incident edges; boundary == 1 if vertex 
        is on boundary, 0 if it is an interior vertex.

        The incidence data for edges and vertices should line up. 
        i.e. the first incident vertex should be connected to the current vertex
        by the first incident edge, and so on.

    	Note: this class does not record vertex position, i.e. (x,y,z)

    """
    def __init__(self):
        self.degree = 0
        self.boundary = 0
        self.incident_vertices = np.empty((MAX_DEGREE), dtype='int')
        self.incident_triangles = np.empty((MAX_DEGREE), dtype='int')
        self.incident_edges = np.empty((MAX_DEGREE), dtype='int')
        self.incident_vertices.fill(-1)
        self.incident_triangles.fill(-1)
        self.incident_edges.fill(-1)

    def __repr__(self):
        return ('<Vertex: V' + str(self.incident_vertices) + '; E' + 
                str(self.incident_edges) + '; T' + str(self.incident_triangles) + '>')

    def __str__(self):
        return ('<Vertex: V' + str(self.incident_vertices) + '; E' + 
                str(self.incident_edges) + '; T' + str(self.incident_triangles) + '>')

class Edge:
    """ Edge class for encoding connectivity structure. Vertices should be a 2-element
    	list of indices to the vertex list  (endpoints of the edge). 

    	Incident triangles are recorded as indices to the list of triangles.

        If vertices are listed in order i, j then incidentTriangle[0] should 
        be on left when traveling from i to j. If edge is on the boundary, then
        incidentTriangle[1] = -1.

    """
    def __init__(self, vertex1, vertex2):
        self.vertices = np.array([vertex1, vertex2], dtype='int')
        self.incident_triangles = np.array([-1,-1], dtype='int')

    def __repr__(self):
        return '<Edge: V' + str(self.vertices) + '; T' + str(self.incident_triangles) + '>'

    def __str__(self):		
        return '<Edge: V' + str(self.vertices) + '; T' + str(self.incident_triangles) + '>'

    def get_other_vertex(self, v):
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
    	listed in the order (a, b, c) means edges should be listed in the order (A, B, C).

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

        Assumption: Triangulation data has 
            - no isolated vertices
            - no bivalent vertices
            - trivalent and quadvalent vertices cannot be easily eliminated

        In other words, assume that data has first been run through
        Triangulation class.
        
    """
    def __init__(self, points, faces, radii=None, edge_info=None):

        self.faces = np.array(faces, dtype='int')		# data structure used to display surface with mlab

        # initialize data arrays in preparation for moving information into
        # the class
        self.x = np.empty(len(points))
        self.y = np.empty(len(points))
        self.z = np.empty(len(points))
        self.vertices = []
        self.edges = []
        self.triangles = []

        # store (x,y,z) coordinates of points, and initialize vertex list
        for index, point in enumerate(points):
            self.vertices.append(Vertex())
            self.x[index] = point[0]
            self.y[index] = point[1]
            self.z[index] = point[2]


        edge_index = 0
        for triangle_index, face in enumerate(faces):
            vertex_list = [face[0], face[1], face[2]]
            triangle = Triangle(vertex_list)
            self.triangles.append(triangle)
            self._add_incident_triangle(triangle_index, vertex_list)
            edge_index = self._create_new_edges(edge_index, triangle_index,
                                                vertex_list)

    def _add_incident_triangle(self, triangle_index, vertex_list):
        """ Adds incidence data of triangle to vertices listed in
            vertex_list.

        """
        for vertex_index in vertex_list:
            vertex = self.vertices[vertex_index]
            i = vertex.degree - vertex.boundary
            vertex.incident_triangles[i] = triangle_index
            vertex.boundary -= 1  # boundary = degree - # of incident triangles
        #for index, vertex_index in enumerate(vertex_list):
        #    vertex = self.vertices[vertex_index]
        #    other_vertices = [vertex_list[((index+1) % 3)],
        #                      vertex_list[((index+2) % 3)]]
        #    for ov in other_vertices:
        #        i = 0
        #        duplicate = False
        #        while vertex.incident_vertices[i] <> -1:
        #            if vertex.incident_vertices[i] == ov:
        #                duplicate = True 
        #                break
        #            i += 1
        #        if not duplicate:
        #            vertex.incident_vertices[i] = ov
        #            vertex.boundary += 1

    def _create_new_edges(self, edge_index, triangle_index, vertex_list):
        """ Creates any new edges coming from the triangle with vertices given
        by vertex_list. It first checks to see if each edge has been
        previously constructed, and updates and returns edge_index reflecting
        the number of new edges created.

        It also records incidence data on all relevant vertices, triangles and
        edges.

        """
        edge_ref = -1
        for i in range(3):
            vertex1_index = vertex_list[i]
            vertex2_index = vertex_list[(i+1) % 3]
            vertex1 = self.vertices[vertex1_index]
            vertex2 = self.vertices[vertex2_index]
            duplicate = False
            for j in range(vertex2.degree):
                if vertex2.incident_vertices[j] == vertex1_index:
                    duplicate = True
                    edge_ref = vertex2.incident_edges[j]
            if not duplicate:
                edge_ref = edge_index
                edge_index += 1
                # create new edge
                edge = Edge(vertex1_index, vertex2_index)
                edge.incident_triangles[0] = triangle_index
                self.edges.append(edge)
                # record incidence data of vertices
                vertex1.incident_vertices[vertex1.degree] = vertex2_index
                vertex1.incident_edges[vertex1.degree] = edge_ref
                vertex1.degree += 1
                vertex1.boundary += 1
                vertex2.incident_vertices[vertex2.degree] = vertex1_index
                vertex2.incident_edges[vertex2.degree] = edge_ref
                vertex2.degree += 1
                vertex2.boundary += 1
            else:
                edge = self.edges[edge_ref]
                edge.incident_triangles[1] = triangle_index
            # add edge ref to triangle
            edge_position = (i+2) % 3
            self.triangles[triangle_index].edges[edge_position] = \
                edge_ref
        return edge_index

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
