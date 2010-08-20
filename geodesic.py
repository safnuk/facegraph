# geodesic.py
import math
import triangulate as tri

class Ant():
    def __init__(self, boundary, edge, order, distanceFromVertex, mesh):
        # initial position
        self.boundary = int(boundary)
        self.startingEdge = int(edge)
        self.distanceFromVertex = float(distanceFromVertex)
        self.order = order
        self.mesh = mesh

        self.distanceTraveled = 0
        self.currentTriangle = mesh.edges[edge].incidentTriangles[0]
        self.config = Geodesic(self.currentTriangle, edge, distanceFromVertex, math.pi/2, 0.0, mesh)

    def __repr__(self):
        return ('<Ant: B' + str(self.boundary) + '; E' + 
                str(self.startingEdge) + '; Dist:' + str(self.distanceTraveled) + '; T: ' + 
                str(self.currentTriangle)+'>')

    def __str__(self):
        return ('<Ant: B' + str(self.boundary) + '; E' + 
                str(self.startingEdge) + '; Dist:' + str(self.distanceTraveled) + '; T: ' + 
                str(self.currentTriangle)+'>')

    def stepForward(self, distance):
        self.config.localDistanceTraveled += distance
        self.distanceTraveled += distance
        trianglesEntered = []
        if self.config.localDistanceTraveled >= self.config.totalLength:
            self.enterNextTriangle(self.config.localDistanceTraveled - self.config.totalLength)
            trianglesEntered.append(self.currentTriangle)
            while self.config.localDistanceTraveled >= self.config.totalLength:
                self.enterNextTriangle(self.config.localDistanceTraveled - self.config.totalLength)
                trianglesEntered.append(self.currentTriangle)
            return trianglesEntered
        else:
            return None

    def enterNextTriangle(self, travelDistance):
        edge = self.config.exitEdge
        triangles = self.mesh.edges[edge].incidentTriangles[:]
        triangles.remove(self.config.triangle)
        if len(triangles) == 0:
            message = "Ant ", self, " ran into boundary edge " + str(edge) + '.'
            raise tri.BoundaryEdgeError(message)
        self.currentTriangle = triangles[0]
        newEntryPoint = self.config.exitPoint
        newEntryAngle = self.config.exitAngle
        #print "Entering triangle ", triangles[0], " through edge ", edge, " with data ", newEntryPoint, newEntryAngle
        self.config = Geodesic(self.currentTriangle, edge, newEntryPoint, newEntryAngle, travelDistance, self.mesh)
        # WARNING - does not account for possibility of stepping entirely through new triangle by excess 
        # step distance

class Geodesic():
    """ Encodes the configuration of a geodesic passing through a triangle.
    	Initialized with the entry edge, distance from vertex v0 (assuming orientation
    	inherited from triangle makes the edge v0 --> v1 ) and entry angle.
    	If entry point is O and exit point is P the angle is the counterclockwise
    	rotation from Ov1 to OP.

    	Encodes the exit edge vi ---> vj, point of exit (measured as distance from vj, 
    	also matched distance measurement for next entry),
    	exit angle  from Pvj to PO (which matches entry angle for next triangle),
    	total length of geodesic within triangle, and current position in triangle.
    """
    def __init__(self, triangleIndex, edgeIndex, distanceFromVertex, entryAngle, 
                 localDistanceTraveled, mesh):
        self.triangle = triangleIndex
        self.entryEdge = edgeIndex
        self.entryAngle = entryAngle
        self.entryPoint = distanceFromVertex
        self.localDistanceTraveled = localDistanceTraveled

        cos = math.cos
        acos = math.acos
        sin = math.sin
        cosh = math.cosh
        sinh = math.sinh
        acosh = math.acosh
        asinh = math.asinh
        triangle = mesh.triangles[triangleIndex]
        edge = mesh.edges[edgeIndex]
        vertices = mesh.findOrientedEdgeVertices(edgeIndex, triangleIndex)
        indices = [triangle.vertices.index(x) for x in vertices]
        angles = [mesh.innerAngles[triangleIndex, i] for i in indices]
        edges = [triangle.edges[i] for i in indices]
        edgeLengths = [mesh.edgeLengths[x] for x in edges]
        #print "t, e, v, i, a, es, eLs"
        #print triangle, edge, vertices
        #print indices, angles, edges, edgeLengths
        #print "entry point, angle = ", distanceFromVertex, entryAngle
        # construct triangle for geodesic intersection edge v1 --> v2
        if (entryAngle + angles[1]) >= math.pi:
            cosGamma = 2
        else:
            cosGamma = (-1*cos(entryAngle) * cos(angles[1]) 
                        + sin(entryAngle)*sin(angles[1])*cosh(edgeLengths[2] - distanceFromVertex))
        #print "cosGamma = ", cosGamma
        if abs(cosGamma) <= 1:
            gamma = acos(cosGamma)
            #print "gamma = ", gamma
            coshR = (cos(entryAngle) + cos(angles[1])*cosGamma) / (sin(angles[1]) * sin(gamma))
            #print "coshR = ", coshR
            R = acosh(coshR)
            if R > edgeLengths[0]:
                gamma = None
        else:
            gamma = None

        # if necessary, use triangle from geodesic intersection v2 --> v0
        if gamma is None:
            self.exitEdge = edges[1]
            cosDelta = (-1*cos(angles[0]) * cos(math.pi - entryAngle) 
                        + sin(angles[0]) * sin(math.pi - entryAngle) * cosh(distanceFromVertex))
            delta = acos(cosDelta)
            coshS = (cos(math.pi - entryAngle) + cos(angles[0]) * cosDelta) / (sin(angles[0]) * sin(delta))
            S = acosh(coshS)
            D = asinh(sinh(S) * sin(angles[0]) / sin(math.pi - entryAngle))
            self.totalLength = D
            self.exitPoint = S
            self.exitAngle = delta

        else:
            self.exitEdge = edges[0]
            D = asinh(sinh(R) * sin(angles[1]) / sin(entryAngle))
            self.totalLength = D
            self.exitPoint = edgeLengths[0] - R
            self.exitAngle = math.pi - gamma


