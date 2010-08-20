# ricci.py

import triangulate as tri
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as spalg
import scipy.linalg as linalg

class Ricci(tri.Triangulation):
    def __init__(self, points, faces, radii=None, edgeInfo=None):
        tri.Triangulation.__init__(self, points, faces, radii, edgeInfo)
        if radii == None:
            self.calcOptimalRadii()

    def gradientRicciFlow(self, threshold=1.01e-5, maxIterations=1500, stepSize=0.1):
        self.calcCurvatureParameters()
        iterations = 0

        while (supNorm(self.curvatures) > threshold) and (iterations < maxIterations) :
            iterations = iterations + 1
            sinhRadii = np.sinh(self.radii)
            step = -1 * stepSize * self.curvatures * sinhRadii
            for k in range(10):
                if min(self.radii + (0.5**k)*step) >= 0 :
                    self.radii = self.radii + (0.5**k)*step
                    break
                else:
                    print "Tried k=" + str(k) + ", failed. Min r = " + str(min(self.radii + ((0.5**k)*step)))
            self.calcCurvatureParameters()
        print "Iterations = " + str(iterations)
        print "Max curvature = " + str(supNorm(self.curvatures))

    def ricciFlow(self, threshold=1.01e-5, maxIterations=15, hessType=0, supTest=True):
        """ Use Newton's conjugate gradient method to find 0 curvature metric.
        """
        # initialize data
        #self.forceRadiiPositive(tolerance=1.01e-12)
        u = np.log(np.tanh(self.radii / 2))
        self.calcCurvatureParameters()
        if hessType==0:
            self.calcSparseHessianLength()
        iterations = 0

        r0 = np.linalg.norm(self.curvatures)
        stopPoint = (r0 * threshold + threshold) / 2.0
        print "Initial max curvature = " + str(supNorm(self.curvatures))
        print "Initial curvature norm = " + str(r0)
        while (np.linalg.norm(self.curvatures) > stopPoint) and (iterations < maxIterations):
            iterations = iterations + 1
            self.calcHessianParameters()		
            if hessType == 0:
                print "Constructing sparse Hessian...."
                self.calcSparseHessian()
                print "Solving conjugate gradient...."  
                solution = spalg.cg(self.sparseHessian, self.curvatures, tol=threshold/10.0)
                step = solution[0]
            elif hessType == 1: # seems to be least efficient method
                H = spalg.LinearOperator((len(u), len(u)), lambda v: Ricci.hessian(self, v), 
                                         dtype='float64')
                print "Solving conjugate gradient...."  
                solution = spalg.cg(H, self.curvatures)	# use conjugate gradient method to solve for Newton step
                step = solution[0]
            else: #only use for small triangulations
                H = self.hessianMatrix()
                print "Solving conjugate gradient...."  
                step = linalg.solve(H, self.curvatures)
            print "Done."
            u = u - step
            print "Step size = " + str(np.linalg.norm(step))
            # check Wolfe condition (or at least the parts that can be checked)
            oldRadii = self.radii.copy()
            oldCurvatures = self.curvatures.copy()
            #wolfeBound = abs(wolfeFactor * np.dot(step, oldCurvatures))
            for k in range(1,10):
                if max(u) > 0:
                    print "u positive when k = " + str(k)
                    u = u + ((0.5**k)*step)
                else:
                    self.radii = 2* np.arctanh(np.exp(u))
                    #self.forceRadiiPositive()
                    self.calcCurvatureParameters()
                    if supTest and supNorm(self.curvatures) > supNorm(oldCurvatures):
                        print "Curvature sup norm increased when k = " + str(k)
                        u = u + ((0.5**k)*step)
                    elif supTest == False and np.linalg.norm(self.curvatures) > np.linalg.norm(oldCurvatures):
                        print "Curvature 2-norm increased when k = " + str(k)
                        u = u + ((0.5**k)*step)
                    else:
                        break
            if k==9:
                print "Wolfe condition failed - optimization cannot continue"
                self.radii = oldRadii
                self.calcCurvatureParameters()
                print "Max curvature = " + str(supNorm(self.curvatures)), " Norm = ", np.linalg.norm(self.curvatures)
                break
            print "Max curvature = " + str(supNorm(self.curvatures)), " Norm = ", np.linalg.norm(self.curvatures)

        # TO DO: check exit status, set appropriate flags


    def calcCurvatureParameters(self):
        """ Assumes self.radii and self.circleAngles are known, then calculates
        	edge lengths, inner angles and vertex curvatures
        """
        self.calcEdgeLengths()
        self.calcInnerAngles()
        self.calcCurvatures()

    def calcHessianParameters(self):
        sinAngles = np.sin(self.innerAngles)
        cosAngles = np.cos(self.innerAngles)
        cosCircleAngles = np.cos(self.circleAngles)
        sinhEdgeLengths = np.sinh(self.edgeLengths)
        sinhRadii = np.sinh(self.radii)
        coshRadii = np.cosh(self.radii)
        # derivatives of edge lengths w.r.t. radii
        edgeDerivatives = np.zeros((len(self.edges), 2))   
        for index, edge in enumerate(self.edges):
            v1, v2 = edge.vertices
            s1c2 = sinhRadii[v1] * coshRadii[v2]
            s2c1 = sinhRadii[v2] * coshRadii[v1]
            D1 = (s1c2 - s2c1 * cosCircleAngles[index]) / sinhEdgeLengths[index]
            D2 = (s2c1 - s1c2 * cosCircleAngles[index]) / sinhEdgeLengths[index]
            edgeDerivatives[index] = [D1, D2]

        self.innerAngleDerivatives = np.zeros((len(self.triangles), 3, 3))
        for index, triangle in enumerate(self.triangles):
            sinA0 = sinAngles[index, 0]
            sinA1 = sinAngles[index, 1]
            cosA0,cosA1, cosA2 = cosAngles[index]
            sinhL1 = sinhEdgeLengths[triangle.edges[1]]
            sinhL2 = sinhEdgeLengths[triangle.edges[2]]
            sinhR0 = sinhRadii[triangle.vertices[0]]
            sinhR1 = sinhRadii[triangle.vertices[1]]
            sinhR2 = sinhRadii[triangle.vertices[2]]
            dAngle_dEdge = np.matrix([ 
                [1 / (sinA1 * sinhL2), -1*cosA2 / (sinA1 * sinhL2), -1*cosA1 / (sinA1 * sinhL2)],
                [-1*cosA2 / (sinA0 * sinhL2), 1 / (sinA0 * sinhL2), -1*cosA0 / (sinA0 * sinhL2)],
                [-1*cosA1 / (sinA0 * sinhL1), -1*cosA0 / (sinA0 * sinhL1), 1 / (sinA0 * sinhL1)]
            ])
            # find vertices adjacent to each edge
            v01, v02 = self.edges[triangle.edges[0]].vertices
            v10, v12 = self.edges[triangle.edges[1]].vertices
            v20, v21 = self.edges[triangle.edges[2]].vertices
            # Line up vertex indices with proper order in triangle
            if v01 == triangle.vertices[1]:
                v01, v02 = 0, 1
            else:
                v01, v02 = 1, 0
            if v10 == triangle.vertices[0]:
                v10, v12 = 0, 1
            else:
                v10, v12 = 1, 0
            if v20 == triangle.vertices[0]:
                v20, v21 = 0, 1
            else:
                v20, v21 = 1, 0
            dEdge_dRadii = np.matrix([
                [0, edgeDerivatives[triangle.edges[0], v01], edgeDerivatives[triangle.edges[0], v02]],
                [edgeDerivatives[triangle.edges[1], v10], 0, edgeDerivatives[triangle.edges[1], v12]],
                [edgeDerivatives[triangle.edges[2], v20], edgeDerivatives[triangle.edges[2], v21], 0]
            ])
            dRadii_du = np.matrix([
                [sinhR0, 0, 0],
                [0, sinhR1, 0],
                [0, 0, sinhR2]
            ])
            self.innerAngleDerivatives[index] = dAngle_dEdge * dEdge_dRadii * dRadii_du

    def calcSparseHessian(self):
        rowIndices = np.zeros(self.sparseHessianLength, dtype='int32')
        columnIndices = np.zeros(self.sparseHessianLength, dtype='int32')
        data = np.zeros(self.sparseHessianLength, dtype='float64')
        location=0
        for index, vertex in enumerate(self.vertices):
            for triangleIndex in vertex.incidentTriangles:
                # find vertex in triple
                triangle = self.triangles[triangleIndex]
                vI = triangle.vertices.index(index)
                rowIndices[location:location+3]= [index, index, index]
                columnIndices[location:location+3] = triangle.vertices
                data[location:location+3] = -1*self.innerAngleDerivatives[triangleIndex, vI]
                location += 3
        self.sparseHessian = sparse.coo_matrix((data, (rowIndices, columnIndices)), 
                                               shape=(len(self.vertices), len(self.vertices)))

    def hessian(self, p):
        """ Calculates the product of the hessian with the vector p
        """
        output = np.zeros(len(self.vertices))
        for index, vertex in enumerate(self.vertices):
            for triangleIndex in vertex.incidentTriangles:
                # find vertex in triple
                triangle = self.triangles[triangleIndex]
                vI = triangle.vertices.index(index)
                pTriangle = np.array([p[triangle.vertices[0]], p[triangle.vertices[1]],
                                      p[triangle.vertices[2]]])
                output[index] -= np.dot(pTriangle, self.innerAngleDerivatives[triangleIndex, vI])
        return output

    def hessianMatrix(self):
        output = np.zeros((len(self.vertices), len(self.vertices)))
        for i, vertex in enumerate(self.vertices):
            for j,  vertex2 in enumerate(self.vertices):
                for triangleIndex in vertex.incidentTriangles:
                    triangle = self.triangles[triangleIndex]
                    vI = triangle.vertices.index(i)
                    for k, rI in enumerate(triangle.vertices):
                        if rI == j:
                            output[i, j] -= self.innerAngleDerivatives[triangleIndex, vI, k]
        return output

    def calcBoundaryLengths(self):
        """ Adds up lengths of edges along each boundary component.

        	Assumes that self.boundaryCycles has been established (e.g., by calling
        	self.calcBoundaryCycles()).

        	Returns a list of values.
        """
        lengths = []
        for cycle in self.boundaryCycles:
            length = 0
            for e in cycle:
                length += self.edgeLengths[e]
            lengths.append(length)
        return lengths

    def calcSparseHessianLength(self):
        total = 0
        for vertex in self.vertices:
            total = total + 3*(len(vertex.incidentTriangles))
        self.sparseHessianLength = total

    def forceRadiiPositive(self, tolerance=1.01e-12):
        for i in range(len(self.vertices)):
            if self.radii[i] < tolerance:
                self.radii[i] = tolerance

def supNorm(v):
    """ Calculates max |v_i| """
    return max([abs(x) for x in v])
