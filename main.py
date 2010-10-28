# main.py

import geodesic
import numpy as np
import triangulate as tri
import waveFront
import ricci

import gridFile

def main(fileName):
    print "Loading ", fileName 
    face = ricci.Ricci(*gridFile.read(fileName))
    face.calcCurvatureParameters()
    print "Running ants..."
    runWaves(face, 4, stepSize = 0.002, steps=2500)

def runWaves(face, antNumber, stepSize=0.01, steps=1000, source=None):
    front = waveFront.WaveFront(face, antNumber, stepSize, steps)
    collisions =  front.calcCollisionGraph(source)
    print "Vertex list"
    for triangle, hitList in collisions.items():
        hitList.sort()
        boundary = hitList[0][0]
        count = 1
        for curBoundary, edgePosition in hitList:
            if curBoundary <> boundary:
                count += 1
                boundary = curBoundary
        if count > 2:
            print "Triangle: ", triangle
    return collisions
         

def run(mesh, antNumber, stepSize=0.01, steps=1000, delay=0.0):
    import enthought.mayavi.mlab as mlab
    antList = []
    for boundary, cycle in enumerate(mesh.boundaryCycles):
        for order, edge in enumerate(cycle):
            for gap in range(antNumber):
                newAnt = geodesic.Ant(boundary, edge, order, (2*gap + 1) * mesh.edgeLengths[edge]/(2.0*antNumber), mesh)
                antList.append(newAnt)
    colors = np.zeros(len(mesh.vertices))
    surface = mesh.display(scalars=colors)
    source = surface.mlab_source
    for k in range(steps):
        """
		if k < 10:
			filename = "./Movie/Face000" + str(k) + ".png"
		elif k < 100:
			filename = "./Movie/Face00" + str(k) + ".png"
		elif k < 1000:
			filename = "./Movie/Face0" + str(k) + ".png"
		else:
			filename = "./Movie/Face" + str(k) + ".png"
		"""
        for pos, a in enumerate(antList):
            e = a.config.entryEdge
            v1, v2 = mesh.edges[e].vertices
            if colors[v1] == 0: 
                colors[v1] = a.boundary + 1
            if colors[v2] == 0:
                colors[v2] = a.boundary + 1
            if colors[v1] <> a.boundary + 1 or colors[v2] <> a.boundary + 1:
                antList.remove(a)
            else:
                try:
                    a.stepForward(stepSize)
                except tri.BoundaryEdgeError as e:
                    print e.value	
                    antList.remove(a)
        if k%30 == 0:
            source.scalars = colors
        #mlab.savefig(filename)
        if len(antList) == 0:
            break

    source.scalars = colors

def removeEyes(mesh, xLeftEye, yLeftEye, xRightEye, yRightEye, radius):
    cutList = cutCircleList(mesh, xLeftEye, yLeftEye, radius) + cutCircleList(mesh, xRightEye, yRightEye, radius) 
    cutList.sort()
    cutList.reverse()
    points = np.transpose(np.vstack((mesh.x, mesh.y, mesh.z)))
    faces = list(mesh.faces)
    for t in cutList:
        faces.pop(t)
    return points, faces

def cutCircleList(mesh, xCenter, yCenter, radius):
    cutList = []
    radiusSquared = radius**2
    for v, vertex in enumerate(mesh.vertices):
        xVertex = mesh.x[v]
        yVertex = mesh.y[v]
        distanceSquared = (xVertex - xCenter)**2 + (yVertex - yCenter)**2
        if distanceSquared < radiusSquared:
            cutList.extend(vertex.incidentTriangles)
    cutList = tri.uniquify(cutList)
    return cutList


def colorEyes(mesh, xLeftEye, yLeftEye, xRightEye, yRightEye, radius):
    colors = np.zeros(len(mesh.vertices))
    radiusSquared = radius**2
    for v, vertex in enumerate(mesh.vertices):
        x = mesh.x[v]
        y = mesh.y[v]
        d1Sq = (x - xLeftEye)**2 + (y - yLeftEye)**2
        d2Sq = (x - xRightEye)**2 + (y - yRightEye)**2
        if d1Sq < radiusSquared or d2Sq < radiusSquared:
            colors[v] = 1
    return colors

if __name__ == "__main__":
    main("./Data/face1Cut-flat2.mesh")
