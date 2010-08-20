# waveFront.py

import geodesic
import numpy
import triangulate as tri

class WaveFront():
    def __init__(self, face, antNumber, stepSize=0.01, steps=1000):
        self.stepSize = stepSize
        self.face = face
        self.steps = steps
        self.antNumber = antNumber
        self.waveList = []
        for boundary, cycle in enumerate(face.boundaryCycles):
            antList = []
            for order, edge in enumerate(cycle):
                for gap in range(antNumber):
                    newAnt = geodesic.Ant(boundary, edge, order,
                                     (2*gap + 1) * face.edgeLengths[edge]/(2.0*antNumber), 
                                     face) 
                    antList.append(newAnt)
            self.waveList.append(antList)
    
    def calcCollisionGraph(self, source=None):
        # initialize the list of triangles hit by an ant
        trianglesHit = []
        collisionList = {}
        for triangle in self.face.triangles:
            trianglesHit.append([])
        for wave in self.waveList:
            for ant in wave:
                trianglesHit[ant.currentTriangle].append((ant.boundary, ant.order))
            
        # start marching and looking for collisions
        for k in range(self.steps):
            triangleCheckList = []
            # make one step forward
            for wave in self.waveList:
                for ant in wave:
                    try:
                        newTriangles = ant.stepForward(self.stepSize)
                        if newTriangles is not None:
                            for newTriangle in newTriangles:
                                triangleCheckList.append((newTriangle, ant))
                                trianglesHit[newTriangle].append((ant.boundary, ant.order))
                    except tri.BoundaryEdgeError as e:
                        print e.value	
                        wave.remove(ant)
           # check new triangles entered this round for collisions
            for triangle, ant in triangleCheckList:
                if self.collisionsFound(trianglesHit[triangle]):
                    try:
                        self.waveList[ant.boundary].remove(ant)
                    except ValueError: # do nothing if ant has already been removed
                        pass           # this is for ants which have overstepped past their collision triangle
                    collisionList[triangle] =  trianglesHit[triangle]                   
            # check for surviving ants
            if sum([len(wave) for wave in self.waveList]) == 0:
                break
        # colorize collision triangles
        if source is not None:
            colors = numpy.zeros(len(self.face.vertices))
            for triangle in collisionList:
                for vertex in self.face.triangles[triangle].vertices:
                    colors[vertex] = 1
            source.scalars = colors
        return collisionList
        
    def collisionsFound(self, hitList, maxJump=5):
        if len(hitList) <= 1:
            return False
        ants = []
        boundary = hitList[0][0]
        for currentBoundary, ant in hitList:
            if currentBoundary <> boundary:
                return True  # different boundary ants always indicates a collision
            ants.append(ant)
        ants.sort()
        # look for gaps, which indicates a same-boundary collision
        jumps = [ants[i+1] - ants[i] for i in range(len(ants)-1)]
        if ants[0] > maxJump:  
            if max(jumps) <= maxJump:
                return False
            else:
                return True    
        # the first edge in the boundary cycle is a special case 
        # (continuous with the last edge in the boundary cycle
        else:
            if ants[0] - ants[-1] + len(self.face.boundaryCycles[boundary]) - 1 <= maxJump:
                jumps.remove(max(jumps))
            try:
                if max(jumps) <= maxJump:
                    return False
                else:
                    return True
            except ValueError:  # catch finding max of empty jump list
                return False
        
            