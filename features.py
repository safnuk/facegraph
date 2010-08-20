# features.py
import triangulate as tri
import numpy as np


def eyePoints(mesh, vIndex, threshold=0.01):
    eyeList = [vIndex]
    curLength = 0
    while curLength < len(eyeList):
        curLength = len(eyeList)
        eyeCopy = eyeList[:]
        for v in eyeCopy:
            height = mesh.z[v]
            vertex = mesh.vertices[v]
            for w in vertex.incidentVertices:
                if height - mesh.z[w] < threshold:
                    eyeList.append(w)
        eyeList = tri.uniquify(eyeList)
    return eyeList

def cutVertices(mesh, vertexList):
    cutList = []
    for v in vertexList:
        vertex = mesh.vertices[v]
        cutList.extend(vertex.incidentTriangles)
    cutList = tri.uniquify(cutList)
    cutList.sort()
    cutList.reverse()
    points = np.transpose(np.vstack((mesh.x, mesh.y, mesh.z)))
    faces = list(mesh.faces)
    for t in cutList:
        faces.pop(t)
    return points, faces
