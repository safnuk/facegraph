# gridCreate.py
import random

def gridTShirt(size):
    points = []
    pD = {}
    index = 0
    faces = []
    for i in range(size):
        for j in range(size):
            points.append([float(i), float(j), random.random()])
            pD[(i, j)] = index
            index += 1
    for i in range(size - 1):
        for j in range(size - 1):
            if isHole(i, j, size) == False:
                faces.append([pD[(i,j)], pD[(i, j+1)], pD[(i+1,j+1)]])
                faces.append([pD[(i,j)], pD[(i+1, j+1)], pD[(i+1,j)]])
    return points, faces

def isHole(i, j, size):
    if (j < size/5) or (j >= 4 * size / 5):
        return False
    if (j >= 2*size / 5) and (j < 3*size/5):
        return False

    if j < 2*size/5:
        if (i > size/5) and (i <= 2*size/5):
            return True
        if (i > 3*size/5) and (i <= 4*size/5):
            return True
        else:
            return False

    if (i > 2*size/5) and (i <= 3*size/5):
        return True

    return False

def gridFace(size, radii, centers):
    """ Lists out the vertices and faces of a disc of radius size centered at the origin.
    
    Also removes all holes with given radii and centers.
    radii = list of radii of all holes cut out
    centers = list of [x, y]-coordinates for centers of holes.
    """
    points = []
    pD = {}
    index = 0
    faces = []
    for i in range(-size, size):
        for j in range(-size, size):
            points.append([float(i), float(j), random.random()])
            pD[(i, j)] = index
            index += 1
    for i in range(-size, size - 1):
        for j in range(-size, size - 1):
            if isFaceHole(i, j, size, radii, centers) == False:
                faces.append([pD[(i,j)], pD[(i, j+1)], pD[(i+1,j+1)]])
                faces.append([pD[(i,j)], pD[(i+1, j+1)], pD[(i+1,j)]])
    return points, faces

def isFaceHole(i, j, size, radii, centers):
    #largeCenter = [size/2.0, size/2.0]
    if ((i)**2 + (j)**2) >= (size)**2:
        return True
    for k, c in enumerate(centers):
        if ((i - c[0])**2 + (j - c[1])**2) <= (radii[k]**2):
            return True
    return False