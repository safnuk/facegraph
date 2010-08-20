# gridFile.py



def read(fileName):
    f = open(fileName)
    vertexRef = {}
    faces = []
    points = []
    radii = []
    edgeInfo = {}
    index = 0
    for line in f:
        lineList = line.split(' ')
        if lineList[0] == 'Vertex':
            vertexRef[int(lineList[1])] = index
            index += 1
            point = [float(lineList[2]), float(lineList[3]), float(lineList[4])]
            if len(lineList) == 6:
                radius = float(lineList[5])
                radii.append(radius)
            points.append(point)
        elif lineList[0] == 'Face':
            face = [vertexRef[int(lineList[2])], vertexRef[int(lineList[3])], 
                    vertexRef[int(lineList[4])]]
            faces.append(face)
        elif lineList[0] == 'Edge':
            edge = (int(lineList[1]), int(lineList[2]))
            angle = float(lineList[3])
            length = float(lineList[4])
            edgeInfo[edge] = [angle, length]
    f.close()
    if len(radii) == 0:
        radii = None
    if len(edgeInfo) == 0:
        edgeInfo = None
    return points, faces, radii, edgeInfo


def write(fileName, mesh):
    f = open(fileName, 'w')
    for i in range(len(mesh.vertices)):
        line = ("Vertex " + str(i) + ' ' + str(mesh.x[i]) + ' ' 
                + str(mesh.y[i]) + ' ' + str(mesh.z[i]) 
                + ' ' + str(mesh.radii[i]) + '\n')
        f.write(line)

    for i, face in enumerate(mesh.faces):
        line = ('Face ' + str(i) + ' ' + str(face[0]) + ' ' +
                str(face[1]) + ' ' + str(face[2]) + '\n')
        f.write(line)

    for i, e in enumerate(mesh.edges):
        v0, v1 = e.vertices
        if v0 > v1:
            v0, v1 = v1, v0
        line = ('Edge ' + str(v0) + ' ' + str(v1) 
                + ' ' + str(mesh.circleAngles[i]) 
                + ' ' + str(mesh.edgeLengths[i]) + '\n')
        f.write(line)

    f.close()