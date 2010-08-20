# plyFile.py

import os
import sys
import struct

def read(filename):
    """ Reads PLY file filename and creates vertex and triangle arrays """
    #Read face data from file START
    faceFile = open(filename, 'rb')
    #faceOutFile = open('C:/Users/Aaron/Desktop/facePlyTest.txt', 'w')
    count=0
    header=[]
    for count in range(13):
        line=faceFile.readline()
        header.append(line.split())

    print header
    vertexCount = int(header[2][2])
    faceCount = int(header[10][2])
    print vertexCount
    print faceCount

    vertexArray=[]
    xArray=[]
    yArray=[]
    zArray=[]
    faceArray=[]
    colorArray=[]

    for count in range(vertexCount):

        faceBytex=faceFile.read(4)
        faceBytey=faceFile.read(4)
        faceBytez=faceFile.read(4)
        faceByteStuff=faceFile.read(4)
        faceByter=faceFile.read(1)
        faceByteg=faceFile.read(1)
        faceByteb=faceFile.read(1)

        vertex = []


        vertex.append(struct.unpack('>f',faceBytex)[0])
        vertex.append(struct.unpack('>f',faceBytey)[0])
        vertex.append(struct.unpack('>f',faceBytez)[0])
        vertexArray.append(vertex)

        tempArray=[]
        tempArray.append(struct.unpack('>B',faceByter)[0])
        tempArray.append(struct.unpack('>B',faceByteg)[0])
        tempArray.append(struct.unpack('>B',faceByteb)[0])
        colorArray.append(tempArray)

    for count in range(faceCount):
        faceByteSize=faceFile.read(1)
        faceByte1=faceFile.read(4)
        faceByte2=faceFile.read(4)
        faceByte3=faceFile.read(4)

        tempArray=[]
        tempArray.append(struct.unpack('>i',faceByte1)[0])
        tempArray.append(struct.unpack('>i',faceByte2)[0])
        tempArray.append(struct.unpack('>i',faceByte3)[0])
        faceArray.append(tempArray)

    faceFile.close()
    return vertexArray, faceArray, colorArray
