# this script converts a obj quad mesh into a tri mesh stl

import sys
import os
import argparse
from pathlib import Path

class Vertex:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __sub__(self, other):
        return Vertex(self.x - other.x, self.y - other.y, self.z - other.z)

    def __add__(self, other):
        return Vertex(self.x + other.x, self.y + other.y, self.z + other.z)


def dot(a, b):
    return a.x * b.x + a.y * b.y + a.z * b.z


def cross(a, b):
    return Vertex(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)



if __name__ == "__main__":
    print("[quad2stl] Starting")

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, help="obj", required=True)
    parser.add_argument("-o", type=Path, help="stl", required=True)
    args = parser.parse_args()

    # read obj
    vertices = []
    faces = []
    with open(str(args.i), "r") as obj:
        lines = obj.readlines()
        for line in lines:
            l = line.strip().split()
            if l[0] == 'v':
                vertices.append(Vertex(float(l[1]), float(l[2]), float(l[3])))
            elif l[0] == 'f':
                faces.append([int(x)-1 for x in l[1:]])
    print("[quad2stl] read", len(faces), "faces and", len(vertices),"vertices")

    # write stl
    with open(str(args.o), 'w') as file:
        file.write("solid " + "terrain" + "\n")
        for face in faces:
            normal = cross(vertices[face[1]] - vertices[face[0]], vertices[face[2]] - vertices[face[1]])
            face_indices = [[0,1,2],[0,2,3]]
            for i in range(2):
                file.write("facet normal %f %f %f\n" % (normal.x, normal.y, normal.z))
                file.write("outer loop\n")
                for index in face_indices[i]:
                    vertex = vertices[face[index]]
                    file.write("vertex %f %f %f\n" % (vertex.x, vertex.y, vertex.z))
                file.write("endloop\n")
                file.write("endfacet\n")
        file.write("endsolid " + "terrain" + "\n")

