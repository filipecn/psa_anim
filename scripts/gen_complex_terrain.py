# this script generates a ramp terrain surface with irregularities along the x axis
# the output is a stl file that can be extruded later to create the simulation mesh
# The mesh can be inclined. The 2D scheme follows
#    .
#        . w
#             .  
#          a (   .
# where
#    w is the mesh width (inclination direction)
#    a is the inclination angle

import math
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


def rotate(vertex, xp, xpw, a):
    # here we rotate in y axis
    t = vertex - Vertex(xp, 0, 0)
    r = Vertex(t.x * math.cos(a) + t.z * math.sin(a), t.y, t.z * math.cos(a) - t.x * math.sin(a))
    return r + Vertex(xpw, 0, 0)


def vertex_index(ny, i, j):
    return i * (ny + 1) + j


def write_obj(filename, vertices, faces):
    with open(filename, 'w') as file:
        for v in vertices:
            file.write('v %f %f %f\n' % (v.x, v.y, v.z))
        for f in faces:
            file.write('f %s\n' % ' '.join([str(i + 1) for i in f]))


def write_stl(filename, vertices, faces, name):
    with open(filename, 'w') as file:
        file.write("solid " + name + "\n")
        for face in faces:
            normal = cross(vertices[face[1]] - vertices[face[0]], vertices[face[2]] - vertices[face[1]])
            file.write("facet normal %f %f %f\n" % (normal.x, normal.y, normal.z))
            file.write("outer loop\n")
            for index in face:
                vertex = vertices[index]
                file.write("vertex %f %f %f\n" % (vertex.x, vertex.y, vertex.z))
            file.write("endloop\n")
            file.write("endfacet\n")
        file.write("endsolid " + name + "\n")


def gaussian(a, c, r):
    return a * math.exp(-r * r / (2 * c * c))


def distance(xa, ya, xb, yb):
    return math.sqrt((xa - xb) ** 2 + (ya - yb) ** 2)

def func0(x, y):
    return 0

def func1(x, y):
    if x < 100:
        return 0
    return 4 * math.sin(math.pi * x / 10)

def func2(x, y):
    # one bump every 50m, intercalation at y +- 15m
    #
    #    y       O         O
    #    |__x    
    #    |            O
    #
    a = 8
    c = 10
    cx = []
    cy = []
    bump_count = 5
    for ce in range(bump_count):
        cx.append(ce * 80 + 100)
        cy.append(0)
        # cx.append(ce * 50 + 320)
        # if ce % 2:
        #    cy.append(-15)
        #else:
        #    cy.append(15)
    value = 0
    for i in range(len(cx)):
        r = distance(x, y, cx[i], cy[i])
        value = max(value, gaussian(a*(0.2*(i + 1)), c, r))
    return value

def func3(x, y):
    
    value = 0

    a = 8
    c = 10
    r = distance(x, y, x, 0)
    value_line = -gaussian(a, c, r)
    # river sin
    a = 8
    c = 10
    river_pos_x = x
    river_pos_y = 50 * math.sin(x / (10 * math.pi))
    r = distance(x, y, river_pos_x, river_pos_y)
    value_sin = -gaussian(a, c, r)
    # river -sin
    a = 8
    c = 10
    river_pos_x = x
    river_pos_y = -50 * math.sin(x / (10 * math.pi))
    r = distance(x, y, river_pos_x, river_pos_y)
    value_cos = -gaussian(a, c, r)

    joint_radius = 5

    if x < 300:
        value = min([value_line, value])
    #elif x < joints[0]:
    elif x < 695:
        value = min([value_sin, value])
    else:
        a = 8
        c = 10
        r = distance(x, y, 700, y)
        value = min([value, gaussian(a, c, r) - 8])
        value = min([value_line, value])

        bumps = [[750, 0], [825, -40], [825, 40], [900, -80], [900, 0], [900, 80]]
        for b in bumps:
            a = 10
            c = 12
            r = distance(x, y, b[0], b[1])
            value = max([gaussian(a, c, r) - 8, value])

    #elif x < joints[1] - joint_radius:
    #    value = min([value_sin, value])
    #    value = min([value_line, value])
    #elif x < joints[1]:
    #    value = min([value_line, value])
    #elif x < joints[2] - joint_radius:
    #    value = min([value, value_cos])
    #    value = min([value_sin, value])
    #elif x < joints[2]:
    #    value = min([value_line, value])
    #    value = min([value, value_cos])
    #    value = min([value_sin, value])
    #else:
    #    value = min([value_line, value])
    return value

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--stl', type=Path, help="output stl file", required=True)
    parser.add_argument('--obj', type=Path, help="output obj file")
    parser.add_argument('-a', type=float, help="angle in degrees", default=0)
    parser.add_argument('-w', type=float, help="width - x direction")
    parser.add_argument('-d', type=float, help="depth - y direction", default=1)
    parser.add_argument('-c', type=float, help="cell size", default=1)
    parser.add_argument('-p', type=str, help="name", default="terrain")
    parser.add_argument('-f', type=int, help="function", default=1)
    args = parser.parse_args()

    # angle in radians
    angle = math.radians(args.a)

    # compute angle elevation
    e = args.w * math.sin(angle)

    # get block dimensions
    w = args.w * math.cos(angle)
    d = args.d
    c = args.c

    print("w = " + str(w))
    print("e = " + str(e))

    # the idea is to
    # 1 - generate the plane with no inclination
    # 2 - displace vertices by a given function
    # 3 - rotate everything by the slope angle

    ny = int(d // c)
    nx = int(args.w // c)
    hy = d / 2

    func = func0
    if args.f == 0:
        func = func0
    elif args.f == 1:
        func = func1
    elif args.f == 2:
        func = func2
    elif args.f == 3:
        func = func3
    else:
        print("[gen_complex_terrain] invalid func option")
        exit()


    vertices = []
    for x in range(nx + 1):
        for y in range(ny + 1):
            vertices.append(rotate(Vertex(x * c, y * c - hy, func(x * c, y * c - hy)), args.w, w, angle))

    faces = []
    quads = []
    for i in range(nx):
        for j in range(ny):
            faces.append([vertex_index(ny, i, j), vertex_index(ny, i + 1, j + 1), vertex_index(ny, i + 1, j)])
            faces.append([vertex_index(ny, i, j), vertex_index(ny, i, j + 1), vertex_index(ny, i + 1, j + 1)])
            quads.append([
                vertex_index(ny, i, j),
                vertex_index(ny, i, j + 1),
                vertex_index(ny, i + 1, j + 1),
                vertex_index(ny, i + 1, j)
            ])

    write_stl(args.stl, vertices, faces, args.p)
    write_obj(args.obj, vertices, quads)
