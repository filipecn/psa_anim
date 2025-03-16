# this script creates a releaseArea openfoam file to describe initial conditions for DSL
# the output directory must be the constant directory of a openfoam project

import os
import argparse
from pathlib import Path

foam_header = "/*--------------------------------*- C++ -*----------------------------------*\n" \
              "| =========                 |                                                 |\n" \
              "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n" \
              "|  \\    /   O peration     | Version:  v2012                                 |\n" \
              "|   \\  /    A nd           | Website:  www.openfoam.com                      |\n" \
              "|    \\/     M anipulation  |                                                 |\n" \
              "\\*---------------------------------------------------------------------------*/\n" \
              "FoamFile\n" \
              "{\n" \
              "    version     2.0;\n" \
              "    format      ascii;\n" \
              "    class       dictionary;\n" \
              "    object      setFieldsDict;\n" \
              "}\n" \
              "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
foam_footer = "// ************************************************************************* //\n"

def write_polygon(path, area_type, value):
    s = ""
    with open(path, "r") as file:
        # read obj
        vertices = []
        faces = []
        lines = file.readlines()
        for line in lines:
            l = line.strip().split()
            if l[0] == 'v':
                vertices.append([l[1], l[2], l[3]])
            elif l[0] == 'f':
                p = []
                for i in range(1, len(l)):
                    p.append(int(l[i]) - 1)
                faces.append(p)
        i = 1
        for face in faces:
            s += "            " + area_type + str(i) + "\n"
            i += 1
            s += "            {\n"
            s += "                type polygon;\n"
            s += "                filltype constant;\n"
            s += "                offset (0 0 0);\n"
            s += "                vertices\n"
            s += "                (\n"
            for v in face:
                s += "                    (%s %s 0)\n" % (vertices[v-1][0], vertices[v-1][1])
            # close vertices
            s += "                );\n"
            s += "                value " + str(value) + ";\n"
            # close region
            s += "            }\n"
    return s

def write_release_area(data):
    value = data[0]
    xmin = data[1]
    ymin = data[2]
    xmax = data[3]
    ymax = data[4]
    s = ""
    s += "            {\n"
    s += "                type polygon;\n"
    s += "                filltype constant;\n"
    s += "                offset (0 0 0);\n"
    s += "                vertices\n"
    s += "                (\n"
    s += "                    (%.2f %.2f 0)\n" % (xmin, ymax)
    s += "                    (%.2f %.2f 0)\n" % (xmin, ymin)
    s += "                    (%.2f %.2f 0)\n" % (xmax, ymin)
    s += "                    (%.2f %.2f 0)\n" % (xmax, ymax)
    # close vertices
    s += "                );\n"
    s += "                value " + str(value) + ";\n"
    # close region
    s += "            }\n"
    return s


def write_regions(regions_data):
    s = ""
    s += "        regions\n"
    s += "        (\n"
    number_of_boxes = len(regions_data) // 5
    for i in range(number_of_boxes):
        j = i * 5
        s += "            releaseArea" + str(i) + "\n"
        s += write_release_area(regions_data[j:j + 5])
    # close regions
    s += "        );\n"
    return s

def write_polygons(paths, area_type, values):
    s = ""
    s += "        regions\n"
    s += "        (\n"
    for path in paths:
        s += write_polygon(path, area_type, values[0])
    s += "        );\n"
    return s

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", type=Path, help="openfoam project constant path")
    parser.add_argument("--area", type=float, nargs='+', help="[value, lower corner(x,y), upper corner, value...]")
    parser.add_argument("--entrain", type=float, nargs='+', help="[value, lower corner(x,y), upper corner, value...]")
    parser.add_argument("--default-entrain", type=float, help="default entrain value", default=0.1)
    parser.add_argument("--default-h", type=float, help="default height value", default=0.0001)
    parser.add_argument("--area-files", type=Path, nargs='+', help="[file1 ...]")
    parser.add_argument("--entrain-files", type=Path, nargs='+', help="[file1 ...]")
    args = parser.parse_args()

    print("set release ++++++++++++++++++++++++++++++++++++++")
    print(args.area)
    print(args.entrain)
    print(args.area_files)
    print(args.entrain_files)
    print("++++++++++++++++++++++++++++++++++++++++++++++++++")

    file_content = foam_header

    file_content += "fields\n"
    file_content += "{\n"

    file_content += "    hentrain\n"
    file_content += "    {\n"
    file_content += "        default default [0 1 0 0 0 0 0] " + str(args.default_entrain) + ";\n"
    if args.entrain is not None:
        file_content += write_regions(args.entrain)
    if args.entrain_files is not None:
        file_content += write_polygons(args.entrain_files, "entrainArea", [args.default_entrain])
    # close entrain
    file_content += "    }\n"

    file_content += "    h\n"
    file_content += "    {\n"
    file_content += "        default default [0 1 0 0 0 0 0] " + str(args.default_h) + ";\n"
    if args.area is not None:
        file_content += write_regions(args.area)
    if args.area_files is not None:
        file_content += write_polygons(args.area_files, "releaseArea", [1.0])
    # close h
    file_content += "    }\n"

    # close fields
    file_content += "}\n"

    file_content += foam_footer

    with open(args.o / "releaseArea", "w") as f:
        f.write(file_content)
