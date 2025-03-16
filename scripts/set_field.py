# this script creates a setFieldsDict openfoam file to describe initial conditions for a given field
# the output directory must be the system directory of a openfoam project

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", type=Path, help="openfoam project system path")
    parser.add_argument("-b", type=float, nargs='+', help="[value, lower corner, upper corner, value...]")
    args = parser.parse_args()
    
    file_content = foam_header
    file_content += "defaultFieldValues\n(\n" 
    file_content += "    volScalarFieldValue alpha.snow 0\n"
    file_content += "    volVectorFieldValue U (0 0 0)\n"
    file_content += ");\n"

    file_content += "regions\n(\n"
    
    if args.b is not None:
        file_content += "    boxToCell\n"
        file_content += "    {\n"
        
        number_of_boxes = len(args.b) // 7
        for i in range(number_of_boxes):
            j = i * 7
            file_content += "        box "
            file_content += "(%.2f %.2f %.2f) " % (args.b[j + 1], args.b[j + 2], args.b[j + 3])
            file_content += "(%.2f %.2f %.2f);\n" % (args.b[j + 4], args.b[j + 5], args.b[j + 6])
            file_content += "        fieldValues\n"
            file_content += "        (\n"
            file_content += "            volScalarFieldValue alpha.snow %.2f\n" % args.b[j]
            file_content += "        );\n"

        file_content += "    }\n"


    file_content += ");\n"
    file_content += foam_footer

    with open(args.o / "setFieldsDict", "w") as f:
        f.write(file_content)
