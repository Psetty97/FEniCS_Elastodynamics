import meshio
import os
import sys
import numpy as np


def convertMeshXML(filename):
        mFile = meshio.read(filename)
        pts, clls = mFile.points, mFile.cells
        msh = meshio.Mesh(pts[:,0:2], clls)
        meshio.write(os.path.splitext(filename)[0]+'.xml',msh)


def convertMeshXDMF(filename):
    msh = meshio.read(filename)
    meshio.write(os.path.splitext(filename)[0]+'.xdmf',meshio.Mesh(points=msh.points[:,0:2],
                                                                   cells={"triangle": msh.cells["triangle"]}))
    meshio.write(os.path.splitext(filename)[0]+'_boundary'+'.xdmf',
                 meshio.Mesh(points=msh.points[:,0:2],
                             cells={"line": msh.cells["line"]},
                             cell_data={"line": {"boundary": msh.cell_data["line"]["gmsh:physical"]}}))
    meshio.xdmf.read(os.path.splitext(filename)[0]+'.xdmf')


def findFiles(path,ext):
    fname = []
    for root,d_names,f_names in os.walk(path):
        for f in f_names:
            if f.endswith(ext):
                fname.append(os.path.join(root, f))
    return fname

def convertMesh(file, convert_type):
    input_type = os.path.splitext(file)[1]
    if input_type == ".msh":
        if convert_type == "xml":
            convertMeshXML(file)
            converted_file = os.path.splitext(file)[0]+'.xml'
        elif convert_type == "xdmf":
            convertMeshXDMF(file)
            converted_file = os.path.splitext(file)[0]+'.xdmf'
        else:
            raise ValueError("conversion of mesh to ." +convert_type+" is not supported")
    elif input_type == ".xml" or input_type == ".xdmf":
        converted_file = file
    else:
        raise ValueError("reading input mesh of type " +input_type+" is not supported")
    return converted_file

