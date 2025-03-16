# trace generated using paraview version 5.11.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

import os
import argparse
from pathlib import Path

#### import the simple module from the paraview
from paraview.simple import *

def LOG(m):
    print("[render_with_paraview] " + str(m))

# detect all folders with VTK directories
def findVTKRoots(root):
    r = []
    for subdir, dirs, files in os.walk(root):
        for dir in dirs:
            if dir == "VTK":
                r.append(root / subdir / dir)
    return r

def listVTKFiles(vtk_roots):
    r = {}
    for vtk_root in vtk_roots:
        r[vtk_root] = list(vtk_root.glob('*.vtm'))
    return r

def renderSim(name, root, vtm_files, output, annotation):
    # create a new 'XML MultiBlock Data Reader'
    vtm = XMLMultiBlockDataReader(registrationName='sim.vtm*', FileName=[str(f) for f in vtm_files])
    # Pick fields to load
    vtm.CellArrayStatus = ['U', 'alpha.snow', "Uinj"]
    vtm.PointArrayStatus = []


    # get animation scene
    animationScene1 = GetAnimationScene()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    vtmDisplay = Show(vtm, renderView1, 'GeometryRepresentation')

    # get the material library
    materialLibrary1 = GetMaterialLibrary()
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # set scalar coloring
    ColorBy(vtmDisplay, ('CELLS', 'alpha.snow'))
    
    # show color bar/color legend
    vtmDisplay.SetScalarBarVisibility(renderView1, True)
    
    animationScene1.Play()

    # set scalar coloring
    ColorBy(vtmDisplay, ('CELLS', 'vtkCompositeIndex'))
    ## set scalar coloring
    ColorBy(vtmDisplay, ('CELLS', 'alpha.snow'))

    
    # injection velocity
    # create a new 'Glyph'
    #glyph1 = Glyph(registrationName='Glyph1', Input=vtm, GlyphType='Arrow')
    #glyph1.OrientationArray = ['POINTS', 'Uinj']
    #glyph1.ScaleArray = ['POINTS', 'Uinj']
    #glyph1.ScaleFactor = 1.64100036621094
    #glyph1.GlyphTransform = 'Transform2'
    

    # show data in view
    # glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')
    # show color bar/color legend
    #glyph1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get layout
    layout1 = GetLayout()
    
    # layout/tab size in pixels
    layout1.SetSize(1584, 1403)

    # create a new 'Annotate Time Filter'
    annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=vtm)

    # Properties modified on annotateTimeFilter1
    annotateTimeFilter1.Format = annotation + ': {time:f}'

    # show data in view
    annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1, 'TextSourceRepresentation')


    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [70, -476, 128]
    renderView1.CameraFocalPoint = [158, -45.5, 133]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 140

    # save animation
    SaveAnimation(str(output) + "/" + name + '.avi', renderView1, ImageResolution=[1584, 1400],
        FrameRate=20,
        FrameWindow=[0, 400])

    # destroy annotateTimeFilter1
    Delete(annotateTimeFilter1)
    del annotateTimeFilter1
    Delete(renderView1)
    del renderView1
    Delete(vtmDisplay)
    del vtmDisplay
    Delete(animationScene1)
    del animationScene1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, help="simulations root")
    parser.add_argument("-o", type=Path, help="rendering output dir")
    args = parser.parse_args()

    if not args.i.exists() or not args.o.exists():
        LOG("invalid input dir")
        exit(1)

    param_lines = []
    # read parameters file 
    if (args.i / "parameters.txt").exists():
        with open(args.i / "parameters.txt", "r") as parmFile:
            param_lines = parmFile.readlines()
        
    VTK_folders = findVTKRoots(args.i.absolute())
    VTK_files = listVTKFiles(VTK_folders)
    
    sim_id = 0
    for root in VTK_files:
        # find parameters
        parms = ""
        for parm in param_lines:
            # get path
            path = parm.split()[0]
            if path in str(root):
                parms = parm

        renderSim("sim_" + str(sim_id), root, VTK_files[root], args.o.absolute(), parms)
        sim_id += 1

