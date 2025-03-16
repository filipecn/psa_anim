import bpy, bmesh

obj = bpy.context.object

mesh = obj.data

bm = bmesh.from_edit_mesh(mesh)


def addMaterial(name, color):
    if name not in bpy.data.materials:
        mat = bpy.data.materials.new(name)
    else:
        mat = bpy.data.materials[name]
    mat.diffuse_color = color
    obj.data.materials.append(mat)

obj.data.materials.clear()
addMaterial("good", (0.4,0.4,0.4,1.0))
addMaterial("bad", (1.0,0.0,0.0,1.0))
addMaterial("small", (0.0,1.0,0.0,1.0))


# adding material and setup nodes:
#mat = bpy.data.materials.new("topo_mat")
#mat.use_nodes = True
#node = mat.node_tree.nodes.new("ShaderNodeAttribute")
#mat.node_tree.links.new(node.outputs['Color'], mat.node_tree.nodes['Principled BSDF'].inputs[0])
#node.attribute_name = "color"

material_slots = obj.material_slots
for material_slot in material_slots:
    print(material_slot.material)
    
for face in bm.faces:
    face.select = True
obj.active_material_index = 0
bpy.ops.object.material_slot_assign()
for face in bm.faces:
    face.select = False

# compute the dihedral angle between two faces
# accross all edges
vertices_to_merge = []
ma = 360
Ma = 0
me = 60000
Me = 0
for edge in bm.edges:
    angle = edge.calc_face_angle(0) * 180 / 3.14151696
    size = edge.calc_length()
    ma = min([ma, angle])
    Ma = max([Ma, angle])
    me = min([me, size])
    Me = max([Me, size])
    if angle > 20 or size < 10:
        # flag both faces
        obj.active_material_index = 1
        if size < 10:
            obj.active_material_index = 2
            for vertex in edge.verts:
                vertices_to_merge.append(vertex)
        faces = []
        for face in edge.link_faces:
            face.select = True
        bpy.ops.object.material_slot_assign()    
        for face in edge.link_faces:
            face.select = False

verts = [*set(vertices_to_merge)]
print("angle range     [%.2f %.2f]" % (ma,Ma))
print("edge size range [%.2f %.2f]" % (me,Me))
print(len(verts))
    
    
    
    
    
#    s = ""
#    for vertex in edge.verts:
#        vertex.select = True
#        s += " " + str(vertex.index)
#        s += ", " + str(vertex.calc_shell_factor())
#    print(s)
#    bpy.ops.mesh.merge(type='CENTER', uvs=False)
#    break
    
    
    #    if edge.calc_length() < 5:
        
#    
#    print(len(faces))
    
#    print(faces[0].color)
#    # for manifolds, we should expect
#    # at most 2 faces per edge
#    if len(faces) != 2:
#        continue
#    # before calculating the angle lets
#    # check if both faces have a consistent
#    # vertex ordering
#    edge.
#    print(edge.link_faces)
##    for face in edge.link_faces:
##        s += str(face.index) + " "
##        
##    print(s)


#for face in bm.faces:
#    print("f = " + str(face.index))
#    print(face.normal)
#    for loop in face.loops:
#        print(loop.vert.co[:])
#    for edge in face.edges:
#        print(edge.index)
#        s = ""
#        for f in edge.link_faces:
#            s += " " + str(f.index)
#        print(s)
#    print("")
