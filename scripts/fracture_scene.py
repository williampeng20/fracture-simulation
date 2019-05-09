import bpy
import mathutils
import bmesh
import lattice
import os

# To start from console for console output:
# "/Applications/Blender/blender.app/Contents/MacOS/blender"

def clear_scene():
    for obj in bpy.data.objects:
        obj.select = True
    bpy.ops.object.delete()
    lattice.SCENE_OBJECTS = []

clear_scene()

scene = bpy.context.scene

#bpy.ops.mesh.primitive_cube_add(location=(0,0,10))
cube_vertices = [
    mathutils.Vector((1.0, 1.0, -1.0)),
    mathutils.Vector((1.0, -1.0, -1.0)),
    mathutils.Vector((-1.0, -1.0, -1.0)),
    mathutils.Vector((-1.0, 1.0, -1.0)),
    mathutils.Vector((1.0, 1.0, 1.0)),
    mathutils.Vector((1.0, -1.0, 1.0)),
    mathutils.Vector((-1.0, -1.0, 1.0)),
    mathutils.Vector((-1.0, 1.0, 1.0))]
cube_faces = [
    (0, 1, 2, 3),
    (4, 7, 6, 5),
    (0, 4, 5, 1),
    (1, 5, 6, 2),
    (2, 6, 7, 3),
    (4, 0, 3, 7)]
#cube_vertices = [vertex+mathutils.Vector((-5,-5,10)) for vertex in cube_vertices]
cube_mesh_data = bpy.data.meshes.new("cube_mesh_data")
cube_mesh_data.from_pydata(cube_vertices, [], cube_faces)
cube_mesh_data.update()
cube_obj = bpy.data.objects.new("Cube", cube_mesh_data)
cube_obj.location += mathutils.Vector((0,0,10))
scene.objects.link(cube_obj)

cube = lattice.SceneObject(0.5, [vertex+mathutils.Vector((0,0,10)) for vertex in cube_vertices], cube_obj, False)
lattice.SCENE_OBJECTS.append(cube)

'''obj_path = "/Users/williampeng/Documents/CS-184/fracture-simulation/objects/freedom7.obj"
freedom7_obj = bpy.ops.import_scene.obj(filepath=obj_path)'''

bpy.ops.mesh.primitive_cube_add()
for obj in bpy.context.selected_objects:
    obj.name = "Plane"
plane_vertices = [
    mathutils.Vector((8.0, 8.0, 0.0)),
    mathutils.Vector((8.0, -8.0, 0.0)),
    mathutils.Vector((-8.0, -8.0, 0.0)),
    mathutils.Vector((-8.0, 8.0, 0.0)),
    mathutils.Vector((8.0, 8.0, 0.5)),
    mathutils.Vector((8.0, -8.0, 0.5)),
    mathutils.Vector((-8.0, -8.0, 0.5)),
    mathutils.Vector((-8.0, 8.0, 0.5))]
for i in range(0, len(plane_vertices)):
    bpy.data.objects['Plane'].data.vertices[i].co.x *= 8
    bpy.data.objects['Plane'].data.vertices[i].co.y *= 8
    bpy.data.objects['Plane'].data.vertices[i].co.z = bpy.data.objects['Plane'].data.vertices[i].co.z / 4 + 0.25

plane = lattice.SceneObject(0.5, plane_vertices, bpy.data.objects['Plane'], False)
lattice.SCENE_OBJECTS.append(plane)

#this is the plane
bpy.ops.mesh.primitive_plane_add()
for obj in bpy.context.selected_objects:
    obj.name = "Glass"
glass_obj = bpy.data.objects['Glass']
'''glass_vertices = [
    mathutils.Vector((4.0, 4.0, 0.0)),
    mathutils.Vector((4.0, -4.0, 0.0)),
    mathutils.Vector((-4.0, -4.0, 0.0)),
    mathutils.Vector((-4.0, 4.0, 0.0)),
    mathutils.Vector((4.0, 4.0, 0.5)),
    mathutils.Vector((4.0, -4.0, 0.5)),
    mathutils.Vector((-4.0, -4.0, 0.5)),
    mathutils.Vector((-4.0, 4.0, 0.5))
]'''
for i in range(0, 4):
    glass_obj.data.vertices[i].co.x *= 4
    glass_obj.data.vertices[i].co.y *= 4
    #glass_obj.data.vertices[i].co.z = bpy.data.objects['Plane'].data.vertices[i].co.z / 4 + 0.25
glass_obj.location += mathutils.Vector((0,0,3))
glass = lattice.glassPane(mathutils.Vector((-4.0, -4.0, 0.0)), mathutils.Vector((4.0, 4.0, 0.0)), glass_obj)
lattice.SCENE_OBJECTS.append(glass)

#glass.generate_fractures(mathutils.Vector((0,0,0)), 5)
test_pts = [mathutils.Vector((0.0, 0.0, 0.0)), mathutils.Vector((0.30522313714027405, 0.7625949382781982, 0.0)), mathutils.Vector((0.28216731548309326, -0.2529039978981018, 0.0)), mathutils.Vector((-0.0072681717574596405, -3.709054946899414, 0.0)), mathutils.Vector((-1.6303175687789917, -0.002151413122192025, 0.0)), mathutils.Vector((0.4423271715641022, 1.6230614185333252, 0.0))]
glass.generate_fractures(mathutils.Vector((0,0,0)), 5, test_pts)

'''fragment_verts = [
    mathutils.Vector((1.0, 0.0, 0.0)),
    mathutils.Vector((0.0, 1.0, 0.0)),
    mathutils.Vector((-1.0, 0.0, 0.0)),
    mathutils.Vector((-0.5, -1.0, 0.0)),
    mathutils.Vector((0.5, -1.0, 0.0)),
]
frag_mesh_data = bpy.data.meshes.new("frag_mesh_data")
frag_mesh_data.from_pydata(fragment_verts, [], [(0,1,2,3,4)])
frag_mesh_data.update()
frag_obj = bpy.data.objects.new("frag", frag_mesh_data)
frag_obj.location += mathutils.Vector((0,0,4))
scene.objects.link(frag_obj)'''

#for f in range(0,250):
    #cube.move()
    #cube.obj.keyframe_insert(data_path="location", frame=f)
    