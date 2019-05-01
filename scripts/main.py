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
    mathutils.Vector((1.0, 1.0, 9.0)),
    mathutils.Vector((1.0, -1.0, 9.0)),
    mathutils.Vector((-1.0, -1.0, 9.0)),
    mathutils.Vector((-1.0, 1.0, 9.0)),
    mathutils.Vector((1.0, 1.0, 11.0)),
    mathutils.Vector((1.0, -1.0, 11.0)),
    mathutils.Vector((-1.0, -1.0, 11.0)),
    mathutils.Vector((-1.0, 1.0, 11.0))]
cube_faces = [
    (0, 1, 2, 3),
    (4, 7, 6, 5),
    (0, 4, 5, 1),
    (1, 5, 6, 2),
    (2, 6, 7, 3),
    (4, 0, 3, 7)]

cube_mesh_data = bpy.data.meshes.new("cube_mesh_data")
cube_mesh_data.from_pydata(cube_vertices, [], cube_faces)
cube_mesh_data.update()
cube_obj = bpy.data.objects.new("Cube", cube_mesh_data)
scene.objects.link(cube_obj)

cube = lattice.SceneObject(0.5, cube_vertices, cube_obj)
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

'''plane_mesh_data = bpy.data.meshes.new("plane_mesh_data")
plane_mesh_data.from_pydata(plane_vertices, plane_edges, plane_faces)
plane_mesh_data.update()
plane_obj = bpy.data.objects.new("Plane", plane_mesh_data)
scene.objects.link(plane_obj)'''

plane = lattice.SceneObject(0.5, plane_vertices, bpy.data.objects['Plane'])
lattice.SCENE_OBJECTS.append(plane)

#v = mathutils.Vector((0,0,0))
#a = mathutils.Vector((0,0,-9.8))
for f in range(0,250):
    cube.move()
    cube.obj.keyframe_insert(data_path="location", frame=f)

    #if cube.detect_collision(plane):
        #break
    


'''ob = bpy.context.object
mesh = ob.data

print(mesh.vertices[0].co)'''