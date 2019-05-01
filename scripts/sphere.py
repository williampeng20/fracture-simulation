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

bpy.ops.mesh.primitive_uv_sphere_add(location=(-18,-18,10))
sphere_obj = bpy.data.objects['Sphere']
sphere_vertices = [vertex.co+mathutils.Vector((-18,-18,10)) for vertex in sphere_obj.data.vertices]

sphere = lattice.SceneObject(0.5, sphere_vertices, sphere_obj)
lattice.SCENE_OBJECTS.append(sphere)


bpy.ops.mesh.primitive_cube_add()
for obj in bpy.context.selected_objects:
    obj.name = "Plane"
plane_vertices = [
    mathutils.Vector((16.0, 16.0, 0.0)),
    mathutils.Vector((16.0, -16.0, 0.0)),
    mathutils.Vector((-16.0, -16.0, 0.0)),
    mathutils.Vector((-16.0, 16.0, 0.0)),
    mathutils.Vector((16.0, 16.0, 0.5)),
    mathutils.Vector((16.0, -16.0, 0.5)),
    mathutils.Vector((-16.0, -16.0, 0.5)),
    mathutils.Vector((-16.0, 16.0, 0.5))]
for i in range(0, len(plane_vertices)):
    bpy.data.objects['Plane'].data.vertices[i].co.x *= 16
    bpy.data.objects['Plane'].data.vertices[i].co.y *= 16
    bpy.data.objects['Plane'].data.vertices[i].co.z = bpy.data.objects['Plane'].data.vertices[i].co.z / 4 + 0.25

plane = lattice.SceneObject(0.5, plane_vertices, bpy.data.objects['Plane'])
lattice.SCENE_OBJECTS.append(plane)

#v = mathutils.Vector((0,0,0))
sphere.v += mathutils.Vector((8,8,0))
for f in range(0,100):
    print("frame",f)
    sphere.move()
    sphere.obj.keyframe_insert(data_path="location", frame=f)

    #if cube.detect_collision(plane):
        #break
    


'''ob = bpy.context.object
mesh = ob.data

print(mesh.vertices[0].co)'''