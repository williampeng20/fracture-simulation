import bpy
import mathutils

# To start from console for console output:
# "/Applications/Blender/blender.app/Contents/MacOS/blender"

def clear_scene():
    for obj in bpy.data.objects:
        obj.select = True
    bpy.ops.object.delete()

clear_scene()

bpy.ops.mesh.primitive_cube_add(location=(0,0,10))

bpy.ops.mesh.primitive_plane_add()
bpy.ops.transform.resize(value=(10,10,10))

cube = bpy.data.objects['Cube']

v = mathutils.Vector((0,0,0))
a = mathutils.Vector((0,0,-9.8))
for f in range(0,100):
    t = f / bpy.data.scenes['Scene'].render.fps
    cube.location += v*t + 0.5*a*t**2
    cube.keyframe_insert(data_path="location", frame=f)
    


'''ob = bpy.context.object
mesh = ob.data

print(mesh.vertices[0].co)'''