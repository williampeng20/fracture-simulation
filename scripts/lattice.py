import bpy
import mathutils
import math
import random

PARTICLE_MASS = 0
BREAKING_THRESHOLD = 0.2
BOND_STRENGTH = 1.0
DAMPENING_FACTOR = 10.0
ALPHA = 0.5
TAU = 0.4
COLLISION_THRESHOLD = 0.0001
MAX_ADAPTIVE_LOOPS = 20

SCENE_OBJECTS = []

def is_inside(p, obj):
    # max_dist = 1.84467e+19
    result, point, normal, face = obj.closest_point_on_mesh(p)
    p2 = point-p
    v = p2.dot(normal)
    #print(v)
    return not(v < 0.0)

class CubeParticle:

    def __init__(self, parent, location):
        self.parent = parent
        self.center = location
        self.connections = {}
        self.mass = PARTICLE_MASS
        #normal corresponds to upwards direction of cube

    def add_neighbors(self, neighbors):
        for neighbor in neighbors:
            self.connections[neighbor] = BOND_STRENGTH


class SceneObject:

    #Collection of Particles

    def __init__(self, resolution, vertices, render_object):
        self.v = mathutils.Vector((0,0,0))
        self.rv = mathutils.Vector((0,0,0))
        self.a = mathutils.Vector((0,0,-9.8))
        self.ra = mathutils.Vector((0,0,0))

        self.obj = render_object

        self.vertices = vertices

        self.normal = mathutils.Vector((0,0,1))

        self.particles = []
        # Assuming 3D quadrilateral geometry
        x = [co[0] for co in vertices]
        y = [co[1] for co in vertices]
        z = [co[2] for co in vertices]
        width = max(x) - min(x)
        length = max(y) - min(y)
        height = max(z) - min(z)
        self.generate_mesh(width, length, height, min(x), min(y), min(z), resolution)
        self.res = resolution
        self.mass = sum([p.mass for p in self.particles])

        self.width = 0
        self.height = 0
        self.depth = 0

    #for now we assume that the cube is perfectly xyz-axis-aligned
    def generate_mesh(self, width, length, height, x, y, z, res):
        # width = x-axis, length = y-axis, height = z-axis
        particles = [[[CubeParticle(self, mathutils.Vector((x + i*res, y + j*res, z + k*res))) for k in range(0, math.ceil(height / res))] for j in range(0, math.ceil(length / res))] for i in range(0, math.ceil(width / res))]
        for i in range(0,len(particles)):
            for j in range(0,len(particles[0])):
                for k in range(0,len(particles[0][0])):
                    neighbors = []
                    if i-1 in range(0, len(particles)):
                        neighbors.append(particles[i-1][j][k])
                    if i+1 in range(0, len(particles)):
                        neighbors.append(particles[i+1][j][k])
                    if j-1 in range(0, len(particles[0])):
                        neighbors.append(particles[i][j-1][k])
                    if j+1 in range(0, len(particles[0])):
                        neighbors.append(particles[i][j+1][k])
                    if k-1 in range(0, len(particles[0][0])):
                        neighbors.append(particles[i][j][k-1])
                    if k+1 in range(0, len(particles[0][0])):
                        neighbors.append(particles[i][j][k+1])
                    particles[i][j][k].add_neighbors(neighbors)
                    self.particles.append(particles[i][j][k])

    # Move self to next time frame
    def move(self):
        t = 1 / bpy.data.scenes['Scene'].render.fps
        disp = self.v*t + 0.5*self.a*t**2
        self.v += self.a*t

        for vertex in self.vertices:
            vertex += disp
        # Check for collision
        for obj in SCENE_OBJECTS:
            if obj is not self and self.detect_collision(obj):
                correction = 0.0
                for vertex in self.vertices:
                    if is_inside(vertex, obj.obj):
                        adaptive_disp = disp*0.5
                        collision_point = vertex - adaptive_disp
                        count = MAX_ADAPTIVE_LOOPS
                        while count > 0:
                            _, point, normal, _ = obj.obj.closest_point_on_mesh(collision_point)
                            p2 = point-collision_point
                            v = p2.dot(normal)
                            print(v)
                            if abs(v) < COLLISION_THRESHOLD:
                                break
                            else:
                                adaptive_disp *= 0.5
                                collision_point += adaptive_disp * (1 if v < 0.0 else -1)
                                count -= 1
                        print("done")
                        disp_vector = collision_point - vertex
                        d = disp_vector.dot(disp_vector)**0.5 / disp.dot(disp)**0.5
                        correction = max(correction, d)
                        assert correction < 1
                for vertex in self.vertices:
                    vertex -= correction*disp
                disp *= 1 - correction
                self.v = mathutils.Vector((0,0,0))

        self.obj.location += disp
        for particle in self.particles:
            particle.center += disp

    #Input: List of adjacent particles to split from this scene object
    #Modifies self to account for splitting objects
    #Output: Reference to splitted object
    def split(self, particles):
        return

    # Input: SceneObject to detect collision with self
    # Output: True if collision, False if no collision
    def detect_collision(self, other_obj):
        for vertex in self.vertices:
            if is_inside(vertex, other_obj.obj):
                return True
        return False
