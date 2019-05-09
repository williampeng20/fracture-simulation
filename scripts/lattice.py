import bpy
import mathutils
import math
import random
import copy

PARTICLE_MASS = 0
BREAKING_THRESHOLD = 0.2
BOND_STRENGTH = 1.0
DAMPENING_FACTOR = 10.0
ALPHA = 0.5
TAU = 0.4
COLLISION_THRESHOLD = 0.0001
MAX_ADAPTIVE_LOOPS = 20
BOUNCE_FACTOR = 0.2
FRICTION = 0.3

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

    def __init__(self, resolution, vertices, render_object, fracturable):
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
        self.fracturable = fracturable

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
            '''if obj.__class__.__name__ == 'glassPane':
                print("glass")
                if True:
                    print("collision")
                    obj.generate_fractures(mathutils.Vector((0,0,0.5)), 5)'''
            if obj is not self and self.detect_collision(obj):
                correction = 0.0
                normal = mathutils.Vector((0,0,0))
                for vertex in self.vertices:
                    if is_inside(vertex, obj.obj):
                        adaptive_disp = disp*0.5
                        collision_point = vertex - adaptive_disp
                        count = MAX_ADAPTIVE_LOOPS
                        while count > 0:
                            _, point, normal, _ = obj.obj.closest_point_on_mesh(collision_point)
                            p2 = point-collision_point
                            v = p2.dot(normal)
                            #print(v)
                            if abs(v) < COLLISION_THRESHOLD:
                                break
                            else:
                                adaptive_disp *= 0.5
                                collision_point += adaptive_disp * (1 if v < 0.0 else -1)
                                count -= 1
                        #print("done")
                        disp_vector = collision_point - vertex
                        d = disp_vector.dot(disp_vector)**0.5 / disp.dot(disp)**0.5
                        correction = max(correction, d)
                        assert correction < 1
                for vertex in self.vertices:
                    vertex -= correction*disp
                disp *= 1 - correction
                self.v += (self.v.length * normal * BOUNCE_FACTOR) - (self.v.dot(-normal) / (normal.length) * (-normal).normalized()) - FRICTION * (self.v - self.v.dot(-normal) / (normal.length) * (-normal).normalized() )

        self.obj.location += disp
        for particle in self.particles:
            particle.center += disp

    def cell_fracture(self, impact_points):
        return False

    #Input: List of adjacent particles to split from this scene object
    #Modifies self to account for splitting objects
    #Output: Reference to splitted object
    def split(self, particles):
        return False

    # Input: SceneObject to detect collision with self
    # Output: True if collision, False if no collision
    def detect_collision(self, other_obj):
        for vertex in self.vertices:
            if is_inside(vertex, other_obj.obj):
                return True
        return False

class glassPane: 

    #coordinates in object space
    def __init__(self, bot_left, top_right, obj):
        self.left = bot_left.x
        self.right = top_right.x
        self.top = top_right.y
        self.bottom = bot_left.y

        self.edges = [
            (mathutils.Vector((self.left, self.bot, 0.0)), mathutils.Vector((self.left, self.top, 0.0))),
            (mathutils.Vector((self.left, self.bot, 0.0)), mathutils.Vector((self.right, self.bot, 0.0))),
            (mathutils.Vector((self.left, self.top, 0.0)), mathutils.Vector((self.right, self.top, 0.0))),
            (mathutils.Vector((self.right, self.top, 0.0)), mathutils.Vector((self.right, self.bot, 0.0)))]

        self.obj = obj

    def sample_particle(self, p):

        while True:
            nx = random.uniform(-1, 1)
            ny = random.uniform(-1, 1)
            
            nx_bias = pow(nx, 3)
            ny_bias = pow(ny, 3)

            sx = p.x + nx_bias * abs(self.right - self.left)
            sy = p.y + ny_bias * abs(self.top - self.bot)

            if sx >= self.left and sx <= self.right and sy <= self.top and sy >= self.bottom:
                return mathutils.Vector((sx, sy, p.z))


    #p := impact point
    #k := total number of fracture reference particles
    def generate_fractures(self, p, k):
        pts = []
        for _ in range(k):
            s = self.sample_particle(p)
            pts += [s]
        hull = self.convex_hull(pts)

    def distance(self, u, v): 
        return ((u.x - v.x) ** 2 + (u.y - v.y) ** 2) ** 0.5
    
    def closest_point(self, p, points):
        min_dist = 100
        closest_v = None
        for v in points:
            if p != v:
                if self.distance(p, v) < min_dist:
                    min_dist = self.distance(p, v)
                    closest_v = v
        return closest_v
    
    def new_shard(self, vertices): 
        ordered_vertices = make_ccw(vertices)
        #make new frag
        fragment_verts = ordered_vertices
        frag_mesh_data = bpy.data.meshes.new("frag_mesh_data")
        frag_mesh_data.from_pydata(fragment_verts, [], [(0,1,2,3,4)])
        frag_mesh_data.update()
        frag_obj = bpy.data.objects.new("frag", frag_mesh_data)
        frag_obj.location += mathutils.Vector((0,0,4))
        scene.objects.link(frag_obj)
        
    def make_ccw(self, points):
        vertices = copy.deepcopy(points)
        start = vertices[0]
        vertices.remove(start)
        next_0 = self.closest_point(start, vertices)
        vertices.remove(next_0)
        next_1 = self.closest_point(start, vertices)
        vertices.remove(next_1)
        ordered = []
        if self.direction(start, next_0, next_1) < 0:
            # ordering is ccw
            ordered = [start, next_0, next_1]
        if self.direction(start, next_1, next_0) < 0:
            # ordering is ccw
            ordered = [start, next_1, next_0]

        if self.direction(next_0, start, next_1) < 0:
            # ordering is ccw
            ordered = [next_0, start, next_1]
        if self.direction(next_0, next_1, start) < 0:
            # ordering is ccw
            ordered = [next_0, next_1, start]

        if self.direction(next_1, start, next_0) < 0:
            # ordering is ccw
            ordered = [next_1, start, next_0]
        if self.direction(next_1, next_0, start) < 0:
            # ordering is ccw
            ordered = [next_1, next_0, start]

        while vertices:
            next_closest = self.closest_point(ordered[-1], vertices)
            ordered.append(next_closest)
            vertices.remove(next_closest)
        return ordered

        cvh = self.convex_hull(pts)
        outer = self.outer_points(self.edges, math.ceil(len(cvh) / len(self.edges)))
        while len(cvh) > 2:
            for j in range(len(cvh)):
                cur_frag = []
                cur_frag.append(cvh[j])
                crease, start = self.closest_point(cvh[j], outer)
                cur_frag.append(crease)
                x = j+1 if j+1 < len(cvh) else 0
                cur_frag.append(cvh[x])
                ptsd, end = self.closest_point(cvh[x], outer)
                if start != end:
                    cur_frag.append(ptsd)
                    cur_frag += outer[min(start, end)+1:max(start,end)]
                self.generate_frag(cur_frag)
            outer = cvh
            samples = [sample for sample in pts if sample not in cvh]
            cvh = self.convex_hull(samples)

    def generate_frag(self, points):
        #points = self.make_ccw(points)
        frag_mesh_data = bpy.data.meshes.new("frag_mesh_data")
        frag_mesh_data.from_pydata(points, [], [tuple(range(len(points)))])
        frag_mesh_data.update()
        frag_obj = bpy.data.objects.new("frag", frag_mesh_data)
        frag_obj.location = self.obj.location
        bpy.context.scene.objects.link(frag_obj)


    # k := number of random outer-edge points
    # edges := list of edges to use
    def outer_points(self, edges, k):
        #for _ in range(k):
        return [mathutils.Vector((self.left, self.bot, 0.0)), mathutils.Vector((self.left, self.top, 0.0)),
            mathutils.Vector((self.right, self.top, 0.0)), mathutils.Vector((self.right, self.bot, 0.0))]

    # p := point of reference
    # points := list of points to find closest to p
    def closest_point(self, p, points):
        closest = None
        dist = 999999
        index = -1
        for i in range(len(points)):
            point = points[i]
            length = (p - point).length
            if length < dist:
                dist = length
                closest = point
                index = i
        return closest, index

    #input: points = set of points
    #output: hull points
    def convex_hull(self, points):
        # find the leftmost point
        a =  min(points, key = lambda point: point.x)
        index = points.index(a)
        
        # selection sort
        l = index
        result = []
        result.append(a)
        while (True):
            q = (l + 1) % len(points)
            for i in range(len(points)):
                if i == l:
                    continue
                # find the greatest left turn
                # in case of collinearity, consider the farthest point
                d = self.direction(points[l], points[i], points[q])
                if d > 0:
                    q = i
            l = q
            if l == index:
                break
            result.append(points[q])

        return result

    # calculates the cross product of vector p1 and p2
    # if p1 is clockwise from p2 wrt origin then it returns +ve value
    # if p2 is anti-clockwise from p2 wrt origin then it returns -ve value
    # if p1 p2 and origin are collinear then it returs 0
    def cross_product(self, p1, p2):
        return p1.x * p2.y - p2.x * p1.y

    # returns the cross product of vector p1p3 and p1p2
    # if p1p3 is clockwise from p1p2 it returns +ve value
    # if p1p3 is anti-clockwise from p1p2 it returns -ve value
    # if p1 p2 and p3 are collinear it returns 0
    def direction(self, p1, p2, p3):
        return self.cross_product(p3 - p1, p2 -p1)