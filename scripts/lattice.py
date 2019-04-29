import bpy
import mathutils

class Cube:
    
    def __init__(self, parent, adj, location, vertices):
        self.parent_object = parent
        self.vertices = vertices
        self.center = location
        self.connections = {}
        for neighbor in adj:
            self.connections[neighbor] = 1.0
        
class TVertex:

    def __init__(self, location):
        self.location = location
