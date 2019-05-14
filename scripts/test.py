import math

class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        return Vector(self.x * other.x, self.y * other.y, self.z * other.z)

    def __repr__(self):
        return "({},{},{})".format(self.x, self.y, self.z)

    def __str__(self):
        return "({},{},{})".format(self.x, self.y, self.z)

def cross_product(p1, p2):
    return p1.x * p2.y - p2.x * p1.y

# returns the cross product of vector p1p3 and p1p2
# if p1p3 is clockwise from p1p2 it returns +ve value
# if p1p3 is anti-clockwise from p1p2 it returns -ve value
# if p1 p2 and p3 are collinear it returns 0
def direction(p1, p2, p3):
    return cross_product(p3 - p1, p2 -p1)


def convex_hull(points):
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
            d = direction(points[l], points[i], points[q])
            if d > 0:
                q = i
        l = q
        if l == index:
            break
        result.append(points[q])
    return result

# works but we need a tie breaker to bias towards the correct y values.
def make_ccw2(vertices):
    points = vertices[:]
    points.sort(key = lambda p : p.x)
    lm = points[0]
    rm = points[-1]

    cutoff_slope = (rm.y - lm.y) / (rm.x - lm.x)

    ordered = [lm]

    points.remove(lm)

    for i, p in enumerate(points):
        if p.y <= cutoff_slope * (p.x - lm.x) + lm.y:
            # j = i + 1
            # ties = []
            # while j < len(points) and points[j].x == p.x:
            #     ties.append(points[j])
            #
            ordered.append(p)
            points.remove(p)

    points.reverse()
    ordered.extend(points)
    return ordered

def test2():
    points1 = [Vector(-1,0,0),
        Vector(0,1,0),
        Vector(-4,-4,0),
        Vector(4,4,0)]
    result = make_ccw2(points1)
    print(result)

    points2 = [Vector(-4,4,0),
        Vector(0,1,0),
        Vector(-4,-4,0),
        Vector(-1,0,0)]
    result = make_ccw2(points2)
    print(result)

def test1():
    points = [Vector(1,2,0),
        Vector(2,6,0),
        Vector(3,5,0),
        Vector(4,4,0),
        Vector(5,5,0),
        Vector(5,7,0),
        Vector(6,3,0)]
    result = convex_hull(points)
    print(result)
    #for v in result:
    #    print("(", v.x, ", ", v.y, ")")
