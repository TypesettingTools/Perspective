import sys
import argparse
import math
import operator
import re

class Point:
    def __init__(self, x=0, y=0, z=0):
        if type(x) == type(tuple()) or x is list:
            self.x = x[0]
            self.y = x[1]
            self.z = x[2] if len(x)>2 else z
        else:
            self.x = x
            self.y = y
            self.z = z
    
    def __repr__(self):
        return str((self.x, self.y, self.z) if abs(self.z)>1E5 else (self.x, self.y))
        
    def length(self):
        return math.sqrt(self.x**2+self.y**2+self.z**2)
    
    def add(self, p):
        return Point(self.x + p.x, self.y + p.y, self.z + p.z)
    
    def sub(self, p):
        return Point(self.x - p.x, self.y - p.y, self.z - p.z)
    
    def rot_y(self, a):
        rot_v = Point(self.x*math.cos(a) - self.z * math.sin(a),
                      self.y,
                      self.x*math.sin(a) + self.z * math.cos(a))
        return rot_v

    def rot_x(self, a):
        rot_v = Point(self.x, 
                        self.y*math.cos(a) + self.z * math.sin(a),
                       -self.y*math.sin(a) + self.z * math.cos(a))
        return rot_v
    
    def rot_z(self, a):
        rot_v = Point(self.x*math.cos(a) + self.y * math.sin(a),
                       -self.x*math.sin(a) + self.y * math.cos(a),
                        self.z)
        return rot_v
    
    def mul(self, m):
        return Point(self.x * m, self.y * m, self.z * m)

def vector(a, b):
    return Point(b.x-a.x, b.y-a.y, b.z-a.z)

def vec_pr(a, b):
    return Point(a.y*b.z - a.z*b.y, - a.x*b.z + a.z*b.x, a.x*b.y - a.y*b.x)

def sc_pr(a, b):
    return a.x*b.x + a.y*b.y + a.z*b.z

def dist(a, b):
    return a.sub(b).length()

def intersect(l1, l2):
    if vec_pr(vector(*l1), vector(*l2)).length()==0:
        return l1[0].add(vector(*l1).mul(1E30))
    else:
        d = ((l1[0].x - l1[1].x)*(l2[0].y-l2[1].y)
            -(l1[0].y - l1[1].y)*(l2[0].x-l2[1].x))
        x = (vec_pr(*l1).z*(l2[0].x-l2[1].x)
           - vec_pr(*l2).z*(l1[0].x-l1[1].x))
        y = (vec_pr(*l1).z*(l2[0].y-l2[1].y)
           - vec_pr(*l2).z*(l1[0].y-l1[1].y))
        x/=d
        y/=d
        return Point(x, y)


def unrot(coord_in, org, diag=True, get_rot=False):
    screen_z = 312.5
    shift = org.mul(-1)
    coord = [c.add(shift) for c in coord_in]
    center = intersect(coord[0::2], coord[1::2])
    center = Point(center.x, center.y, screen_z)
    rays = [Point(c.x, c.y, screen_z) for c in coord]
    f = []
    for i in xrange(2):
        vp1 = vec_pr(rays[0+i], center).length()
        vp2 = vec_pr(rays[2+i], center).length()
        a = rays[0+i]
        c = rays[2+i].mul(vp1/vp2)
        m = a.add(c).mul(0.5)
        r = center.z/m.z;
        a = a.mul(r)
        c = c.mul(r)
        f.append((a, c))
    (a, c), (b, d) = f[0], f[1]
    ratio = abs(dist(a, b)/dist(a, d))
    diag_diff = ((dist(a, c)-dist(b, d)))/(dist(a, c)+dist(b, d))
    n = vec_pr(vector(a, b), vector(a, c))
    n0 = vec_pr(vector(rays[0], rays[1]), vector(rays[0], rays[2]))
    flip = 1 if sc_pr(n, n0)>0 else -1
    if not get_rot:
        return diag_diff if diag else ratio#*flip
    if flip<0:
        return None
    fry = math.atan(n.x/n.z)
    s = ""
    s+= "\\fry%.2f" % (-fry/math.pi*180)
    rot_n = n.rot_y(fry)
    frx = -math.atan(rot_n.y/rot_n.z)
    if n0.z < 0:
        frx += math.pi
    s+= "\\frx%.2f" % (-frx/math.pi*180)
    n = vector(a, b)
    ab_unrot = vector(a, b).rot_y(fry).rot_x(frx)
    ac_unrot = vector(a, c).rot_y(fry).rot_x(frx)
    ad_unrot = vector(a, d).rot_y(fry).rot_x(frx)
    frz = math.atan2(ab_unrot.y, ab_unrot.x)
    s += "\\frz%.2f" % (-frz/math.pi*180)
    ad_unrot = ad_unrot.rot_z(frz)
    fax = ad_unrot.x/ad_unrot.y
    if abs(fax)>0.01:
        s += "\\fax%.2f" % (fax)
    return s

def binary_search(f, l, r, eps):
    fl = f(l)
    fr = f(r)
    if fl<=0 and fr>=0:
        op = operator.gt
    elif fl>=0 and fr<=0:
        op = operator.lt
    else:
        return None
    l, r = map(float, (l, r))
    while (r - l > eps):
        c = (l + r) / 2
        if op(f(c), 0):
            r = c
        else:
            l = c
    return (l + r) / 2

def find_ex(f):
    w_center = [0, 0]
    w_size = 100000.0
    iterations = int(math.log(w_size*100, 4))
    s = 4
    for k in xrange(iterations):
        res = []
        for i in xrange(-s, s):
            x = w_center[0] + w_size*i/10
            for j in xrange(-s, s):
                y = w_center[1] + w_size*j/10
                res.append((unrot(coord, Point(x, y)), x, y))
        ex = f(res)
        w_center = [ex[1], ex[2]]
        w_size/=3
    return Point(ex[1], ex[2])

def parse_args(args):
    dsc = "Script finds tags for .ass subtitles that can map a rectangle to a given tetragon"
    parser = argparse.ArgumentParser(description=dsc)

    parser.add_argument(type=str, dest="coord", metavar="<Coordinates>",
                        help="Coordinates of tetragon")
    parser.add_argument("-o", "--origin", default=False,
                        type=str, metavar="O", dest="target_origin",
                        help="Desired origin coordinates. Coordinates for multiple origin points should be separated by semicolon.")
    parser.add_argument("-r", "--ratio", default=False,
                        type=str, metavar="R", dest="target_ratio",
                        help="Output file")
    parser.add_argument("-s", "--scale", default=1,
                        type=float, metavar="S", dest="scale",
                        help="Scaling factor. Simply multiples input coordinates.")
    return parser.parse_args(args)

def getfloats(s):
    return list(map(float, re.findall('[\d\.]+', s)))

args = parse_args(sys.argv[1:])
coord = getfloats(args.coord)[:8]
coord = zip(coord[::2], coord[1::2])
coord = list(map(Point, coord))
target_org = args.target_origin
scale = args.scale
if args.target_origin:
    target_orgs = []
    for o in args.target_origin.split(";"):
        target_orgs.append(Point(*getfloats(o)[:2]))
else:
    target_orgs = False
if scale != 1:
    for i in xrange(len(coord)):
        coord[i] = coord[i].mul(scale)
    if target_orgs:
        for i in xrange(len(target_orgs)):
            target_orgs[i] = target_orgs[i].mul(scale)
if args.target_ratio:
    target_ratio = getfloats(args.target_ratio)[0]
else:
    target_ratio = dist(coord[0], coord[1])/dist(coord[0], coord[3])

if target_orgs:
    print("\nTransform for target org:")
    for target_org in target_orgs:
        tf_tags = unrot(coord, target_org, get_rot=True)
        if tf_tags is None:
            print(tf_tags)
        else:
            print("\\org(%.1f, %.1f)\n" % (target_org.x, target_org.y) + tf_tags)

mn_point = find_ex(min)
mx_point = find_ex(max)
c = mn_point.add(mx_point).mul(0.5)
v = mn_point.sub(mx_point).rot_z(math.pi/2).mul(100000)
inf_p = c.add(v)
if unrot(coord, inf_p)>0:
    mn_center = True
    center = mn_point
    other = mx_point
else:
    mn_center = False
    center = mx_point
    other = mn_point
v = other.sub(center)

def zero_on_ray(center, v, a, eps):
    vrot = v.rot_z(a)
    def f(x):
        p=vrot.mul(x).add(center)
        return unrot(coord, p, diag=True)
    l = binary_search(f, 0, (center.length()+1000000)/v.length(), eps)
    if l is None:
        return None
    p = vrot.mul(l).add(center)
    ratio = unrot(coord, p, diag=False)
    r = unrot(coord, p, get_rot=False)
    if r is None:
        return None
    else:
        return p, ratio

rots = []
steps = 100
for i in xrange(steps):
    a = 2*math.pi*i/steps
    zero = zero_on_ray(center, v, a, 1E-02)
    if zero is None:
        continue
    p, ratio = zero
    rots.append((ratio, p, a))

if len(rots)==0:
    print("\nNo proper perspective found.")
    exit()
print("\nTransforms near center of tetragon:")
t_center = coord[0].add(coord[1]).add(coord[2]).add(coord[3]).mul(0.25)
ratio, p, a = min(rots, key = lambda x:dist(t_center, x[1]))
tf_tags = unrot(coord, p, get_rot=True)
if tf_tags is None:
    print(tf_tags)
else:
    print("Ratio=%f" % ratio)
    print("\\org(%.1f, %.1f)" % (p.x, p.y) + tf_tags)
segs = []
for i in xrange(len(rots)):
    if (rots[i-1][0]-target_ratio)*(rots[i][0]-target_ratio)<=0:
        segs.append((rots[i-1][2], rots[i][2]))
orgs = []
got_tf = False
if len(segs)>0:
    print("\nTransforms with target ratio:")
    for seg in segs:
        def f(a):
            res = zero_on_ray(center, v, a, 1E-05)
            if res is None:
                return 1E7
            else:
                p, ratio = res
                return ratio-target_ratio
        a = binary_search(f, seg[0], seg[1], eps=1E-04)
        if a is None:
            a = seg[0]
        p, ratio = zero_on_ray(center, v, a, 1E-05)
        tf_tags = unrot(coord, p, get_rot=True)
        if tf_tags is None: continue
        print("Ratio=%f" % ratio)
        print("\\org(%.1f, %.1f)" % (p.x, p.y) + tf_tags)
        got_tf = True
if not got_tf and len(rots)>0:
    print("\nTransforms close to target ratio:")
    ratio, p, a = min(rots, key = lambda x:abs(target_ratio - x[0]))
    tf_tags = unrot(coord, p, get_rot=True)
    if tf_tags is not None:
        print("Ratio=%f" % ratio)
        print("\\org(%.1f, %.1f)" % (p.x, p.y) + tf_tags)
