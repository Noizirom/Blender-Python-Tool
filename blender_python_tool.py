import bpy, bmesh, numpy as np, mathutils, math, os, json as js
from copy import deepcopy as dc


###############################################################
###Mesh Operations
###############################################################
#get selected verts, edges, and faces information
def get_sel():
    '''
    gs[0] = vert ##[0] = co, [1] = index, [2] = uv, [3] = normal, [4] = undeformed_co
    gs[1] = edge ##[0] = index, [1] = vert indexes
    gs[2] = face ##[0] = area, [1] = index, [2] = normal, [3] = center, [4] = vert indexes
    gs[3] = new from selected ##[0] = new vert index, [1] = new edge vert indexes,
    ##[2] = new face vert indexes, [3] = vert count, [4] = edge count, [5] = face count
    '''
    #vert, edge, face, new = get_sel()
    mode = bpy.context.active_object.mode
    vt = bpy.context.object.data.vertices
    ed = bpy.context.object.data.edges
    fa = bpy.context.object.data.polygons
    bpy.ops.object.mode_set(mode='OBJECT')
    #vertices
    countv = len(vt)
    selv = np.empty(countv, dtype=np.bool)
    vt.foreach_get('select', selv)
    co = np.empty(countv * 3, dtype=np.float32)
    vt.foreach_get('co', co)
    co.shape = (countv, 3)
    vidx = np.empty(countv, dtype=np.int32)
    vt.foreach_get('index', vidx)
    vnorm = np.empty(countv * 3, dtype=np.float32)
    vt.foreach_get('normal', vnorm)
    vnorm.shape = (countv, 3)
    und_co = np.empty(countv * 3, dtype=np.float32)
    vt.foreach_get('undeformed_co', und_co)
    und_co.shape = (countv, 3)
    uv_dict = dc({loop.vertex_index: bpy.context.object.data.uv_layers.active.data[loop.index].uv for loop in bpy.context.object.data.loops})
    uv_co = np.array([uv_dict[i] for i in vidx[selv]])
    #edges
    counte = len(ed)
    sele = np.empty(counte, dtype=np.bool)
    ed.foreach_get('select', sele)
    eidx = np.empty(counte, dtype=np.int32)
    ed.foreach_get('index', eidx)
    edg = np.array([i.vertices[:] for i in ed])
    #faces
    countf = len(fa)
    selfa = np.empty(countf, dtype=np.bool)
    fa.foreach_get('select', selfa)
    farea = np.empty(countf, dtype=np.float32)
    fa.foreach_get('area', farea)
    fidx = np.empty(countf, dtype=np.int32)
    fa.foreach_get('index', fidx)
    fnorm = np.empty(countf * 3, dtype=np.float32)
    fa.foreach_get('normal', fnorm)
    fnorm.shape = (countf, 3)
    fcnt = np.empty(countf * 3, dtype=np.float32)
    fa.foreach_get('center', fcnt)
    fcnt.shape = (countf, 3)
    fac = np.array([i.vertices[:] for i in fa])
    #New indexes
    v_count = len(vidx[selv])
    e_count = len(eidx[sele])
    f_count = len(fidx[selfa])
    new_idx = [i for i in range(v_count)]
    nv_Dict = {o: n for n, o in enumerate(vidx[selv].tolist())}
    new_e = [[nv_Dict[i] for i in nest] for nest in edg[sele]]
    new_f = [[nv_Dict[i] for i in nest] for nest in fac[selfa]]
    return dc([[co[selv], vidx[selv], uv_co, vnorm[selv], und_co[selv]], [eidx[sele], edg[sele]], [farea[selfa], fidx[selfa], fnorm[selfa], fcnt[selfa], fac[selfa]], [new_idx, new_e, new_f, v_count, e_count, f_count]])

#get a dictionary of {index: [coordinate]}
def get_co_idx(coords):
    mode = bpy.context.active_object.mode
    vt = bpy.context.object.data.vertices
    bpy.ops.object.mode_set(mode='OBJECT')
    #vertices
    coords = np.array(coords)
    coords = coords.tolist()
    countv = len(vt)
    co = np.empty(countv * 3, dtype=np.float32)
    vt.foreach_get('co', co)
    co.shape = (countv, 3)
    co = co.tolist()
    co_Dict = {idx: coor for idx, coor in enumerate(co)}
    return dc([keyfinder(co_Dict, i) for i in coords])

#get coordinates for greater than / less than for n axis
def get_gt_lt(colist, axis, var):
    if axis == 'x':
        axis = 0
    elif axis == 'y':
        axis = 1
    elif axis == 'z':
        axis = 2
    np.array(colist)
    colist.tolist()
    gt = [c for c in colist if c[axis] >= var]
    gti = get_co_idx(gt)
    gti = np.array(gti)
    lt = [c for c in colist if c[axis] <= var]
    lti = get_co_idx(lt)
    lti = np.array(lti)
    return(gti, lti)

#intersection of arrays
def array_intersect(a, b):
    #items in a that match in b
    return np.intersect1d(a, b)

#get range of verts in between 2 points
def vert_limits(colist, axis, hi, lo):
    up = get_gt_lt(colist, axis, hi)
    dn = get_gt_lt(colist, axis, lo)
    return np.intersect1d(up[1], dn[0])

#find nearest points
def v_nearest_range(co, dist):
    #return [location, normal, index, distance]
    bvht = mathutils.bvhtree.BVHTree()
    o_bvh = bvht.FromObject(bpy.context.object, bpy.context.scene)
    nearest_range = o_bvh.find_nearest_range(co, dist)
    return np.array(dc(nearest_range))

#find nearest point
def v_nearest(co, dist):
    #return [location, normal, index, distance]
    bvht = mathutils.bvhtree.BVHTree()
    o_bvh = bvht.FromObject(bpy.context.object, bpy.context.scene)
    nearest = o_bvh.find_nearest(co, dist)
    return np.array(dc(nearest))

#select vertices from a list
def sel_co(coords):
    coord = get_co_idx(coords)
    bpy.ops.object.mode_set(mode='OBJECT')
    for v in coord:
        bpy.context.object.data.vertices[v].select = True
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type="VERT")

#select vert
def vert_sel(idx):
    assert isinstance(idx, (int, list, tuple)), "Index must be int, list, or tuple."
    setMode()
    if isinstance(idx, (list, tuple)):
        for i in idx:
            obj[Name].data.vertices[i].select = True
    else:
        obj[Name].data.vertices[idx].select = True
    setMode('EDIT')

#all vert selection
def vert_sel_all():
    bpy.ops.object.mode_set(mode="OBJECT")
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.object.mode_set(mode="EDIT")

#all selected object vert deselection
def vert_sel_none():
    bpy.ops.object.mode_set(mode="OBJECT")
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode="EDIT")

#selected object vert delete
def vert_del_sel():
    bpy.ops.object.mode_set(mode="OBJECT")
    bpy.ops.mesh.delete(type='VERT')
    bpy.ops.object.mode_set(mode="EDIT")

#set smooth shading
def set_smooth(bool):
    bpy.ops.object.mode_set(mode='OBJECT')
    smooth = [bool] * len(bpy.context.object.data.polygons)
    bpy.context.object.data.polygons.foreach_set("use_smooth", smooth)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.context.scene.update()


###############################################################
###Vertex Group Operations
###############################################################
#vertex group index list
def vg_idx_list(vgn):
    return np.array([v.index for v in bpy.context.object.data.vertices if bpy.context.object.vertex_groups[vgn].index in [vg.group for vg in v.groups]])

#vertex group {name: [indexes]} dictionary
def vg_idx_dict():
    vn = [v.name for v in bpy.context.object.vertex_groups[:]]
    vd = {n: vg_idx_list(n) for n in vn}
    return dc(vd)

#vertex group weights {name: {indexes: weight}} dictionary
def vg_idx_weights():
    vid = vg_idx_dict()
    return dc({i: {v : bpy.context.object.vertex_groups[i].weight(v) for v in vid[i]} for i in vid})

#Create vertex group
def new_vg(newGroup, idx, weight=1.0, type="ADD"):
    #type = "ADD", "REPLACE", "REMOVE"
    nvg = bpy.context.object.vertex_groups.new(newGroup)
    nvg[newGroup].add(idx, weight, type)


###############################################################
#Object Operations
###############################################################
#Set mode of selected object
def setMode(m="OBJECT"):
    '''EDIT, OBJECT, POSE, SCULPT, VERTEX_PAINT,
    WEIGHT_PAINT, TEXTURE_PAINT, PARTICLE_EDIT,
    GPENCIL_EDIT '''
    bpy.ops.object.mode_set(mode=m)

#Object selection
def selectOb(Name):
    assert isinstance(Name, (str, list)), "Name must be string or list."
    if isinstance(Name, list):
        for i in Name:        
            bpy.data.objects[i].select = True
    else:
        bpy.data.objects[Name].select = True

#Object deselection
def deselectOb(Name):
    assert isinstance(Name, (str, list)), "Name must be string or list."
    if isinstance(Name, list):
        for i in Name:        
            bpy.data.objects[i].select = False
    else:
        bpy.data.objects[Name].select = False

#deselect all objects
def deselectAll():
    for o in bpy.data.objects[:]:
        deselectOb(o.name)

#Active Object
def activeOb(o):
    bpy.context.scene.objects.active = bpy.data.objects[o]

#Delete Object/s
def deleteOb():
    bpy.ops.object.delete(use_global=False)

#Create Object
def obj_mesh(co, faces):
    mesh = bpy.data.meshes.new("Obj")
    mesh.from_pydata(co, [], faces)
    mesh.validate()
    mesh.update(calc_edges = True) #calc_tessface = True
    Object = bpy.data.objects.new("Obj", mesh)
    Object.data = mesh
    bpy.context.scene.objects.link(Object)
    bpy.context.scene.objects.active = Object#
    Object.select = True

def obj_new(Name, co, faces):
    obj_mesh(co, faces)
    bpy.data.objects["Obj"].name = Name
    bpy.data.meshes[bpy.data.objects[Name].data.name].name = Name


#Create Edge    
def edge_mesh(verts, edges):
    mesh = bpy.data.meshes.new("Obj")
    mesh.from_pydata(verts, edges, [])
    mesh.validate()
    mesh.update()
    Object = bpy.data.objects.new("Obj", mesh)
    Object.data = mesh
    bpy.context.scene.objects.link(Object)
    bpy.context.scene.objects.active = Object#
    Object.select = True

def edge_new(Name, verts, edges):
    edge_mesh(verts, edges)
    bpy.data.objects["Obj"].name = Name
    bpy.data.meshes[bpy.data.objects[Name].data.name].name = Name

#Create Vert    
def vert_mesh(verts):
    mesh = bpy.data.meshes.new("Obj")
    mesh.from_pydata(verts, [], [])
    mesh.validate()
    mesh.update()
    Object = bpy.data.objects.new("Obj", mesh)
    Object.data = mesh
    bpy.context.scene.objects.link(Object)
    bpy.context.scene.objects.active = Object#
    Object.select = True

def vert_new(Name, verts):
    vert_mesh(verts)
    bpy.data.objects["Obj"].name = Name
    bpy.data.meshes[bpy.data.objects[Name].data.name].name = Name


###############################################################
###Collection Operations
###############################################################
#find dictionary key for value
def keyfinder(dictionary, val):
    for k, v in dictionary.items():
        if v == val:
            return k

#check if all items of collection are in other collection
def if_all(id, collection):
    Li = []
    for i in collection:
        if all(x in id for x in i) == True:
            Li.append(list(i))
    return dc(Li)


###############################################################
###Math Operations
###############################################################
#distance
def get_dist(a, b):
    #coordinate a and coordinate b
    a = np.array(a)
    b = np.array(b)
    return np.linalg.norm(a-b)

#angle
def get_angle(a, b, acute_bool):
    #coordinate a and coordinate b
    a = np.array(a)
    b = np.array(b)
    ang = np.arccos(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    if (acute_bool == True):
        return ang
    else:
        return 2 * np.pi - ang

#line median point
def get_median(a, b):
    #coordinate a and coordinate b
    a = list(a)
    b = list(b)
    vec = [[a], [b]]
    vec = np.array(vec)
    return np.median(vec, axis=0)

#radius from height and width
def radius_from(height, width):
    r = ((height)**2 + (width)**2)/2*(height)
    return(r)

#Quadratic Equations
###############################################################

def quad_split(a, b, c):
    '''
    #ax2 + bx +c = 0
    (x-x1)(x-x2) = 0
    x2 + 3x – 4 = (x + 4)(x – 1) = 0
    x = 1 and x = –4
    '''
    d = ((b**2) - (4 * a * c)) # discriminant
    #validate
    if d < 0:
        print ("This equation has no real solution")
    elif d == 0:
        x1 = (-b + math.sqrt(d)) / (2 * a)
        x2 = None
    else:
        x1 = (-b + math.sqrt(d)) / (2 * a)
        x2 = (-b - math.sqrt(d)) / (2 * a)
    return(x1, x2)

def quad(a, b, c, fr=0):
    '''
    #ax2 + bx +c = 0
    a = speed
    b = initial velocity
    c = initial height
    fr = obj.frame_current - obj.frame_start
    '''
    #quad expression fo z height
    obj_z = (a**2 * fr)/2 + (b * fr) + c
    return(obj_z)

#Trigonometry
###############################################################

def findAdjacent(opp, ang):
    #opp = opposite side
    #ang = angle #radians(ang)
    adj = (opp / (tan(ang)))
    return(adj)

def findOpposite(adj, ang):
    #adj = adjacent side
    #ang = angle #radians(ang)
    opp = (adj * (tan(ang)))
    return(opp)

def hyp(adj, opp):
    #hypotenuse
    #adj = adjacent side
    #opp = opposite side
    h = sqrt((adj**2) + (opp**2))
    return(h)

def ang_fr_2_sides(adj, opp):
    #angle from adjacent and opposite
    #ang = angle #radians(ang)
    ang = atan(opp / adj)
    return(ang)

def rem_rt_ang(ang):
    #remaining right angle in degrees
    #ang = angle #radians(ang)
    a = 180 - 90 - ang
    return(a)
    
#physics
###############################################################

def phys_time(dist, vel):
    #dist = distance total
    #vel = velocity
    t = (dist / vel)
    return(t) #time total #delta_x

def phys_dist_drop(time, grav=-9.81):
    #time = current frames #/ 30 #/ 60
    #grav = gravity #world gravity
    d = ((grav / 2) * (time**2))
    return(d) #delta_y for x

def traj(x, v0, ang, g=9.81):
    #trajectory
    #get y at current x
    '''
    y = vertical position (m)
    x = horizontal position (m)
    v0 = initial velocity (combined components, m/s)
    g = acceleration due to gravity (9.80 m/s2)
    ang = angle of the initial velocity from the horizontal plane (radians or degrees)
    '''
    y = ((x * tan(ang))- ((g * (x**2)) / ((2 * v0) * ((cos(ang))**2))))
    return(y) #delta_y for x



###############################################################
###File Operations
###############################################################
#write to json file
def write_js(fileName, path, data):
    fpath = path + fileName + '.json'
    with open(fpath, 'w') as f:
        js.dump(data, f)

#append to json file
def append_js(fileName, path, data):
    fpath = path + fileName + '.json'
    with open(fpath, 'a') as f:
        js.dump(data, f)


###############################################################
### Main ###
###############################################################

#v_, e_, f_, n_ = get_sel()

