import bpy, numpy as np
from math import*


def dome_ops(rad, seg, rings):
    count = seg * rings
    iseg = seg-1
    icount = count-1
    last_ring = count-seg
    verts, edges = [], []
    s_ang = 360 / seg
    r_ang = 90 / rings
    idx = [i for i in range(seg) if i % 2 == 0]
    for r in range(rings):
        r_theta = r * r_ang
        r_theta = radians(r_theta)
        r_rad = cos(r_theta) * rad
        step = r * seg
        e_list = [[i + step, (i + 1) + step] for i in range((seg - 1))]
        e_end = [e_list[-1][1], e_list[0][0]]
        e_copy = e_list.copy()
        e_copy.append(e_end.copy())
        if r != 0 and r % 2 != 0:
            e_list2 = e_list[::-1].copy()
            e_list2 = [[i[1], i[0]] for i in e_list2]
            e_end = [e_list2[-1][1], e_list2[0][0]]
            edges.append(e_end)
            edges = edges + e_list2
        else:
            edges = edges + e_list
            edges.append(e_end)
        if r > 0:
            sides = []
            for i in e_copy:
                side = [[i[0], i[0] - seg], [i[1] - seg, i[1]]]
                sides.append(side)
            for s in idx:
                edges = edges + sides[s]
        if r == (rings - 1):
            e_top =[[[count, i[0]], [i[1], count]] for i in e_copy]
            i_list = []
            for i in idx:
                i_list = i_list + e_top[i]
            edges = edges + i_list
        for s in range(seg):
            s_theta = s * s_ang
            s_theta = radians(s_theta)
            z = float(round(((sin(r_theta)) * rad), ndigits=5))
            co = [float(round((cos(s_theta) * (rad*cos(r_theta))), ndigits=5)),
                float(round((sin(s_theta) * (rad*cos(r_theta))), ndigits=5)), z]
            verts.append(co)
    if rings > 1:
        verts.append([0.,0., rad])
    else:
        verts.append([0.,0.,0.])
    f_close = []
    for r in range(rings-1):
        segr = seg*r
        f_close.append([iseg+segr, 0+segr, seg+segr, seg+iseg+segr])
    f_top_idx = [i for i in range((last_ring),count)]
    ft = [[i, i+1, count] for i in f_top_idx if i != (icount)]
    f_top = ft + [[icount,last_ring,count]]
    f_sides = [[i, i + 1, i + 1 + seg, i + seg] for i in range(last_ring) if i not in [(r*seg-1) for r in range(rings) if r!=0]]
    faces = f_sides + f_close + f_top
    return verts, edges, faces



#Create Object
def obj_mesh(co, edges, faces):
    mesh = bpy.data.meshes.new("Obj")
    mesh.from_pydata(co, edges, faces)
    mesh.validate()
    mesh.update()
    Object = bpy.data.objects.new("Obj", mesh)
    Object.data = mesh
    bpy.context.collection.objects.link(Object)
    bpy.context.view_layer.objects.active = Object
    Object.select_set(True)

def obj_new(Name, co, edges, faces):
    obj_mesh(co, edges, faces)
    bpy.data.objects["Obj"].name = Name
    bpy.data.meshes[bpy.data.objects[Name].data.name].name = Name



def create_dome(Name, rad, seg, rings):
    co, edges, faces = dome_ops(rad, seg, rings)
    obj_new(Name, co, edges, faces)
    print(Name, co, edges, faces)


def roty_matrix(theta):
    ry = np.array(
                    [ [np.cos(theta), 0, np.sin(theta)],
                    [0, 1, 0],
                    [-np.sin(theta), 0, np.cos(theta)] ]
                    )
    return ry

def create_sphere(Name, rad, seg, rings):
    co, edges, faces = dome_ops(rad, seg, rings)
    obj_new(Name+'2', co, edges, faces)
    mirror = roty_matrix(radians(180))
    co2 = [(v @ mirror) for v in co]
    obj_new(Name, co2, edges, faces)
    bpy.ops.object.join()
    bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='MEDIAN')

def create_cylinder(Name, rad, seg, rings, length):
    if rings == 1:
        step = 1
    else:
        step = length / (rings-1)
    count = seg*rings
    co, e_none, f_one = dome_ops(rad, seg, 1)
    co_none, edges, faces = dome_ops(rad, seg, rings)
    r1 = co[:-1]
    edges = edges[:-seg]
    verts = [[c[0],c[1],c[2]+(step*r)] for r in [j for j in range(rings+2)] for c in co[:-1]]
    f_total = faces[:-seg]
    if rings == 1:
        verts = r1+[[0,0,0]]
        f_total = [[i,i+1] for i in range(seg-1)]+[[seg-1,0 ]]
    obj_new(Name, verts, edges, f_total)
    
def create_capsule(Name, rad, seg, rings, length, sections):
    obj = bpy.data.objects
    new_length = (length / 2) - rad
    co, edges, faces = dome_ops(rad, seg, rings)
    obj_new(Name+'2', co, edges, faces)
    mirror = roty_matrix(radians(180))
    co2 = [(v @ mirror) for v in co]
    obj_new(Name+'3', co2, edges, faces)
    create_cylinder(Name, rad, seg, sections, new_length*2)
    obj[Name+'2'].location.z = new_length
    obj[Name+'3'].location.z = -new_length
    obj[Name].location.z = -new_length
    bpy.ops.object.join()
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS', center='MEDIAN')    



# main #

#Create Capsule
#create_capsule("Capsule", .5, 12, 3, 2, 2)

#Create Sphere
create_sphere("Sphere", 1, 48, 6)
bpy.data.objects["Sphere"].location = (0,0,1.5)
#Create Dome
create_dome("Dome", 1, 4, 3)
bpy.data.objects["Dome"].location = (-1,0,0)
create_dome("Dome", 1, 4, 3)
bpy.data.objects["Dome"].location = (1,0,0)

