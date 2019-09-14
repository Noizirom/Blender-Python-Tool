import bpy, numpy as np, os, json as js
from copy import deepcopy as dc
from mathutils import Matrix, Vector, Euler
from math import radians, degrees

cap_co = np.array([[0.38268, 0.0, 1.92388], [0.2706, 0.2706, 1.92388], [0.0, 0.38268, 1.92388], [-0.2706, 0.2706, 1.92388], [-0.38268, 0.0, 1.92388], [-0.2706, -0.2706, 1.92388], [-0.0, -0.38268, 1.92388], [0.2706, -0.2706, 1.92388], [0.70711, 0.0, 1.70711], [0.5, 0.5, 1.70711], [0.0, 0.70711, 1.70711], [-0.5, 0.5, 1.70711], [-0.70711, 0.0, 1.70711], [-0.5, -0.5, 1.70711], [-0.0, -0.70711, 1.70711], [0.5, -0.5, 1.70711], [0.92388, 0.0, 1.38268], [0.65328, 0.65328, 1.38268], [0.0, 0.92388, 1.38268], [-0.65328, 0.65328, 1.38268], [-0.92388, 0.0, 1.38268], [-0.65328, -0.65328, 1.38268], [-0.0, -0.92388, 1.38268], [0.65328, -0.65328, 1.38268], [1.0, 0.0, 1.0], [0.70711, 0.70711, 1.0], [0.0, 1.0, 1.0], [-0.70711, 0.70711, 1.0], [-1.0, 0.0, 1.0], [-0.70711, -0.70711, 1.0], [-0.0, -1.0, 1.0], [0.70711, -0.70711, 1.0], [0.0, 0.0, 2.0], [0.38268, 0.0, -1.92388], [0.2706, 0.2706, -1.92388], [0.0, 0.38268, -1.92388], [-0.2706, 0.2706, -1.92388], [-0.38268, 0.0, -1.92388], [-0.2706, -0.2706, -1.92388], [-0.0, -0.38268, -1.92388], [0.2706, -0.2706, -1.92388], [0.70711, 0.0, -1.70711], [0.5, 0.5, -1.70711], [0.0, 0.70711, -1.70711], [-0.5, 0.5, -1.70711], [-0.70711, 0.0, -1.70711], [-0.5, -0.5, -1.70711], [-0.0, -0.70711, -1.70711], [0.5, -0.5, -1.70711], [0.92388, 0.0, -1.38268], [0.65328, 0.65328, -1.38268], [0.0, 0.92388, -1.38268], [-0.65328, 0.65328, -1.38268], [-0.92388, 0.0, -1.38268], [-0.65328, -0.65328, -1.38268], [-0.0, -0.92388, -1.38268], [0.65328, -0.65328, -1.38268], [1.0, 0.0, -1.0], [0.70711, 0.70711, -1.0], [0.0, 1.0, -1.0], [-0.70711, 0.70711, -1.0], [-1.0, 0.0, -1.0], [-0.70711, -0.70711, -1.0], [-0.0, -1.0, -1.0], [0.70711, -0.70711, -1.0], [0.0, 0.0, -2.0]])
cap_fa = [(24, 25, 17, 16), (58, 57, 49, 50), (25, 26, 18, 17), (59, 58, 50, 51), (26, 27, 19, 18), (60, 59, 51, 52), (27, 28, 20, 19), (61, 60, 52, 53), (28, 29, 21, 20), (62, 61, 53, 54), (29, 30, 22, 21), (63, 62, 54, 55), (30, 31, 23, 22), (64, 63, 55, 56), (31, 24, 16, 23), (57, 64, 56, 49), (16, 17, 9, 8), (50, 49, 41, 42), (17, 18, 10, 9), (51, 50, 42, 43), (18, 19, 11, 10), (52, 51, 43, 44), (19, 20, 12, 11), (53, 52, 44, 45), (20, 21, 13, 12), (54, 53, 45, 46), (21, 22, 14, 13), (55, 54, 46, 47), (22, 23, 15, 14), (56, 55, 47, 48), (23, 16, 8, 15), (49, 56, 48, 41), (8, 9, 1, 0), (42, 41, 33, 34), (9, 10, 2, 1), (43, 42, 34, 35), (10, 11, 3, 2), (44, 43, 35, 36), (11, 12, 4, 3), (45, 44, 36, 37), (12, 13, 5, 4), (46, 45, 37, 38), (13, 14, 6, 5), (47, 46, 38, 39), (14, 15, 7, 6), (48, 47, 39, 40), (15, 8, 0, 7), (41, 48, 40, 33), (25, 24, 57, 58), (26, 25, 58, 59), (27, 26, 59, 60), (28, 27, 60, 61), (29, 28, 61, 62), (30, 29, 62, 63), (31, 30, 63, 64), (24, 31, 64, 57), (0, 1, 32), (1, 2, 32), (2, 3, 32), (3, 4, 32), (4, 5, 32), (5, 6, 32), (6, 7, 32), (7, 0, 32), (34, 33, 65), (35, 34, 65), (36, 35, 65), (37, 36, 65), (38, 37, 65), (39, 38, 65), (40, 39, 65), (33, 40, 65)]

arm = bpy.context.object
bones = bpy.data.armatures[0].bones
pb = arm.pose.bones
amw = arm.matrix_world.copy()
infl = 1

ext = 'UCX_'
tgt = '_target'

xtrb = ['twist', 'root', 'muscle', 'breast']
hand =['hand', 'index', 'middle', 'ring', 'thumb', 'pinky']
hdb = ['head']
pel = ['pelvis']
torso = ['neck', 'spine', 'clavicle']
legs = ['thigh', 'calf']
foot = ['foot', 'toe']



#

rd_dict = {
            'head': [-37, 22, -45, 45, -30, 30, 'torso'],
            'neck': [-37, 22, -45, 45, -30, 30, ''],
            'clavicle_L': [0, 0, -30, 10, -30, 30, ''],
            'clavicle_R': [0, 0, -10, 30, -30, 30, ''],
            'torso': [-45, 68, -45, 45, -30, 30, 'pelvis'],
            'pelvis': [0, 0, 0, 0, 0, 0, ''],
            'spine01': [0, 0, 0, 0, 0, 0, ''],
            'spine02': [-45, 68, -45, 45, -30, 30,''],
            'spine03': [-45, 22, -45, 45, -30, 30, ''],
            'upperarm_L': [-135, 90, -98, 105, -97, 91, 'torso'],
            'upperarm_R': [-135, 90, -105, 98, -91, 97, 'torso'],
            'lowerarm_L': [-90, 79, 0, 0, -146, 0, 'upperarm_L'],
            'lowerarm_R': [-90, 79, 0, 0, 0, 146, 'upperarm_R'],
            #'hand_L': [-45, 45, -90, 86, -25, 36, 'lowerarm_L'],
            #'hand_R': [-45, 45, -86, 90, -36, 25, 'lowerarm_R'],
            'thigh_L': [-155, 45, -105, 85, -17, 88, 'pelvis'],
            'thigh_R': [-155, 45, -85, 105, -88, 17, 'pelvis'],
            'calf_L': [0, 150, 0, 0, 0, 0, 'thigh_L'],
            'calf_R': [0, 150, 0, 0, 0, 0, 'thigh_R'],
            'foot_L': [-31, 63, -26, 26, -15, 74, 'calf_L'],
            'foot_R': [-31, 63, -26, 26, -74, 15, 'calf_R']
            }

###############################################################################################################################

###############################################################################################################################
#ADD RIGIDBODY CONSTRAINTS
try:
    bpy.ops.rigidbody.world_remove()
except:
    pass
bpy.ops.rigidbody.world_add()


def add_rbc(part):
    if not rd_dict[part][6]:
        pass
    else:
        #rg_list
        rig = part + tgt
        p1 = ext + rd_dict[part][6]
        p2 = ext + part
        #rd_dict[part]
        xl = rd_dict[part][0]
        xu = rd_dict[part][1]
        yl = rd_dict[part][2]
        yu = rd_dict[part][3]
        zl = rd_dict[part][4]
        zu = rd_dict[part][5]
        bpy.ops.object.select_all(action='DESELECT')
        bpy.data.objects[rig].select_set(state=True)
        bpy.context.view_layer.objects.active = bpy.data.objects[rig]
        bpy.ops.rigidbody.constraint_add(type = 'GENERIC')
        #
        rbc = bpy.context.object.rigid_body_constraint
        rbc.object1 = bpy.data.objects[p1]
        rbc.object2 = bpy.data.objects[p2]
        ###use limits###
        #angle
        rbc.use_limit_ang_x = True
        rbc.use_limit_ang_y = True
        rbc.use_limit_ang_z = True
        rbc.use_limit_lin_x = True
        rbc.use_limit_lin_y = True
        rbc.use_limit_lin_z = True
        #Linear
        rbc.limit_lin_x_lower = 0
        rbc.limit_lin_x_upper = 0
        rbc.limit_lin_y_lower = 0
        rbc.limit_lin_y_upper = 0
        rbc.limit_lin_z_lower = 0
        rbc.limit_lin_z_upper = 0
        #Rotation Constraints
        rbc.limit_ang_x_lower = radians(xl)
        rbc.limit_ang_x_upper = radians(xu)
        rbc.limit_ang_y_lower = radians(yl)
        rbc.limit_ang_y_upper = radians(yu)
        rbc.limit_ang_z_lower = radians(zl)
        rbc.limit_ang_z_upper = radians(zu)
        bpy.data.objects[rig].select_set(state=False)
        bpy.data.objects[p2].select_set(state=True)
        print(rig+"<<<finished>>>")






def rd_capsule(length, radius):
    sc = np.eye(3)
    sc[0][0] = radius
    sc[1][1] = radius
    sc[2][2] = length/4
    coo = np.array([[i[0], i[1], i[2] + 2] for i in cap_co])
    nc = coo @ sc
    return nc

def obj_mesh(co, faces):
    mesh = bpy.data.meshes.new("Obj")
    mesh.from_pydata(co, [], faces)
    mesh.validate()
    mesh.update(calc_edges = True)
    Object = bpy.data.objects.new("Obj", mesh)
    Object.data = mesh
    bpy.context.collection.objects.link(Object)
    bpy.context.view_layer.objects.active = Object
    Object.select_set(True)

def obj_new(Name, co, faces):
    obj_mesh(co, faces)
    bpy.data.objects["Obj"].name = Name
    bpy.data.meshes[bpy.data.objects[Name].data.name].name = Name

def add_capsule(Name, length, radius):
    cor = rd_capsule(length, radius)
    obj_new(Name, cor, cap_fa)
    try:
        bpy.ops.rigidbody.objects_add(type='ACTIVE')
        #bpy.context.view_layer.objects.active.rigid_body.type=True
        bpy.ops.rigidbody.mass_calculate(material='Beans (Soy)', density=721)
        bpy.context.view_layer.objects.active = bpy.data.objects[p2]
        bpy.data.objects[p2].select_set(state=False)
    except:
        pass
        
def bone_dict():
    '''
    {name: [ #0 index, #1 head, #2 head_local, #3 head_radius, 
    #4 center, #5 tail, #6 tail_local, #7 tail_radius, #8 length, 
    #9 matrix, #10 matrix_local, #11 x_axis, #12 y_axis, #13 z_axis,
    #14 vector, #15 parent ]}
    '''
    bd = {bone.name: 
    [idx, bone.head, bone.head_local, bone.head_radius, 
    bone.center, bone.tail, bone.tail_local, bone.tail_radius, bone.length, 
    bone.matrix, bone.matrix_local ,bone.x_axis, bone.y_axis, bone.z_axis,
    bone.vector, bone.parent] 
    for idx, bone in enumerate(bpy.data.armatures[0].bones)}
    return bd

def flat_list(List):
    return [i[j] for i in List for j in range(len(i))]


def list_inclusion(List, Ref):
    List, Ref = list(List), list(Ref)
    return [[i for i in List if j in i] for j in Ref]

def list_exclusion(List, Ref):
    List, Ref = list(List), list(Ref)
    return [i for i in List if not any(j in i for j in Ref)]

def name_list(Ob):
    return [i.name for i in Ob]

def name_rm_ext(Name, ext):
    return Name.replace(ext, '')

def rotation_matrix(xrot, yrot, zrot):
    rot_mat = np.array(
                        [ [np.cos(xrot)*np.cos(yrot), -np.sin(xrot)*np.cos(zrot) + np.cos(xrot)*np.sin(yrot)*np.sin(zrot), np.sin(xrot)*np.sin(zrot) + np.cos(xrot)*np.sin(yrot)*np.cos(zrot)],
                        [ np.sin(xrot)*np.cos(yrot), np.cos(xrot)*np.cos(zrot) + np.sin(xrot)*np.sin(yrot)*np.sin(zrot), -np.cos(xrot)*np.sin(zrot) + np.sin(xrot)*np.sin(yrot)*np.cos(zrot)],
                        [-np.sin(yrot), np.cos(yrot)*np.sin(zrot), np.cos(yrot)*np.cos(zrot)] ]
                        )
    return rot_mat

def apply_mod(Ref):
    act = bpy.context.view_layer.objects.active
    for o in bpy.context.view_layer.objects:
        for m in o.modifiers:
            if Ref in m.name:
                bpy.context.view_layer.objects.active = o
                bpy.ops.object.modifier_apply(modifier=m.name)
    bpy.context.view_layer.objects.active = act

def obj_del(List):
    bpy.ops.object.select_all(action='DESELECT')
    for o in List:
        obj[o].select_set(state=True)
        bpy.ops.object.delete()

def vec_dist(v1, v2):
    v1 = np.array(v1)
    v2 = np.array(v2)
    return np.sqrt(np.sum((v1 - v2)**2))

def torso_obj(ext):
    n = ext + 'torso'
    l = vec_dist(bd['spine01'][1], bd['neck'][5])
    r = bd['clavicle_L'][8]
    add_capsule(n, l*2.5, r*1.25)
    bpy.data.objects[n].matrix_world = pb['spine01'].matrix @ Matrix(rotation_matrix(0,0,radians(-90))).to_4x4()
    bpy.data.objects[n].display_type = 'WIRE'

def head_obj(ext):
    n = ext + 'head'
    add_capsule(n, bd['head'][8]*1.25, bd['head'][8]/2)
    bpy.data.objects[n].matrix_world = pb['head'].matrix @ Matrix(rotation_matrix(0,0,radians(-90))).to_4x4()
    bpy.data.objects[n].display_type = 'WIRE'

def pelvis_obj(ext):
    n = ext + 'pelvis'
    add_capsule(n, bd['pelvis'][8], bd['pelvis'][8]*.7)
    bpy.data.objects[n].matrix_world = pb['pelvis'].matrix @ Matrix(rotation_matrix(0,0,radians(-90))).to_4x4()
    bpy.data.objects[n].display_type = 'WIRE'

def ragdoll_obj(List, ext):
    for b in List:
        rdn = ext + b
        add_capsule(rdn, bd[b][8], bd[b][7])
        bpy.data.objects[rdn].matrix_world = pb[b].matrix @ Matrix(rotation_matrix(0,0,radians(-90))).to_4x4()
        bpy.data.objects[rdn].display_type = 'WIRE'

def add_emp(pbn):
    tmp_final = pb[pbn].matrix.copy() @ Matrix(rotation_matrix(0,0,radians(-90))).to_4x4()
    matrix_final = tmp_final @ Matrix(rotation_matrix(0,0,radians(90))).to_4x4()
    obj_empty = bpy.data.objects.new(pbn + tgt, None)
    bpy.context.collection.objects.link(obj_empty)
    #draw size
    obj_empty.empty_display_size = 0.1
    obj_empty.matrix_world = matrix_final

def add_parent(rd):
    objects = bpy.data.objects
    bpy.ops.object.select_all(action='DESELECT')
    a = objects[ext + rd]
    b = objects[rd + tgt]
    a.select_set(state=True)
    b.select_set(state=True)
    bpy.context.view_layer.objects.active = a
    bpy.ops.object.parent_set(type='OBJECT', keep_transform=True)



   
#main
bd = bone_dict()
rl = xtrb + hand + foot
b_list = name_list(pb)
bl = list_exclusion(b_list, rl)
b_xtrb = list_inclusion(b_list, xtrb)
b_head = list_inclusion(b_list, ['head'])
b_torso = list_inclusion(b_list, torso)
b_arms = [i for i in bl if 'arm' in i]
b_hand = list_inclusion(b_list, hand)
b_leg = list_inclusion(bl, legs)
b_foot = list_inclusion(b_list, foot)
coll = bpy.data.collections
col_nme = "Ragdoll"
#col_items = [i.name for i in coll[:]]
rdl = b_arms + b_leg[0] + b_leg[1] + b_foot[0]
emp_list = rdl 


cen_mass = ["pelvis", "spine01", "spine02", "spine03", "clavicle_L", "clavicle_R", "neck", "head"]



obj = bpy.data.objects
ragdoll_obj(rdl, ext)
pelvis_obj(ext)
head_obj(ext)
torso_obj(ext)


for b in emp_list:
    add_emp(b)

add_emp("pelvis" )
add_emp("spine01" )
obj["spine01" + tgt].name = "torso" + tgt
add_emp("spine01" )
add_emp("spine02" )
add_emp("spine03" )
add_emp("neck" )
add_emp("clavicle_L" )
add_emp("clavicle_R" )
add_emp("head")

for bn in emp_list:
    add_parent(bn)

add_parent("pelvis")
add_parent("torso")
add_parent("head")



bpy.ops.object.select_all(action='DESELECT')
par = obj[ext + "torso"]
c1 = obj["spine01" + tgt]
c2 = obj["spine02" + tgt]
c3 = obj["spine03" + tgt]
c4 = obj["neck" + tgt]
c5 = obj["clavicle_L" + tgt]
c6 = obj["clavicle_R" + tgt]
par.select_set(state=True)
c1.select_set(state=True)
c2.select_set(state=True)
c3.select_set(state=True)
c4.select_set(state=True)
c5.select_set(state=True)
c6.select_set(state=True)
bpy.context.view_layer.objects.active = par
bpy.ops.object.parent_set(type='OBJECT', keep_transform=True)




for o in emp_list:
    obj[ext + o].select_set(state=True)
    obj[o + tgt].select_set(state=True)

obj[ext + "torso"].select_set(state=True)
obj["torso" + tgt].select_set(state=True)
obj[ext + "head"].select_set(state=True)
obj["head" + tgt].select_set(state=True)
obj[ext + "pelvis"].select_set(state=True)
obj["pelvis" + tgt].select_set(state=True)


arm.select_set(state=False)

try:
    if bpy.data.collections[col_nme]:
        bpy.data.collections.remove(bpy.data.collections[col_nme])
except:
    pass
bpy.ops.object.move_to_collection(collection_index=0, is_new=True, new_collection_name=col_nme)


#Add 'COPY_TRANSFORMS' constraint to all target bones
#bpy.ops.object.mode_set(mode='POSE')
for bo in pb:
    if bo.name in emp_list:
        nc3 = bo.constraints.new(type='COPY_TRANSFORMS')
        nc3.target = bpy.data.objects[bo.name + tgt]
        nc3.influence = infl
    if bo.name in cen_mass:
        nc3 = bo.constraints.new(type='COPY_TRANSFORMS')
        nc3.target = bpy.data.objects[bo.name + tgt]
        nc3.influence = infl




for p in rd_dict:
    if not rd_dict[p][6]:
        pass
    else:
        add_rbc(p)

bpy.ops.object.select_all(action='DESELECT')

