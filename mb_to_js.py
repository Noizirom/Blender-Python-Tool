import bpy, numpy as np, os, json as js
from copy import deepcopy as dc

desktop = os.path.expanduser("~/Desktop")
ob = bpy.context.object
fname = ob.name
vt = bpy.context.object.data.vertices
fa = bpy.context.object.data.polygons
vg = bpy.context.object.vertex_groups
countv = len(vt)
countf = len(fa)
vn = [v.name for v in bpy.context.object.vertex_groups[:]]
#wt = [[i, vt[i].groups[0].weight] for i in range(countv)]


#vertex group index list
def vg_idx_list(vgn):
    return([[v.index, v.groups[0].weight] for v in vt if vg[vgn].index in [vg.group for vg in v.groups]])

#vertex group {name: [indexes]} dictionary
def vg_idx_dict():
    vd = {n: vg_idx_list(n) for n in vn}
    return dc(vd)

def get_info():
    mode_cur = bpy.context.active_object.mode
    bpy.ops.object.mode_set(mode='OBJECT')
    #vertices
    countv = len(vt)
    co = np.empty(countv * 3)
    vt.foreach_get('co', co)
    co.shape = (countv, 3)
    fidx = np.empty(countf, dtype=np.int32)
    fa.foreach_get('index', fidx)
    fac = np.array([i.vertices[:] for i in fa])
    return co, fidx

def write_js(fileName, ext, path, data):
    cwd = os.getcwd()
    os.chdir(path)
    fpath = fileName + '_' + ext + '.json'
    with open(fpath, 'w') as f:
        js.dump(data, f)
    os.chdir(cwd)    


#main-------------------------
co, fidx = get_info()
co = np.around(co, decimals=5)
co = co.tolist()
fidx = fidx.tolist()
vid = vg_idx_dict()

write_js(fname, 'verts', desktop, co)
write_js(fname, 'polygons', desktop, fidx)
write_js(fname, 'vgroups', desktop, vid)



