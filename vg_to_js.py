import bpy, numpy as np, os, json as js
from copy import deepcopy as dc

desktop = os.path.expanduser("~/Desktop")
ob = bpy.context.object
fname = ob.name

#vertex group index list
def vg_idx_list(vgn):
    return([[v.index, v.groups[0].weight] for v in bpy.context.object.data.vertices if bpy.context.object.vertex_groups[vgn].index in [vg.group for vg in v.groups]])

#vertex group {name: [indexes]} dictionary
def vg_idx_dict():
    vn = [v.name for v in bpy.context.object.vertex_groups[:]]
    vd = {n: vg_idx_list(n) for n in vn}
    return dc(vd)

def write_js(fileName, path, data):
    cwd = os.getcwd()
    os.chdir(path)
    fpath = fileName + '_vertex_groups.json'
    with open(fpath, 'w') as f:
        js.dump(data, f)
    os.chdir(cwd)    
        
vid = vg_idx_dict()


write_js(fname, desktop, vid)










