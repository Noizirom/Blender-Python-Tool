
import bpy, numpy as np, os, json as js


desktop = os.path.expanduser("~/Desktop")
ob = bpy.context.object
fname = ob.name


def get_co():
    mode = bpy.context.active_object.mode
    vt = bpy.context.object.data.vertices
    bpy.ops.object.mode_set(mode='OBJECT')
    #vertices
    countv = len(vt)
    co = np.empty(countv * 3)
    vt.foreach_get('co', co)
    co.shape = (countv, 3)
    return co


def write_js(fileName, path, data):
    cwd = os.getcwd()
    os.chdir(path)
    fpath = fileName + '.json'
    with open(fpath, 'w') as f:
        js.dump(data, f)
    os.chdir(cwd)    
        
gc = get_co()
gc = np.around(gc, decimals=5)
gc = gc.tolist()


write_js(fname, desktop, gc)
        
