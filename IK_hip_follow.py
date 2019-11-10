import bpy

arm = bpy.context.object
pb = arm.pose.bones
obj = bpy.data.objects


emp = obj.new("Hip_Follow", None)
bpy.context.collection.objects.link(emp)
emp.empty_display_size = 0.01

bc = emp.constraints.new(type='COPY_LOCATION')
bc.use_z = False
bc.target = arm
bc.subtarget = "IK_control_hip_pos"

ac = arm.constraints.new(type='COPY_LOCATION')
ac.target = emp

