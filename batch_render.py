import bpy 
import os 

# folder path to save rendered image
in_dir = "D:\\program\\sim_render\\image\\"
lst = os.listdir(in_dir)

# folder path for importing data files
in_dir_ply ="D:\\program\\sim_render\\data\\"
lst_ply = os.listdir(in_dir_ply)

rendered_img = []
for item in lst:
    fileName, fileExtension = os.path.splitext(item)
    if 'png' in fileExtension:
        rendered_img.append(fileName)

# folder path to save rendered animation
out_dir = "D:\\program\\sim_render\\animation\\test.avi"

# Filter file list by valid file types.
candidates = []
candidates_name = []
c = 0
for item in lst_ply:
    fileName, fileExtension = os.path.splitext(lst_ply[c])
    if fileExtension == ".ply":
        candidates.append(item)
        candidates_name.append(fileName)
    c = c + 1

file = [{"name":i} for i in candidates]   
n = len(file)
print(n) 

n = 60
# To import mesh.ply in batches
bpy.data.objects['Cube'].hide = True
bpy.data.objects['Cube'].hide_render = True

# for i in range (0,n):
#     if candidates_name[i] in rendered_img:
#         continue
#     bpy.ops.import_mesh.ply(filepath=candidates[i], files=[file[i]], directory=in_dir_ply, filter_glob="*.ply")
#     bpy.data.objects[candidates_name[i]].hide = True
#     bpy.data.objects[candidates_name[i]].hide_render = True

#     bpy.data.objects[candidates_name[i]].rotation_euler = (3.14/2, 0, 0)
#     bpy.data.objects[candidates_name[i]].scale = (0.25, 0.25, 0.25)
#     bpy.data.objects[candidates_name[i]].location = (0, 0, 0.75)

# obj.rotation_euler = (3.14/2, 0, 0)
# obj.scale = (0.25, 0.25, 0.25)
# Set file_format for render images
bpy.data.scenes["Scene"].render.image_settings.file_format = 'PNG'

# To render and save rendered images
for i in range (0,n):
    if candidates_name[i] in rendered_img:
        continue
    bpy.ops.import_mesh.ply(filepath=candidates[i], files=[file[i]], directory=in_dir_ply, filter_glob="*.ply")
    
    bpy.data.objects[candidates_name[i]].rotation_euler = (3.14/2, 0, 0)
    bpy.data.objects[candidates_name[i]].scale = (0.25, 0.25, 0.25)
    bpy.data.objects[candidates_name[i]].location = (0, 0, 0.75)
    #objects must be visible to use modifier
    bpy.data.objects[candidates_name[i]].hide = False    
    #objects must be renderable to export render image
    bpy.data.objects[candidates_name[i]].hide_render = False    
    #get object
    bpy.context.scene.objects.active = bpy.data.objects[candidates_name[i]] 
    #add modifier as particle_system
    #bpy.ops.object.modifier_add(type='PARTICLE_SYSTEM')    
    #assign particle settings to object's particle system
    #bpy.data.objects[candidates_name[i]].particle_systems['ParticleSystem'].settings = bpy.data.particles['ParticleSettings'] 
    #add material
    bpy.data.objects[candidates_name[i]].data.materials.append(bpy.data.materials["Material.001"])
    #set save filepath
    bpy.data.scenes["Scene"].render.filepath = in_dir + candidates_name[i]   
    #render and save
    bpy.ops.render.render( write_still=True )    
    #hide again for next image rendering
    bpy.data.objects[candidates_name[i]].hide = True    
    #hide again for next image rendering
    bpy.data.objects[candidates_name[i]].hide_render = True
    bpy.data.objects.remove(bpy.data.objects[candidates_name[i]])    
