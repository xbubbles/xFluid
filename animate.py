import bpy
import os

# folder path to save rendered image
in_dir = "D:\\program\\sim_render\\image\\demo2\\"
lst = os.listdir(in_dir)
lst.sort(key = lambda x:int(x[:-4]))
# folder path to save rendered animation
out_dir = "D:\\program\\sim_render\\animation\\test3.avi"

# Active VSE to generate rendering animation
bpy.data.scenes["Scene"].render.use_sequencer = True

# Filter file list by valid file types.
re_image = []
re_image_name = []
c = 0
for item in lst:
    fileName, fileExtension = os.path.splitext(lst[c])
    if fileExtension == ".png":
        re_image.append(item)
        re_image_name.append(fileName)
    c = c + 1

# Create the sequencer data
bpy.context.scene.sequence_editor_create()

n = len(lst)
# Add strip into VSE by importing new image
for i in range (0, n):
    bpy.context.scene.sequence_editor.sequences.new_image(
        name=re_image[i],
        filepath=os.path.join(in_dir, re_image[i]),
        channel=1, frame_start=i)

# Resolution settings for animation
resx = 720; #1920 
resy = 480; #1080 
bpy.data.scenes["Scene"].render.resolution_x = resx 
bpy.data.scenes["Scene"].render.resolution_y = resy 
bpy.data.scenes["Scene"].render.resolution_percentage = 100 

bpy.data.scenes["Scene"].frame_end = n 
bpy.data.scenes["Scene"].render.image_settings.file_format = 'AVI_JPEG' 
bpy.data.scenes["Scene"].render.filepath = out_dir
bpy.ops.render.render( animation=True )