from tkinter import Tk, Label, Entry, Button, Canvas, Frame, filedialog
from tkinter import *

import os
import re
import subprocess

import cv2
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from PIL import Image, ImageTk

def draw_image_uv_disps(image_path, points_ref, points_def, show_disp, scale = 1., text = None, filename = None) :
    
    assert (points_ref.shape) == (points_def.shape), 'The shape of reference points array must be the same as the shape of deformed points array!'
    
    image = cv2.imread(image_path)
	
    if text is not None :
        image = cv2.putText(image, text, (50,50), cv2.FONT_HERSHEY_SIMPLEX, 1,(255,255,255),4)
		
    if (show_disp) :
      for i, pt0_z in enumerate(points_ref[:, 2]):
          pt1_z = points_def[:, 2][i]
          if np.isnan(pt0_z)==False and np.isnan(pt1_z)==False :
              x = (points_ref[i, 0]) #int
              y = (points_ref[i, 1])
              deltat_z = pt1_z - pt0_z
              image = cv2.putText(image, str(deltat_z), (x,y), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (50,50,255), 2)

    else :
      for pt in points_ref[:, :2]:
          if not np.isnan(pt[0]) and not np.isnan(pt[1]):
              x = int(pt[0])
              y = int(pt[1])
              image = cv2.circle(image, (x, y), 4, (0, 255, 255), -1)
          
      for i, pt0 in enumerate(points_ref[:, :2]):
          pt1 = points_def[:, :2][i]
          if np.isnan(pt0[0])==False and np.isnan(pt0[1])==False and np.isnan(pt1[0])==False and np.isnan(pt1[1])==False :
              disp_x = (pt1[0]-pt0[0])*scale
              disp_y = (pt1[1]-pt0[1])*scale
              image = cv2.line(image, (int(pt0[0]), int(pt0[1])), (int(pt0[0]+disp_x), int(pt0[1]+disp_y)), (255, 120, 255) , 2)
				
    if filename is not None:
        cv2.imwrite(filename, image)
        return

    img_window = Tk()
    img_window.title("Displacements")
    orig_image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    # orig_image = cv2.pyrDown(orig_image, dstsize=(orig_image.shape[0] / 2, orig_image.shape[1] / 2) )
    orig_aspect_ratio = image.shape[1] / image.shape[0]

    resized_image = cv2.resize(orig_image, (int(1200 * orig_aspect_ratio), 1200), interpolation = cv2.INTER_LINEAR)

    canvas = Canvas(img_window, width = 1200 * orig_aspect_ratio, height = 1200, bg = "#000000")
    image_tk = ImageTk.PhotoImage(master = canvas, image=Image.fromarray(resized_image))
    canvas.pack(side="top", fill=BOTH, expand = YES)

    image_id = canvas.create_image(0, 0, image = image_tk, state = "normal", anchor = NW)

    img_window.mainloop()

def load_array_from_file(path, arr_len, arr_dim) :
  arr = np.full((arr_len, arr_dim), 0)
  with open(path, 'r') as f:

    for i in range(0, arr_len) :
      line = f.readline()
      vals =[float(j) for j in line.split()]
      arr[i, :] = np.array(vals)

  return arr


def find_numbers_in_filenames(filenames, prefix):
  """
  This function finds all unique numbers following a given prefix in filenames.

  Args:
      filenames: A list of filenames (strings).
      prefix: The prefix string to search for (string).

  Returns:
      A set of unique extracted numbers (integers).
  """
  numbers = set()  # Use a set to ensure unique elements
  pattern = rf"{prefix}(\d+)"  # Regular expression pattern with f-string
  for filename in filenames:
    match = re.search(pattern, filename)
    if match:
      number = int(match.group(1))
      numbers.add(number)  # Add number to set for unique elements
  return numbers

def get_filenames(path):
  """
  This function retrieves all filenames from a directory path.

  Args:
      path: The directory path (string).

  Returns:
      A list of filenames (strings) within the directory, or an empty list if the path is invalid.
  """
  try:
    return [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
  except FileNotFoundError:
    print(f"Error: Directory '{path}' not found.")
    return []
  

# Initialize the main window
root = Tk()
root.title("3D DIC")

def browse_directory_images_path():
  directory = filedialog.askdirectory(title="Select a Directory", mustexist=True)
  
  if directory:
    images_path_label_entry.delete(0, END)
    images_path_label_entry.insert(END,directory + '/')

def browse_directory_results_path():
  directory = filedialog.askdirectory(title="Select a Directory", mustexist=True)
  
  if directory:
    results_path_entry.delete(0, END)
    results_path_entry.insert(END,directory + '/')

def browse_directory_settings_path():
  directory = filedialog.askdirectory(title="Select a Directory", mustexist=True)
  
  if directory:
    settings_path_entry.delete(0, END)
    settings_path_entry.insert(END,directory + '/')

ignore_1st_layer_value = 0
def change_ignore_1st_layer():
  global ignore_1st_layer_value
  ignore_1st_layer_value = 1 - ignore_1st_layer_value
  ignore_1st_layer_button.config(text=str(bool(ignore_1st_layer_value)))

backward_value = 0
def change_backward():
  global backward_value
  backward_value = 1 - backward_value
  if (backward_value) :
    button_text = "Backward"
  else :
    button_text = "Forward"
  backward_value_button.config(text=button_text)

case_value = 0
def change_case():
  global case_value
  case_value = 1 - case_value
  if (case_value) :
    button_text = "2D"
  else :
    button_text = "3D"
  case_value_button.config(text=button_text)

show_disp = 0
def change_disp():
  global show_disp
  show_disp = 1 - show_disp
  if (show_disp) :
    button_text = "z"
  else :
    button_text = "uv"
  show_choose_disp_button.config(text=button_text)

def load_settings():
  filetypes = (
        ('text files', '*.txt'),
        ('All files', '*.*')
    )
  directory = filedialog.askopenfilename(title="Select a File", filetypes=filetypes)
  global ignore_1st_layer_value, backward_value, case_value
  
  # file_name = os.path.splitext(os.path.basename(directory))[0]
  path_file = os.path.split(os.path.abspath(directory))

  settings_path_entry.delete(0, END)
  settings_path_entry.insert(END, path_file[0] + '/')
  settings_name_entry.delete(0, END)
  settings_name_entry.insert(END, os.path.splitext(path_file[1])[0])

  if directory:
    with open(directory, 'r') as f:
      line = f.readline()
      subset_size_label_entry.delete(0, END)
      subset_size_label_entry.insert(END, line.split())
      line = f.readline()
      subset_offset_label_entry.delete(0, END)
      subset_offset_label_entry.insert(END, line.split())
      line = f.readline()
      z_spacing_label_entry.delete(0, END)
      z_spacing_label_entry.insert(END, line.split())
      line = f.readline()
      z_search_label_entry.delete(0, END)
      z_search_label_entry.insert(END, line.split())
      line = f.readline()
      images_path_label_entry.delete(0, END)
      images_path_label_entry.insert(END, line.split())
      line = f.readline()
      results_path_entry.delete(0, END)
      results_path_entry.insert(END, line.split())
      line = f.readline()
      images_name_prefix_label_entry.delete(0, END)
      images_name_prefix_label_entry.insert(END, line.split())
      line = f.readline()
      images_name_postfix_label_entry.delete(0, END)
      images_name_postfix_label_entry.insert(END, line.split())
      
      f.readline()
      f.readline()

      line = f.readline()
      ROI_x_min_label_entry.delete(0, END)
      ROI_x_min_label_entry.insert(END, line.split()[0])
      ROI_y_min_label_entry.delete(0, END)
      ROI_y_min_label_entry.insert(END, line.split()[1])
      line = f.readline()
      ROI_x_max_label_entry.delete(0, END)
      ROI_x_max_label_entry.insert(END, line.split()[0])
      ROI_y_max_label_entry.delete(0, END)
      ROI_y_max_label_entry.insert(END, line.split()[1])
      
      line = f.readline()
      downsampling_label_entry.delete(0, END)
      downsampling_label_entry.insert(END, line.split())

      f.readline()
      f.readline()

      line = f.readline()
      ignore_1st_layer_value = int(line.split()[0])
      ignore_1st_layer_button.config(text=str(bool(ignore_1st_layer_value)))

      line = f.readline()
      backward_value = int(line.split()[0])
      if (backward_value) :
        button_text = "Backward"
      else :
        button_text = "Forward"
      backward_value_button.config(text=button_text)

      line = f.readline()
      case_value = int(line.split()[0])
      if (case_value) :
        button_text = "2D"
      else :
        button_text = "3D"
      case_value_button.config(text=button_text)


def show_error(message) :
  error_window = Tk()
  error_window.title("Error")
  error_window.geometry("200x100")
  error_label = Label(error_window, text=message)
  error_label.grid(row=0, column=0, sticky="ew")
  error_window.mainloop()


def save_settings():
  print("Saving current settings...")

  if(images_path_label_entry.get() == "") :
    show_error("Images path is empty!")
    return
  
  if(results_path_entry.get() == "") :
    show_error("Results path is empty!")
    return
  
  if(settings_path_entry.get() == "") :
    show_error("Settins path is empty!")
    return
  
  if(settings_name_entry.get() == "") :
    show_error("Settins name is empty!")
    return
  
  if(subset_size_label_entry.get() == "") :
    show_error("Subset size is empty!")
    return
  
  if(subset_offset_label_entry.get() == "") :
    show_error("Subset offset is empty!")
    return
  
  if(z_search_label_entry.get() == "") :
    show_error("z-search radius is empty!")
    return
  
  if(z_spacing_label_entry.get() == "") :
    show_error("z-spacing is empty!")
    return
  
  if(downsampling_label_entry.get() == "") :
    show_error("downsampling is empty!")
    return
  
  if(ROI_x_min_label_entry.get() == "") :
    show_error("ROI x_min is empty!")
    return
  
  if(ROI_y_min_label_entry.get() == "") :
    show_error("ROI y_min is empty!")
    return
  
  if(ROI_x_max_label_entry.get() == "") :
    show_error("ROI x_max is empty!")
    return
  
  if(ROI_y_max_label_entry.get() == "") :
    show_error("ROI y_max is empty!")
    return
  
  if(images_name_prefix_label_entry.get() == "") :
    show_error("Images' name prefix is empty!")
    return
  
  if(images_name_postfix_label_entry.get() == "" and case_value != 0) :
    show_error("Images' name postfix is empty!")
    return
  
  image_filenames = get_filenames(images_path_label_entry.get())
  _, image_extension = os.path.splitext(image_filenames[0])
  times = find_numbers_in_filenames(image_filenames, images_name_prefix_label_entry.get())
  if (case_value) :
    stack_h = 1
  else :
    stack_hs = find_numbers_in_filenames(image_filenames, images_name_prefix_label_entry.get() + str(list(times)[0]) + images_name_postfix_label_entry.get())
    stack_h = max(list(stack_hs))

  file = settings_path_entry.get() + settings_name_entry.get() + ".txt"
  with open(file, 'w') as f:
    f.write(subset_size_label_entry.get() + "\n")
    f.write(subset_offset_label_entry.get() + "\n")
    f.write(z_spacing_label_entry.get() + "\n")
    f.write(z_search_label_entry.get() + "\n")
    f.write(images_path_label_entry.get() + "\n")
    f.write(results_path_entry.get() + "\n")
    f.write(images_name_prefix_label_entry.get() + "\n")
    if (case_value) :
      f.write("EMPTY" + "\n")
    else:
      f.write(images_name_postfix_label_entry.get() + "\n")
    f.write(str(len(times)) + "\n")
    for number in times :
      f.write(str(number) + " ")
    f.write("\n")
    f.write(ROI_x_min_label_entry.get() + " ")
    f.write(ROI_y_min_label_entry.get() + "\n")
    f.write(ROI_x_max_label_entry.get() + " ")
    f.write(ROI_y_max_label_entry.get() + "\n")
    f.write(downsampling_label_entry.get() + "\n")
    f.write(str(stack_h) + "\n")
    f.write(image_extension + "\n")
    f.write(str(ignore_1st_layer_value) + "\n")
    f.write(str(backward_value) + "\n")
    f.write(str(case_value) + "\n")


def run_dic():
  if(settings_path_entry.get() == "") :
    show_error("Settins path is empty!")
    return
  
  if(settings_name_entry.get() == "") :
    show_error("Settins name is empty!")
    return
  
  print("Running 3D DIC...")
  setting_file = settings_path_entry.get() + settings_name_entry.get() + ".txt"
  p = subprocess.Popen(["./DIC_3D", setting_file])

def show_disps():
  if(images_path_label_entry.get() == "") :
    show_error("Images path is empty!")
    return
  
  if(results_path_entry.get() == "") :
    show_error("Results path is empty!")
    return
  
  if(show_z_label_entry.get() == "") :
    show_error("Select z to show!")
    return
  
  if(show_ref_time_label_entry.get() == "") :
    show_error("Select ref time to show!")
    return
  
  if(show_def_time_label_entry.get() == "") :
    show_error("Select def time to show!")
    return
  
  if(ROI_x_min_label_entry.get() == "") :
    show_error("ROI x_min is empty!")
    return
  
  if(ROI_y_min_label_entry.get() == "") :
    show_error("ROI y_min is empty!")
    return
  
  if(ROI_x_max_label_entry.get() == "") :
    show_error("ROI x_max is empty!")
    return
  
  if(ROI_y_max_label_entry.get() == "") :
    show_error("ROI y_max is empty!")
    return
  
  if(z_spacing_label_entry.get() == "") :
    show_error("z-spacing is empty!")
    return
  
  if(subset_size_label_entry.get() == "") :
    show_error("Subset size is empty!")
    return
  
  if(subset_offset_label_entry.get() == "") :
    show_error("Subset offset is empty!")
    return
  
  if(images_name_prefix_label_entry.get() == "") :
    show_error("Images' name prefix is empty!")
    return
  
  if(images_name_postfix_label_entry.get() == "" and case_value != 0) :
    show_error("Images' name postfix is empty!")
    return
  
  global show_disp
  
  image_filenames = get_filenames(images_path_label_entry.get())
  _, image_extension = os.path.splitext(image_filenames[0])
  times = find_numbers_in_filenames(image_filenames, images_name_prefix_label_entry.get())
  if (case_value) :
    stack_h = 1
    index_to_show = 0
    z_bounce = 1
    bounced_z_id = 0
    s_z = 1
    ref_image_path = images_path_label_entry.get() + images_name_prefix_label_entry.get() + str(show_ref_time_label_entry.get()) + image_extension
    show_case = 0
  else :
    stack_hs = find_numbers_in_filenames(image_filenames, images_name_prefix_label_entry.get() + str(list(times)[0]) + images_name_postfix_label_entry.get())
    stack_h = max(list(stack_hs))
    index_to_show = int(show_z_label_entry.get())
    z_bounce = int(z_spacing_label_entry.get())
    bounced_z_id = index_to_show // z_bounce
    s_z = len(range(0, stack_h, z_bounce))
    ref_image_path = images_path_label_entry.get() + images_name_prefix_label_entry.get() + str(show_ref_time_label_entry.get()) + images_name_postfix_label_entry.get() + '{:03}'.format(int(show_z_label_entry.get())) + image_extension
    show_case = show_disp

  subset_size, subset_offset = int(subset_size_label_entry.get()), int(subset_offset_label_entry.get())
  roi_xy_min = [int(ROI_x_min_label_entry.get()), int(ROI_y_min_label_entry.get())]
  roi_xy_max = [int(ROI_x_max_label_entry.get()), int(ROI_y_max_label_entry.get())]

  s_x = int((roi_xy_max[0] - roi_xy_min[0] - subset_size) // (subset_offset)) + 1
  s_y = int((roi_xy_max[1] - roi_xy_min[1] - subset_size) // (subset_offset)) + 1

  ref_res_path = results_path_entry.get() + 'ref_' + str(show_ref_time_label_entry.get()) + '_' + str(show_def_time_label_entry.get()) + '.txt'
  def_res_path = results_path_entry.get() + 'def_' + str(show_ref_time_label_entry.get()) + '_' + str(show_def_time_label_entry.get()) + '.txt'
  # coefs_res_path = results_path_entry.get() + 'coefs_' + str(show_ref_time_label_entry.get()) + '_' + str(show_def_time_label_entry.get()) + '.txt'

  result_ref = load_array_from_file(ref_res_path, s_x*s_y*s_z, 3)
  result_def = load_array_from_file(def_res_path, s_x*s_y*s_z, 3)

  if ((index_to_show % z_bounce) != 0) :
    show_error("Result for chosen z doesn't exist!")
    return
  
  points_ref, points_def = result_ref[bounced_z_id::s_z], result_def[bounced_z_id::s_z]
  draw_image_uv_disps(ref_image_path, points_ref, points_def, show_case)



  

    

# Create frames 
text_frame = Frame(root)
text_frame.grid(row=0, column=0, sticky="ew")

settings_name_frame = Frame(root)
settings_name_frame.grid(row=1, column=0, sticky="ew")

DIC_set_button_frame = Frame(root)
DIC_set_button_frame.grid(row=2, column=0, sticky="ew")

ROI_frame = Frame(root)
ROI_frame.grid(row=3, column=0, sticky="ew")

images_name_settings_frame = Frame(root)
images_name_settings_frame.grid(row=4, column=0, sticky="ew")

buttons_frame = Frame(root)
buttons_frame.grid(row=5, column=0, sticky="ew")

run_frame = Frame(root)
run_frame.grid(row=2, column=1, sticky="ew")

show_frame = Frame(root)
show_frame.grid(row=4, column=1, sticky="ew")

# Create labels and text entry fields in the text frame
images_path_label = Label(text_frame, text="Images path:")
images_path_label.grid(row=0, column=0, sticky="ew")
images_path_label_entry = Entry(text_frame, width=100)
images_path_label_entry.grid(row=0, column=1, sticky="ew")
images_path_browse_button = Button(text_frame, text="browse", command=browse_directory_images_path)
images_path_browse_button.grid(row=0, column=2, sticky="ew")

results_path = Label(text_frame, text="Results path:")
results_path.grid(row=1, column=0, sticky="ew")
results_path_entry = Entry(text_frame, width=100)
results_path_entry.grid(row=1, column=1, sticky="ew")
results_path_browse_button = Button(text_frame, text="browse", command=browse_directory_results_path)
results_path_browse_button.grid(row=1, column=2, sticky="ew")

settings_path = Label(text_frame, text="Settings path:")
settings_path.grid(row=2, column=0, sticky="ew")
settings_path_entry = Entry(text_frame, width=100)
settings_path_entry.grid(row=2, column=1, sticky="ew")
settings_path_browse_button = Button(text_frame, text="browse", command=browse_directory_settings_path)
settings_path_browse_button.grid(row=2, column=2, sticky="ew")

settings_name = Label(settings_name_frame, text="Settings filename:")
settings_name.grid(row=3, column=0, sticky="ew")
settings_name_entry = Entry(settings_name_frame, width=30)
settings_name_entry.grid(row=3, column=1, sticky="ew")
settings_name_extension = Label(settings_name_frame, text=".txt")
settings_name_extension.grid(row=3, column=2, sticky="ew")
settings_save_button = Button(settings_name_frame, text="Save", command=save_settings)
settings_save_button.grid(row=3, column=3, sticky="ew")


subset_size_label = Label(DIC_set_button_frame, text="Subset size:")
subset_size_label.grid(row=0, column=0, sticky="ew")
subset_size_label_entry = Entry(DIC_set_button_frame, width=5)
subset_size_label_entry.grid(row=0, column=1, sticky="ew")
subset_size_end_label = Label(DIC_set_button_frame, text="px")
subset_size_end_label.grid(row=0, column=2, sticky="ew")

subset_offset_label = Label(DIC_set_button_frame, text="Subset offset:")
subset_offset_label.grid(row=1, column=0, sticky="ew")
subset_offset_label_entry = Entry(DIC_set_button_frame, width=5)
subset_offset_label_entry.grid(row=1, column=1, sticky="ew")
subset_offset_end_label = Label(DIC_set_button_frame, text="px")
subset_offset_end_label.grid(row=1, column=2, sticky="ew")

z_search_label = Label(DIC_set_button_frame, text="z-search radius:")
z_search_label.grid(row=2, column=0, sticky="ew")
z_search_label_entry = Entry(DIC_set_button_frame, width=5)
z_search_label_entry.grid(row=2, column=1, sticky="ew")
z_search_end_label = Label(DIC_set_button_frame, text="px")
z_search_end_label.grid(row=2, column=2, sticky="ew")

z_spacing_label = Label(DIC_set_button_frame, text="z-spacing:")
z_spacing_label.grid(row=3, column=0, sticky="ew")
z_spacing_label_entry = Entry(DIC_set_button_frame, width=5)
z_spacing_label_entry.grid(row=3, column=1, sticky="ew")
z_spacing_end_label = Label(DIC_set_button_frame, text="px")
z_spacing_end_label.grid(row=3, column=2, sticky="ew")

downsampling_label = Label(DIC_set_button_frame, text="downsampling:")
downsampling_label.grid(row=4, column=0, sticky="ew")
downsampling_label_entry = Entry(DIC_set_button_frame, width=5)
downsampling_label_entry.grid(row=4, column=1, sticky="ew")
downsampling_end_label = Label(DIC_set_button_frame, text="px")
downsampling_end_label.grid(row=4, column=2, sticky="ew")

ignore_1st_layer_label = Label(buttons_frame, text="ignore 1st layer:")
ignore_1st_layer_label.grid(row=0, column=0, sticky="ew")
ignore_1st_layer_button = Button(buttons_frame, text=str(bool(ignore_1st_layer_value)), command=change_ignore_1st_layer)
ignore_1st_layer_button.grid(row=0, column=1, sticky="ew")

backward_value_label = Label(buttons_frame, text="Calculation:")
backward_value_label.grid(row=1, column=0, sticky="ew")
backward_value_button = Button(buttons_frame, text="Forward", command=change_backward)
backward_value_button.grid(row=1, column=1, sticky="ew")

case_value_label = Label(buttons_frame, text="Case:")
case_value_label.grid(row=2, column=0, sticky="ew")
case_value_button = Button(buttons_frame, text="3D", command=change_case)
case_value_button.grid(row=2, column=1, sticky="ew")

# SHOW
show_text_label = Label(show_frame, text="Show displacements")
show_text_label.grid(row=0, column=0, sticky="ew")
show_ref_time_label = Label(show_frame, text="ref time:")
show_ref_time_label.grid(row=1, column=0, sticky="ew")
show_ref_time_label_entry = Entry(show_frame, width=5)
show_ref_time_label_entry.grid(row=1, column=1, sticky="ew")

show_def_time_label = Label(show_frame, text=" def time:")
show_def_time_label.grid(row=1, column=2, sticky="ew")
show_def_time_label_entry = Entry(show_frame, width=5)
show_def_time_label_entry.grid(row=1, column=3, sticky="ew")

show_z_label = Label(show_frame, text="z:")
show_z_label.grid(row=2, column=0, sticky="ew")
show_z_label_entry = Entry(show_frame, width=5)
show_z_label_entry.grid(row=2, column=1, sticky="ew")

show_choose_disp_label = Label(show_frame, text="Disp:")
show_choose_disp_label.grid(row=3, column=0, sticky="ew")
show_choose_disp_button = Button(show_frame, text="uv", command=change_disp)
show_choose_disp_button.grid(row=3, column=1, sticky="ew")

show_button = Button(show_frame, text="Show", command=show_disps)
show_button.grid(row=4, column=0, sticky="ew")

# ROI
ROI_label = Label(ROI_frame, text="ROI   ")
ROI_label.grid(row=0, column=0, sticky="ew")
ROI_x_min_label = Label(ROI_frame, text="x_min:")
ROI_x_min_label.grid(row=0, column=1, sticky="ew")
ROI_x_min_label_entry = Entry(ROI_frame, width=5)
ROI_x_min_label_entry.grid(row=0, column=2, sticky="ew")
ROI_y_min_label = Label(ROI_frame, text="y_min:")
ROI_y_min_label.grid(row=0, column=3, sticky="ew")
ROI_y_min_label_entry = Entry(ROI_frame, width=5)
ROI_y_min_label_entry.grid(row=0, column=4, sticky="ew")
ROI_x_max_label = Label(ROI_frame, text="x_max:")
ROI_x_max_label.grid(row=0, column=5, sticky="ew")
ROI_x_max_label_entry = Entry(ROI_frame, width=5)
ROI_x_max_label_entry.grid(row=0, column=6, sticky="ew")
ROI_y_max_label = Label(ROI_frame, text="y_max:")
ROI_y_max_label.grid(row=0, column=7, sticky="ew")
ROI_y_max_label_entry = Entry(ROI_frame, width=5)
ROI_y_max_label_entry.grid(row=0, column=8, sticky="ew")

images_name_prefix_label = Label(images_name_settings_frame, text="Image time prefix:")
images_name_prefix_label.grid(row=0, column=0, sticky="ew")
images_name_prefix_label_entry = Entry(images_name_settings_frame, width=7)
images_name_prefix_label_entry.grid(row=0, column=1, sticky="ew")
images_name_postfix_label = Label(images_name_settings_frame, text="Image slice postfix:")
images_name_postfix_label.grid(row=0, column=2, sticky="ew")
images_name_postfix_label_entry = Entry(images_name_settings_frame, width=7)
images_name_postfix_label_entry.grid(row=0, column=3, sticky="ew")


# RUN DIC
load_settings_button = Button(run_frame, text="Load Settings", command=load_settings)  # Define browse_pt_image function later
load_settings_button.grid(row=0, column=0, sticky="ew")
run_label_button = Button(run_frame, text="Run DIC", command=run_dic)  # Define browse_pt_image function later
run_label_button.grid(row=1, column=0, sticky="ew")

# Run the main loop
root.mainloop()
