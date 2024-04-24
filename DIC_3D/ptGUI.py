from tkinter import Tk, Label, Entry, Button, Canvas, Frame, filedialog
from tkinter import *

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rc('figure',  figsize=(10, 6))
mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"

mpl.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import os
import re
import numpy as np
import pandas as pd
import pickle
import pims
import trackpy as tp

from PIL import Image, ImageTk

import warnings
warnings.filterwarnings("ignore", module="matplotlib")

from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
############################# helpers ##############################
class LinearNDInterpolatorExt(object):
    
  def __init__(self, points,values):
    self.funcinterp  = LinearNDInterpolator(points,values)
    self.funcnearest = NearestNDInterpolator(points,values)
    
  def __call__(self,*args):
    t = self.funcinterp(*args)
    if not np.isnan(t).any():
      return t
    else:
      return self.funcnearest(*args)
    
# def load_array_from_file(path, arr_len, arr_dim) :
#   arr = np.full((arr_len, arr_dim), 0)
#   with open(path, 'r') as f:

#     for i in range(0, arr_len) :
#       line = f.readline()
#       vals =[float(j) for j in line.split()]
#       arr[i, :] = np.array(vals)

#   return arr

def closest_node(node, nodes, m): # Nearest m particels
    nodes = np.array(nodes)
    dist_2 = np.sqrt(np.sum((nodes - node)**2, axis=1))
    return np.argsort(dist_2)[1:m+1]

def RemoveWrongMatch(grid_1, grid_0,corr_max, k, q1, q3, corr_th, max_ux):
    grid_0 = grid_0[~np.isnan(corr_max),:]
    grid_1 = grid_1[~np.isnan(corr_max),:]
    corr_max = corr_max[~np.isnan(corr_max)]
    disp = grid_1 - grid_0
    
    mask = np.full(len(grid_0),True)
    mask[corr_max > corr_th] = False
        
    mask[np.abs(disp[:,0])>max_ux] = False
    # for node_id, node in enumerate(grid_0):
    #     neighbors = closest_node(node, grid_0, k) # compare the nearest n neighbors
    #     Q1 = np.quantile(disp[neighbors,0], q1)
    #     Q3 = np.quantile(disp[neighbors,0], q3)
    #     IQR = Q3-Q1
    #     if disp[node_id,0]>Q3+1.5*IQR or disp[node_id,0]<Q1-1.5*IQR:
    #         mask[node_id]=False
    
    grid_0 = grid_0[mask]
    grid_1 = grid_1[mask]
    corr_max = corr_max[mask]
    return grid_0, grid_1, corr_max

####################################################################


def browse_directory_images_path():
  directory = filedialog.askdirectory(title="Select a Directory", mustexist=True)
  
  if directory:
    images_path_label_entry.delete(0, END)
    images_path_label_entry.insert(END,directory + '/')

def browse_directory_dic_results_path():
  directory = filedialog.askdirectory(title="Select a Directory", mustexist=True)
  
  if directory:
    dic_results_path_entry.delete(0, END)
    dic_results_path_entry.insert(END,directory + '/')

def browse_directory_locations_path():
  directory = filedialog.askdirectory(title="Select a Directory", mustexist=True)
  
  if directory:
    locations_path_entry.delete(0, END)
    locations_path_entry.insert(END,directory + '/')

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
  
  if(dic_results_path_entry.get() == "") :
    show_error("DIC results path is empty!")
    return
  
  if(locations_path_entry.get() == "") :
    show_error("Particles' locations path is empty!")
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
  
  if(scale_x_size_label_entry.get() == "") :
    show_error("Scale x is empty!")
    return
  
  if(scale_y_size_label_entry.get() == "") :
    show_error("Scale y is empty!")
    return
  
  if(scale_z_size_label_entry.get() == "") :
    show_error("Scale z is empty!")
    return
  
  if(scale_z_size_label_entry.get() == "") :
    show_error("Scale z is empty!")
    return
  
  if(pt_search_range_label_entry.get() == "") :
    show_error("PT seacrh range is empty!")
    return
  
  if(images_name_prefix_label_entry.get() == "") :
    show_error("Images' name prefix is empty!")
    return
  
  # if(images_name_postfix_label_entry.get() == "") :
  #   show_error("Images' name postfix is empty!")
  #   return

  file = settings_path_entry.get() + settings_name_entry.get() + ".txt"
  with open(file, 'w') as f:
    f.write(images_path_label_entry.get() + "\n")
    f.write(dic_results_path_entry.get() + "\n")
    f.write(locations_path_entry.get() + "\n")
    f.write(results_path_entry.get() + "\n")
    f.write(images_name_prefix_label_entry.get() + "\n")
    f.write(images_name_postfix_label_entry.get() + "\n")
    f.write(scale_x_size_label_entry.get() + "\n")
    f.write(scale_y_size_label_entry.get() + "\n")
    f.write(scale_z_size_label_entry.get() + "\n")
    f.write(pt_search_range_label_entry.get() + "\n")

def load_settings():
  filetypes = (
        ('text files', '*.txt'),
        ('All files', '*.*')
    )
  directory = filedialog.askopenfilename(title="Select a File", filetypes=filetypes)
  
  # file_name = os.path.splitext(os.path.basename(directory))[0]
  path_file = os.path.split(os.path.abspath(directory))

  settings_path_entry.delete(0, END)
  settings_path_entry.insert(END, path_file[0] + '/')
  settings_name_entry.delete(0, END)
  settings_name_entry.insert(END, os.path.splitext(path_file[1])[0])

  if directory:
    with open(directory, 'r') as f:
      
      line = f.readline()
      images_path_label_entry.delete(0, END)
      images_path_label_entry.insert(END, line.split())

      line = f.readline()
      dic_results_path_entry.delete(0, END)
      dic_results_path_entry.insert(END, line.split())

      line = f.readline()
      locations_path_entry.delete(0, END)
      locations_path_entry.insert(END, line.split())

      line = f.readline()
      results_path_entry.delete(0, END)
      results_path_entry.insert(END, line.split())

      line = f.readline()
      images_name_prefix_label_entry.delete(0, END)
      images_name_prefix_label_entry.insert(END, line.split())
      line = f.readline()
      images_name_postfix_label_entry.delete(0, END)
      images_name_postfix_label_entry.insert(END, line.split())

      line = f.readline()
      scale_x_size_label_entry.delete(0, END)
      scale_x_size_label_entry.insert(END, line.split())

      line = f.readline()
      scale_y_size_label_entry.delete(0, END)
      scale_y_size_label_entry.insert(END, line.split())

      line = f.readline()
      scale_z_size_label_entry.delete(0, END)
      scale_z_size_label_entry.insert(END, line.split())

      line = f.readline()
      pt_search_range_label_entry.delete(0, END)
      pt_search_range_label_entry.insert(END, line.split())

show_disp = 0
def change_disp():
  global show_disp
  show_disp = 1 - show_disp
  if (show_disp) :
    button_text = "z"
  else :
    button_text = "uv"
  show_choose_disp_button.config(text=button_text)

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

def load_array_from_file(path, arr_dim) :
  arr_list = []
  with open(path, 'r') as f:
    lines = f.readlines()
    for line in lines :
      vals =[float(j) for j in line.split()]
      arr_list.append(np.array(vals))
  
  return np.squeeze(np.array(arr_list))

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

  
  if(images_name_prefix_label_entry.get() == "") :
    show_error("Images' time prefix is empty!")
    return
  
  # if(images_name_postfix_label_entry.get() == "") :
  #   show_error("Images' slice postfix is empty!")
  #   return
  
  global show_disp
  results_filenames = get_filenames(results_path_entry.get())

  linked3D = pd.read_pickle(results_path_entry.get() + results_filenames[0])

  ref_fn = int(show_ref_time_label_entry.get())
  ev_fn = int(show_def_time_label_entry.get())

  F0 = linked3D[linked3D['frame']==ref_fn].set_index('particle')
  F1 = linked3D[linked3D['frame']==ev_fn].set_index('particle')

  F0 = F0[F0.index.isin(F1.index)].sort_index()
  F1 = F1[F1.index.isin(F0.index)].sort_index()

  image_filenames = get_filenames(images_path_label_entry.get())
  _, image_extension = os.path.splitext(image_filenames[0])
  frames = pims.ImageSequenceND(images_path_label_entry.get()+'*'+ image_extension, axes_identifiers = [images_name_postfix_label_entry.get(), images_name_prefix_label_entry.get()])
  
  is_2d_case = False
  if (images_name_postfix_label_entry.get() == "") :
    is_2d_case = True
  
  if (is_2d_case) :
    frames.bundle_axes = ['x', 'y']
    X0_img = F0[['x', 'y']].to_numpy()
    X1_img = F1[['x', 'y']].to_numpy()
  else :
    frames.bundle_axes = ['x', 'y', 'z']
    X0_img = F0[['x', 'y', 'z']].to_numpy()
    X1_img = F1[['x', 'y', 'z']].to_numpy()
  frames.iter_axes = 'w'

  # Calibration
  scale_x = float(scale_x_size_label_entry.get()) # um / pixel
  scale_y = float(scale_y_size_label_entry.get()) # um / pixel
  scale_z = float(scale_z_size_label_entry.get()) # um / pixel
  px_to_um = np.array([scale_x, scale_y, scale_z])

  X0_um = X0_img * px_to_um
  X1_um = X1_img * px_to_um
  disp = X1_um - X0_um

  ev_slice = int(show_z_label_entry.get())

  if (show_disp) :
    index = 2
  else :
    index = 0 # add u and v separately

  fig0, ax0 = plt.subplots(figsize=(6,3),dpi=300)
  if (is_2d_case) :
    ax0.imshow(frames[ev_fn],cmap='gray')
    sc0 = ax0.scatter(X1_img[0],X1_img[1], c=disp[0]-np.nanmean(disp[0]),cmap='jet',s=1)
  else :
    pt_filter = np.logical_and(X0_img[:,2]>=40, X0_img[:,2]<60)
    ax0.imshow(frames[ev_fn][ev_slice],cmap='gray')
    sc0 = ax0.scatter(X1_img[pt_filter,0],X1_img[pt_filter,1], c=disp[pt_filter,index]-np.nanmean(disp[pt_filter,index]),cmap='jet',s=1)

  sc0.set_clim(-50,50)
  cbar0 = fig0.colorbar(sc0,ax=ax0)
  ax0.set_aspect('equal', adjustable='box')
  ax0.axis('off')

  image_window = Tk()
  image_window.title("PT Result")

  canvas = FigureCanvasTkAgg(fig0, master=image_window)
  canvas.draw()
  canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
  image_window.mainloop()

def sort_by_postfix(filenames):
  """
  Sorts a list of filenames by the postfix before the extension.

  Args:
      filenames: A list of filenames.

  Returns:
      A new list of filenames sorted by postfix.
  """
  def get_postfix(filename):
    name, ext = os.path.splitext(filename)
    parts = name.split(".")  # Split by any dots in the filename (not just extension)
    return int(parts[-1]) if parts[-1].isdigit() else parts[-1]  # Extract last part and convert to int if digit

  return sorted(filenames, key=get_postfix)

def run_pt():
  if(settings_path_entry.get() == "") :
    show_error("Settins path is empty!")
    return
  
  if(settings_name_entry.get() == "") :
    show_error("Settins name is empty!")
    return
  
  print("Running 3D PT...")

  # Calibration
  scale_x = float(scale_x_size_label_entry.get()) # um / pixel
  scale_y = float(scale_y_size_label_entry.get()) # um / pixel
  scale_z = float(scale_z_size_label_entry.get()) # um / pixel
  px_to_um = np.array([scale_x, scale_y, scale_z])

  is_2d_case = False
  if (images_name_postfix_label_entry.get() == "") :
    is_2d_case = True

  image_filenames = get_filenames(images_path_label_entry.get())
  _, image_extension = os.path.splitext(image_filenames[0])

  frames = pims.ImageSequenceND(images_path_label_entry.get()+'*'+ image_extension, axes_identifiers = [images_name_postfix_label_entry.get(), images_name_prefix_label_entry.get()])
  if (is_2d_case) :
    frames.bundle_axes = ['x', 'y']
  else :
    frames.bundle_axes = ['x', 'y', 'z']
  frames.iter_axes = 'w'

  time_size = frames.shape[0]
  # print("times: ", time_size)

  # pt locations
  locations_filenames = get_filenames(locations_path_entry.get())
  locations_filenames = sort_by_postfix(locations_filenames)
  _, image_extension = os.path.splitext(locations_filenames[0])

  pt_locations_list = []
  if (image_extension == '.csv') :
    for location_f in locations_filenames :
      df = pd.read_csv(locations_path_entry.get() + location_f)
      columns_names = df.columns.values.tolist()

      if (is_2d_case) :
        dim = 2
        df = df.rename(columns={columns_names[0]: 'x', columns_names[1]: 'y'})
        df['xum'] = df.x * scale_x
        df['yum'] = df.y * scale_y
        pt_pos_columns = ['xum', 'yum']
      else :
        dim = 3
        df = df.rename(columns={columns_names[0]: 'x', columns_names[1]: 'y', columns_names[2]: 'z'})
        df['xum'] = df['x'] * scale_x
        df['yum'] = df['y'] * scale_y
        df['zum'] = df['z'] * scale_z
        pt_pos_columns = ['xum', 'yum', 'zum']

      pt_locations_list.append(df)
  else :
    show_error("Particles' locations files must have .csv extension!")
    return

  # pt_loc = pd.concat(pt_locations_list)
  # pt_loc.iter_axes = 'frame'
  # print(pt_loc[pt_loc['frame']==1].to_string())
  pt_search_range = int(pt_search_range_label_entry.get())
  disp_interps = []

  dic_results_filenames = get_filenames(dic_results_path_entry.get())
  dic_results_filenames = sort_by_postfix(dic_results_filenames)
  for tt0 in np.arange(time_size-1):
    tt1 = tt0 + 1
    grid_0 = load_array_from_file(dic_results_path_entry.get() + 'ref_'+str(tt0)+r'_'+str(tt1)+'.txt', dim)
    grid_1 = load_array_from_file(dic_results_path_entry.get() + 'def_'+str(tt0)+r'_'+str(tt1)+'.txt', dim)
    corr = load_array_from_file(dic_results_path_entry.get() + 'coefs_'+str(tt0)+r'_'+str(tt1)+'.txt', dim)
    grid_0, grid_1, corr = RemoveWrongMatch(grid_1, grid_0,corr, 10, 0.25, 0.75, corr_th=30, max_ux=500)
    if (is_2d_case) :
      grid_0 = grid_0[:2]
      grid_1 = grid_1[:2]
    grid_0 = grid_0 * px_to_um
    grid_1 = grid_1 * px_to_um
    disp_01 = grid_1 - grid_0
    print("tt0: ", tt0)
    disp_interps.append(LinearNDInterpolatorExt(grid_0, disp_01))

  pt_link_path = results_path_entry.get() + str(pt_search_range) + r'um_corr.pkl'

  @tp.predict.predictor
  def pred01(t1, particle):
    pos_pred_t1 = particle.pos + disp_interps[particle.t](particle.pos)[0]
    return pos_pred_t1
    
  # For data structure
  fTuple = ()
  for pt_loc in pt_locations_list:
    fTuple = fTuple + (pt_loc,)
  # Link particles
  pt_link = pd.concat(tp.link_df_iter(fTuple, search_range = pt_search_range, pos_columns = pt_pos_columns,
                                     adaptive_stop = 0.01, adaptive_step = 0.95, predictor=pred01))
    
  # Save linked
  pt_link.to_pickle(pt_link_path)



# Initialize the main window
root = Tk()
root.title("3D PT")

# Create frames 
text_frame = Frame(root)
text_frame.grid(row=0, column=0, sticky="ew")

settings_name_frame = Frame(root)
settings_name_frame.grid(row=1, column=0, sticky="ew")

images_name_settings_frame = Frame(root)
images_name_settings_frame.grid(row=2, column=0, sticky="ew")

run_frame = Frame(root)
run_frame.grid(row=3, column=0, sticky="ew")

show_frame = Frame(root)
show_frame.grid(row=3, column=1, sticky="ew")


images_path_label = Label(text_frame, text="Images path:")
images_path_label.grid(row=0, column=0, sticky="ew")
images_path_label_entry = Entry(text_frame, width=100)
images_path_label_entry.grid(row=0, column=1, sticky="ew")
images_path_browse_button = Button(text_frame, text="browse", command=browse_directory_images_path)
images_path_browse_button.grid(row=0, column=3, sticky="ew")

dic_results_path = Label(text_frame, text="DIC results path:")
dic_results_path.grid(row=1, column=0, sticky="ew")
dic_results_path_entry = Entry(text_frame, width=100)
dic_results_path_entry.grid(row=1, column=1, sticky="ew")
dic_results_path_browse_button = Button(text_frame, text="browse", command=browse_directory_dic_results_path)
dic_results_path_browse_button.grid(row=1, column=3, sticky="ew")

locations_path = Label(text_frame, text="Particles' locations path:")
locations_path.grid(row=2, column=0, sticky="ew")
locations_path_entry = Entry(text_frame, width=100)
locations_path_entry.grid(row=2, column=1, sticky="ew")
locations_path_browse_button = Button(text_frame, text="browse", command=browse_directory_locations_path)
locations_path_browse_button.grid(row=2, column=3, sticky="ew")

results_path = Label(text_frame, text="Results path:")
results_path.grid(row=3, column=0, sticky="ew")
results_path_entry = Entry(text_frame, width=100)
results_path_entry.grid(row=3, column=1, sticky="ew")
results_path_browse_button = Button(text_frame, text="browse", command=browse_directory_results_path)
results_path_browse_button.grid(row=3, column=3, sticky="ew")

settings_path = Label(text_frame, text="Settings path:")
settings_path.grid(row=4, column=0, sticky="ew")
settings_path_entry = Entry(text_frame, width=100)
settings_path_entry.grid(row=4, column=1, sticky="ew")
settings_path_browse_button = Button(text_frame, text="browse", command=browse_directory_settings_path)
settings_path_browse_button.grid(row=4, column=3, sticky="ew")

settings_name = Label(text_frame, text="Settings filename:")
settings_name.grid(row=5, column=0, sticky="ew")
settings_name_entry = Entry(text_frame, width=30)
settings_name_entry.grid(row=5, column=1, sticky="ew")
settings_name_extension = Label(text_frame, text=".txt")
settings_name_extension.grid(row=5, column=2, sticky="ew")
settings_save_button = Button(text_frame, text="Save", command=save_settings)
settings_save_button.grid(row=5, column=3, sticky="ew")


scale_size_label = Label(settings_name_frame, text="Scale ")
scale_size_label.grid(row=0, column=0, sticky="ew")
scale_x_size_label = Label(settings_name_frame, text="x:")
scale_x_size_label.grid(row=0, column=1, sticky="ew")
scale_x_size_label_entry = Entry(settings_name_frame, width=5)
scale_x_size_label_entry.grid(row=0, column=2, sticky="ew")
scale_y_size_label = Label(settings_name_frame, text="y:")
scale_y_size_label.grid(row=0, column=3, sticky="ew")
scale_y_size_label_entry = Entry(settings_name_frame, width=5)
scale_y_size_label_entry.grid(row=0, column=4, sticky="ew")
scale_z_size_label = Label(settings_name_frame, text="z:")
scale_z_size_label.grid(row=0, column=5, sticky="ew")
scale_z_size_label_entry = Entry(settings_name_frame, width=5)
scale_z_size_label_entry.grid(row=0, column=6, sticky="ew")
scale_size_label_end_label = Label(settings_name_frame, text="um/px")
scale_size_label_end_label.grid(row=0, column=7, sticky="ew")

pt_search_range_label = Label(settings_name_frame, text="PT search range:")
pt_search_range_label.grid(row=1, column=0, sticky="ew")
pt_search_range_label_entry = Entry(settings_name_frame, width=5)
pt_search_range_label_entry.grid(row=1, column=2, sticky="ew")

images_name_prefix_label = Label(images_name_settings_frame, text="Image time prefix:")
images_name_prefix_label.grid(row=0, column=0, sticky="ew")
images_name_prefix_label_entry = Entry(images_name_settings_frame, width=7)
images_name_prefix_label_entry.grid(row=0, column=1, sticky="ew")
images_name_postfix_label = Label(images_name_settings_frame, text="Image slice postfix:")
images_name_postfix_label.grid(row=0, column=2, sticky="ew")
images_name_postfix_label_entry = Entry(images_name_settings_frame, width=7)
images_name_postfix_label_entry.grid(row=0, column=3, sticky="ew")

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


# RUN PT
load_settings_button = Button(run_frame, text="Load Settings", command=load_settings)  # Define browse_pt_image function later
load_settings_button.grid(row=0, column=0, sticky="ew")
run_label_button = Button(run_frame, text="Run PT", command=run_pt)  # Define browse_pt_image function later
run_label_button.grid(row=1, column=0, sticky="ew")

# Run the main loop
root.mainloop()

    

