from tkinter import Tk, Label, Entry, Button, Canvas, Frame, filedialog
from tkinter import *

import os
import re

import cv2
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from PIL import Image, ImageTk

import warnings
warnings.filterwarnings("ignore", module="matplotlib")

from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator


def load_array_from_file(path) :
  arr_list = []
  with open(path, 'r') as f:
    lines = f.readlines()
    for line in lines :
      vals =[float(j) for j in line.split()]
      arr_list.append(np.array(vals))
  
  return np.squeeze(np.array(arr_list))

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


def draw_image_uv_disps(image_path, points_ref, points_def, show_disp, scale = 1., text = None, filename = None) :
    
    assert (points_ref.shape) == (points_def.shape), 'The shape of reference points array must be the same as the shape of deformed points array!'
    
    image = cv2.imread(image_path)
	
    if text is not None :
        image = cv2.putText(image, text, (50,50), cv2.FONT_HERSHEY_SIMPLEX, 1,(255,255,255),4)
		
    if (show_disp) :
      for i, pt0_z in enumerate(points_ref[:, 2]):
          pt1_z = points_def[:, 2][i]
          if np.isnan(pt0_z)==False and np.isnan(pt1_z)==False :
              print("points_ref[i, 0]", points_ref[i, 0])
              x = int(points_ref[i, 0]) #int
              y = int(points_ref[i, 1])
              deltat_z = int(pt1_z - pt0_z)
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


def browse_directory(tkinter_entry):
  directory = filedialog.askdirectory(title="Select a Directory", mustexist=True)
  
  if directory:
    tkinter_entry.delete(0, END)
    tkinter_entry.insert(END,directory + '/')

def show_error(message) :
  error_window = Tk()
  error_window.title("Error")
  error_window.geometry("200x100")
  error_label = Label(error_window, text=message)
  error_label.grid(row=0, column=0, sticky="ew")
  error_window.mainloop()

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
