from helpers import *
import subprocess
import pandas as pd
import trackpy as tp
import pims

import ast

import json
  

# Initialize the main window
root = Tk()
root.title("3D DIC")

# ignore_1st_layer, calculation_order, 2D_or_3D, show_uv_or_z, show_dic_or_pt
buttons_vals = [0, 0, 0, 0, 0]

def change_backward():
  buttons_vals[1] = 1 - buttons_vals[1]
  if (buttons_vals[1]) :
    button_text = "Backward"
  else :
    button_text = "Forward"
  backward_value_button.config(text=button_text)

def change_case():
  buttons_vals[2] = 1 - buttons_vals[2]
  if (buttons_vals[2]) :
    button_text = "2D"
  else :
    button_text = "3D"
  case_value_button.config(text=button_text)
#   update_GUI_fields()

def change_disp():
  buttons_vals[3] = 1 - buttons_vals[3]
  if (buttons_vals[3]) :
    button_text = "z"
  else :
    button_text = "uv"
  show_choose_disp_button.config(text=button_text)

def show_dic_or_pt():
  buttons_vals[4] = 1 - buttons_vals[4]
  if (buttons_vals[4]) :
    button_text = "PT"
  else :
    button_text = "DIC"
  show_dic_or_pt_button.config(text=button_text)

def load_dic_settings():
  filetypes = (
        ('json', '*.json'),
        ('All files', '*.*')
    )
  directory = filedialog.askopenfilename(title="Select a File", filetypes=filetypes)
  
  # file_name = os.path.splitext(os.path.basename(directory))[0]
  path_file = os.path.split(os.path.abspath(directory))

  dic_settings_path_entry.delete(0, END)
  dic_settings_path_entry.insert(END, path_file[0] + '/')
  dic_settings_name_entry.delete(0, END)
  dic_settings_name_entry.insert(END, os.path.splitext(path_file[1])[0])

  if directory:
    with open(directory, 'r') as f:
      json_dic_set = json.load(f)

      subset_size_label_entry.delete(0, END)
      subset_size_label_entry.insert(END, str(json_dic_set['Subset size']))

      subset_offset_label_entry.delete(0, END)
      subset_offset_label_entry.insert(END,  str(json_dic_set['Subset offset']))

      z_spacing_label_entry.delete(0, END)
      z_spacing_label_entry.insert(END,  str(json_dic_set['Z spacing']))

      z_search_label_entry.delete(0, END)
      z_search_label_entry.insert(END,  str(json_dic_set['Z search']))

      images_path_label_entry.delete(0, END)
      images_path_label_entry.insert(END, json_dic_set['Images path'])

      dic_results_path_entry.delete(0, END)
      dic_results_path_entry.insert(END, json_dic_set['DIC Results'])

      images_name_prefix_label_entry.delete(0, END)
      images_name_prefix_label_entry.insert(END, json_dic_set['Images Time Prefix'])

      images_name_postfix_label_entry.delete(0, END)
      images_name_postfix_label_entry.insert(END, json_dic_set['Images Slice Postfix'])
      
      ROI_x_min_label_entry.delete(0, END)
      ROI_x_min_label_entry.insert(END, str(json_dic_set['ROI x min']))
      ROI_y_min_label_entry.delete(0, END)
      ROI_y_min_label_entry.insert(END, str(json_dic_set['ROI y min']))
      ROI_x_max_label_entry.delete(0, END)
      ROI_x_max_label_entry.insert(END, str(json_dic_set['ROI x max']))
      ROI_y_max_label_entry.delete(0, END)
      ROI_y_max_label_entry.insert(END, str(json_dic_set['ROI y max']))
      
      downsampling_label_entry.delete(0, END)
      downsampling_label_entry.insert(END, str(json_dic_set['Downsampling']))

      threshold_label_entry.delete(0, END)
      threshold_label_entry.insert(END, str(json_dic_set['coef threshold']))

      delta_threshold_label_entry.delete(0, END)
      delta_threshold_label_entry.insert(END, str(json_dic_set['delta coef threshold']))

      crack_gradient_label_entry.delete(0, END)
      crack_gradient_label_entry.insert(END, str(json_dic_set['crack gradient threshold']))

      buttons_vals[1] = int(json_dic_set['Is Backward'])
      if (buttons_vals[1]) :
        button_text = "Backward"
      else :
        button_text = "Forward"
      backward_value_button.config(text=button_text)

      buttons_vals[2] = int(json_dic_set['Is 2D'])
      if (buttons_vals[2]) :
        button_text = "2D"
      else :
        button_text = "3D"
      case_value_button.config(text=button_text)

      # if 2D
      if (buttons_vals[2]) :
        z_spacing_label_entry.delete(0, END)
        z_search_label_entry.delete(0, END)

    f.close()

def load_pt_settings():
  filetypes = (
        ('json', '*.json'),
        ('All files', '*.*')
    )
  directory = filedialog.askopenfilename(title="Select a File", filetypes=filetypes)
  
  # file_name = os.path.splitext(os.path.basename(directory))[0]
  path_file = os.path.split(os.path.abspath(directory))

  pt_settings_path_entry.delete(0, END)
  pt_settings_path_entry.insert(END, path_file[0] + '/')
  pt_settings_name_entry.delete(0, END)
  pt_settings_name_entry.insert(END, os.path.splitext(path_file[1])[0])

  if directory:
    with open(directory, 'r') as f:
      json_dic_set = json.load(f)
      
      images_path_label_entry.delete(0, END)
      images_path_label_entry.insert(END, json_dic_set['Images path'])

      dic_results_path_entry.delete(0, END)
      dic_results_path_entry.insert(END, json_dic_set['DIC Results'])

      locations_path_entry.delete(0, END)
      locations_path_entry.insert(END, json_dic_set['Locations path'])

      pt_results_path_entry.delete(0, END)
      pt_results_path_entry.insert(END, json_dic_set['PT Results'])

      images_name_prefix_label_entry.delete(0, END)
      images_name_prefix_label_entry.insert(END, json_dic_set['Images Time Prefix'])

      images_name_postfix_label_entry.delete(0, END)
      images_name_postfix_label_entry.insert(END, json_dic_set['Images Slice Postfix'])

      scale_x_size_label_entry.delete(0, END)
      scale_x_size_label_entry.insert(END, str(json_dic_set['Scale x']))
      scale_y_size_label_entry.delete(0, END)
      scale_y_size_label_entry.insert(END, str(json_dic_set['Scale y']))
      scale_z_size_label_entry.delete(0, END)
      scale_z_size_label_entry.insert(END, str(json_dic_set['Scale z']))

      pt_search_range_label_entry.delete(0, END)
      pt_search_range_label_entry.insert(END, str(json_dic_set['PT search range']))

      buttons_vals[1] = int(json_dic_set['Is Backward'])
      if (buttons_vals[1]) :
        button_text = "Backward"
      else :
        button_text = "Forward"
      backward_value_button.config(text=button_text)

      buttons_vals[2] = int(json_dic_set['Is 2D'])
      if (buttons_vals[2]) :
        button_text = "2D"
      else :
        button_text = "3D"
      case_value_button.config(text=button_text)

      # if 2D
      if (buttons_vals[2]) :
        scale_z_size_label_entry.delete(0, END)

    f.close()
      


def dic_save_settings():
  print("Saving current settings...")

  if(images_path_label_entry.get() == "") :
    show_error("Images path is empty!")
    return
  if(dic_results_path_entry.get() == "") :
    show_error("Results path is empty!")
    return
  if(dic_settings_path_entry.get() == "") :
    show_error("Settins path is empty!")
    return
  if(dic_settings_name_entry.get() == "") :
    show_error("Settins name is empty!")
    return
  if(subset_size_label_entry.get() == "") :
    show_error("Subset size is empty!")
    return
  if(subset_offset_label_entry.get() == "") :
    show_error("Subset offset is empty!")
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
  if(threshold_label_entry.get() == "") :
    show_error("coef threshold is empty!")
    return
  if(delta_threshold_label_entry.get() == "") :
    show_error("delta coef threshold is empty!")
    return
  if(crack_gradient_label_entry.get() == "") :
    show_error("crack gradient threshold is empty!")
    return
  
  image_filenames = get_filenames(images_path_label_entry.get())
  _, image_extension = os.path.splitext(image_filenames[0])
  times = find_numbers_in_filenames(image_filenames, images_name_prefix_label_entry.get())
  times = sorted(times)

  image0 = cv2.imread(images_path_label_entry.get() + image_filenames[0])

  if not (buttons_vals[2]):
    if(z_search_label_entry.get() == "") :
      show_error("z-search radius is empty!")
      return
    if(z_spacing_label_entry.get() == "") :
      show_error("z-spacing is empty!")
      return
    if(images_name_postfix_label_entry.get() == "") :
      show_error("Images' name postfix is empty!")
      return
    
    stack_hs = find_numbers_in_filenames(image_filenames, images_name_prefix_label_entry.get() + str(list(times)[0]) + images_name_postfix_label_entry.get())
    stack_h = max(list(stack_hs))
    z_spacing = int(z_spacing_label_entry.get())
    z_search = int(z_search_label_entry.get())
  else :
    stack_h = 0
    z_spacing = -1
    z_search = -1

  start_layer = 0
  if (calc_layers_from_entry.get() != "") :
    start_layer = int(calc_layers_from_entry.get())
  end_layer = stack_h
  if (calc_layers_to_entry.get() != "") :
    end_layer = int(calc_layers_to_entry.get())
  layers_to_calculate = list(range(start_layer, end_layer + 1))

  file = dic_settings_path_entry.get() + dic_settings_name_entry.get() + ".json"
  with open(file, 'w') as f:
    json.dump({
      "Image size x": int(image0.shape[1]),
      "Image size y": int(image0.shape[0]),
      "Subset size": int(subset_size_label_entry.get()), 
      "Subset offset": int(subset_offset_label_entry.get()), 
      "Z spacing": z_spacing, 
      "Z search": z_search, 
      "Images path": images_path_label_entry.get(), 
      "DIC Results": dic_results_path_entry.get(), 
      "Images Time Prefix": images_name_prefix_label_entry.get(),   
      "Images Slice Postfix": images_name_postfix_label_entry.get(), 
      "Times": times, 
      "ROI x min": int(ROI_x_min_label_entry.get()), 
      "ROI y min": int(ROI_y_min_label_entry.get()), 
      "ROI x max": int(ROI_x_max_label_entry.get()), 
      "ROI y max": int(ROI_y_max_label_entry.get()), 
      "Downsampling": int(downsampling_label_entry.get()), 
      "Stack height": int(stack_h), 
      "Image extension": image_extension, 
      "Is 2D" : int(buttons_vals[2]), 
      "Is Backward" : int(buttons_vals[1]), 
      "layers to calculate" : layers_to_calculate,
      "coef threshold" : float(threshold_label_entry.get()),
      "delta coef threshold" : float(delta_threshold_label_entry.get()),
      "crack gradient threshold" : float(crack_gradient_label_entry.get())
    }, f, indent=2)
  
  f.close()

def pt_save_settings():
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
  if(pt_results_path_entry.get() == "") :
    show_error("Results path is empty!")
    return
  if(pt_settings_path_entry.get() == "") :
    show_error("Settins path is empty!")
    return
  if(pt_settings_name_entry.get() == "") :
    show_error("Settins name is empty!")
    return
  if(scale_x_size_label_entry.get() == "") :
    show_error("Scale x is empty!")
    return
  if(scale_y_size_label_entry.get() == "") :
    show_error("Scale y is empty!")
    return
  if(pt_search_range_label_entry.get() == "") :
    show_error("PT seacrh range is empty!")
    return
  if(images_name_prefix_label_entry.get() == "") :
    show_error("Images' name prefix is empty!")
    return
  #if not 2D
  if not (buttons_vals[2]) :
    if(scale_z_size_label_entry.get() == "") :
      show_error("Scale z is empty!")
      return
    if(images_name_postfix_label_entry.get() == "") :
      show_error("Images' name postfix is empty!")
      return
    scale_z = int(scale_z_size_label_entry.get())
  else :
    scale_z = -1

  file = pt_settings_path_entry.get() + pt_settings_name_entry.get() + ".json"
  with open(file, 'w') as f:
    json.dump({
      "Images path": images_path_label_entry.get(),
      "DIC Results": dic_results_path_entry.get(),
      "Locations path": locations_path_entry.get(),
      "PT Results": pt_results_path_entry.get(),
      "Images Time Prefix": images_name_prefix_label_entry.get(),   
      "Images Slice Postfix": images_name_postfix_label_entry.get(),
      "Scale x": int(scale_x_size_label_entry.get()),
      "Scale y": int(scale_y_size_label_entry.get()),
      "Scale z": scale_z,
      "PT search range": int(pt_search_range_label_entry.get()),
      "Is 2D" : int(buttons_vals[2]), 
      "Is Backward" : int(buttons_vals[1]), 
    }, f, indent=2)
  f.close()

def run_dic():
  dic_save_settings()
  
  print("Running 3D DIC...")
  setting_file = dic_settings_path_entry.get() + dic_settings_name_entry.get() + ".json"
  p = subprocess.Popen(["./DIC_3D", setting_file])

def show_disps_dic():
  if(images_path_label_entry.get() == "") :
    show_error("Images path is empty!")
    return
  if(dic_results_path_entry.get() == "") :
    show_error("Results path is empty!")
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
  if(subset_size_label_entry.get() == "") :
    show_error("Subset size is empty!")
    return
  if(subset_offset_label_entry.get() == "") :
    show_error("Subset offset is empty!")
    return
  if(images_name_prefix_label_entry.get() == "") :
    show_error("Images' name prefix is empty!")
    return
  
  image_filenames = get_filenames(images_path_label_entry.get())
  _, image_extension = os.path.splitext(image_filenames[0])
  times = find_numbers_in_filenames(image_filenames, images_name_prefix_label_entry.get())
  #if not 2D
  if not buttons_vals[2] :
    if(z_spacing_label_entry.get() == "") :
      show_error("z-spacing is empty!")
      return
    if(show_z_label_entry.get() == "") :
      show_error("Select z to show!")
      return
    if(images_name_postfix_label_entry.get() == "" and buttons_vals[2] != 0) :
      show_error("Images' name postfix is empty!")
      return
    stack_hs = find_numbers_in_filenames(image_filenames, images_name_prefix_label_entry.get() + str(list(times)[0]) + images_name_postfix_label_entry.get())
    stack_h = max(list(stack_hs))
    index_to_show = int(show_z_label_entry.get())
    z_bounce = int(z_spacing_label_entry.get())
    bounced_z_id = index_to_show // z_bounce
    s_z = len(range(0, stack_h, z_bounce))
    ref_image_path = images_path_label_entry.get() + images_name_prefix_label_entry.get() + str(show_ref_time_label_entry.get()) + images_name_postfix_label_entry.get() + '{:03}'.format(int(show_z_label_entry.get())) + image_extension
    show_case = buttons_vals[3]
  else :
    stack_h = 1
    index_to_show = 0
    z_bounce = 1
    bounced_z_id = 0
    s_z = 1
    ref_image_path = images_path_label_entry.get() + images_name_prefix_label_entry.get() + str(show_ref_time_label_entry.get()) + image_extension
    show_case = 0

  subset_size, subset_offset = int(subset_size_label_entry.get()), int(subset_offset_label_entry.get())
  roi_xy_min = [int(ROI_x_min_label_entry.get()), int(ROI_y_min_label_entry.get())]
  roi_xy_max = [int(ROI_x_max_label_entry.get()), int(ROI_y_max_label_entry.get())]

  ref_res_path = dic_results_path_entry.get() + 'ref_' + str(show_ref_time_label_entry.get()) + '_' + str(show_def_time_label_entry.get()) + '.txt'
  def_res_path = dic_results_path_entry.get() + 'def_' + str(show_ref_time_label_entry.get()) + '_' + str(show_def_time_label_entry.get()) + '.txt'
  # coefs_res_path = results_path_entry.get() + 'coefs_' + str(show_ref_time_label_entry.get()) + '_' + str(show_def_time_label_entry.get()) + '.txt'

  result_ref = load_array_from_file(ref_res_path)
  result_def = load_array_from_file(def_res_path)

  if ((index_to_show % z_bounce) != 0) :
    show_error("Result for chosen z doesn't exist!")
    return
  
  points_ref, points_def = result_ref[bounced_z_id::s_z], result_def[bounced_z_id::s_z]
  draw_image_uv_disps(ref_image_path, points_ref, points_def, show_case)


def show_disps_pt():
  if(images_path_label_entry.get() == "") :
    show_error("Images path is empty!")
    return
  if(pt_results_path_entry.get() == "") :
    show_error("Results path is empty!")
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
  #if not 2D
  
  results_filenames = get_filenames(pt_results_path_entry.get())

  linked3D = pd.read_pickle(pt_results_path_entry.get() + results_filenames[0])

  ref_fn = int(show_ref_time_label_entry.get())
  ev_fn = int(show_def_time_label_entry.get())

  F0 = linked3D[linked3D['frame']==ref_fn].set_index('particle')
  F1 = linked3D[linked3D['frame']==ev_fn].set_index('particle')

  F0 = F0[F0.index.isin(F1.index)].sort_index()
  F1 = F1[F1.index.isin(F0.index)].sort_index()

  # Calibration
  scale_x = float(scale_x_size_label_entry.get()) # um / pixel
  scale_y = float(scale_y_size_label_entry.get()) # um / pixel

  image_filenames = get_filenames(images_path_label_entry.get())
  _, image_extension = os.path.splitext(image_filenames[0])
  if not buttons_vals[2] :
    if(show_z_label_entry.get() == "") :
      show_error("Select z to show!")
      return
    if(images_name_postfix_label_entry.get() == "" and buttons_vals[2] != 0) :
      show_error("Images' name postfix is empty!")
      return
    X0_img = F0[['x', 'y', 'z']].to_numpy()
    X1_img = F1[['x', 'y', 'z']].to_numpy()
    image_postfix = images_name_postfix_label_entry.get().lstrip('_')
    frames = pims.ImageSequenceND(images_path_label_entry.get()+'*'+ image_extension, axes_identifiers = [image_postfix, images_name_prefix_label_entry.get()])
    frames.bundle_axes = ['x', 'y', 'z']
    px_to_um = np.array([scale_x, scale_y, float(scale_z_size_label_entry.get())])
    show_component = buttons_vals[3] + 1 # y or z
  else :
    X0_img = F0[['x', 'y']].to_numpy()
    X1_img = F1[['x', 'y']].to_numpy()
    frames = pims.ImageSequenceND(images_path_label_entry.get()+'*'+ image_extension, axes_identifiers = [images_name_prefix_label_entry.get()])
    frames.bundle_axes = ['x', 'y']
    px_to_um = np.array([scale_x, scale_y])
    show_component = 1 # y

  frames.iter_axes = 'w'

  X0_um = X0_img * px_to_um
  X1_um = X1_img * px_to_um
  disp = X1_um - X0_um

  fig0, ax0 = plt.subplots(figsize=(6,3), dpi = 300)
  if (buttons_vals[2]) :
    ax0.imshow(frames[ev_fn],cmap='gray')
    sc0 = ax0.scatter(X1_img[0],X1_img[1], c=disp[0]-np.nanmean(disp[0]),cmap='jet',s=1)
  else :
    ev_slice = int(show_z_label_entry.get())
    pt_filter = np.logical_and(X0_img[:,2]>=40, X0_img[:,2]<60)
    ax0.imshow(frames[ev_fn][ev_slice],cmap='gray')
    sc0 = ax0.scatter(X1_img[pt_filter,0],X1_img[pt_filter,1], c=disp[pt_filter,show_component]-np.nanmean(disp[pt_filter,show_component]),cmap='jet',s=1)

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

def generate_locations_window() :
  generate_window = Tk()
  generate_window.title("Generate PT locations")
  # generate_window.geometry("200x100")
  generate_label = Label(generate_window, text="No PT locations found! Do you want to generate locations?")
  generate_label.grid(row=0, column=0, sticky="ew")

  diameter_label = Label(generate_window, text="Diameter ")
  diameter_label.grid(row=1, column=0, sticky="ew")
  diameter_x_size_label = Label(generate_window, text="x:")
  diameter_x_size_label.grid(row=1, column=1, sticky="ew")
  diameter_x_size_label_entry = Entry(generate_window, width=5)
  diameter_x_size_label_entry.grid(row=1, column=2, sticky="ew")
  diameter_y_size_label = Label(generate_window, text="y:")
  diameter_y_size_label.grid(row=1, column=3, sticky="ew")
  diameter_y_size_label_entry = Entry(generate_window, width=5)
  diameter_y_size_label_entry.grid(row=1, column=4, sticky="ew")
  diameter_z_size_label = Label(generate_window, text="z:")
  diameter_z_size_label.grid(row=1, column=5, sticky="ew")
  diameter_z_size_label_entry = Entry(generate_window, width=5)
  diameter_z_size_label_entry.grid(row=1, column=6, sticky="ew")

  separation_label = Label(generate_window, text="Separation ")
  separation_label.grid(row=2, column=0, sticky="ew")
  separation_x_size_label = Label(generate_window, text="x:")
  separation_x_size_label.grid(row=2, column=1, sticky="ew")
  separation_x_size_label_entry = Entry(generate_window, width=5)
  separation_x_size_label_entry.grid(row=2, column=2, sticky="ew")
  separation_y_size_label = Label(generate_window, text="y:")
  separation_y_size_label.grid(row=2, column=3, sticky="ew")
  separation_y_size_label_entry = Entry(generate_window, width=5)
  separation_y_size_label_entry.grid(row=2, column=4, sticky="ew")
  separation_z_size_label = Label(generate_window, text="z:")
  separation_z_size_label.grid(row=2, column=5, sticky="ew")
  separation_z_size_label_entry = Entry(generate_window, width=5)
  separation_z_size_label_entry.grid(row=2, column=6, sticky="ew")

  locations_args = Label(generate_window, text="Optinal args")
  locations_args.grid(row=3, column=0, sticky="ew")
  locations_args_entry = Entry(generate_window, width=100)
  locations_args_entry.grid(row=3, column=7, sticky="ew")

  generate_button = Button(generate_window, text="Generate", command=lambda: generate_lc( \
    [int(diameter_x_size_label_entry.get()), int(diameter_y_size_label_entry.get()), int(diameter_z_size_label_entry.get())], \
    [int(separation_x_size_label_entry.get()), int(separation_y_size_label_entry.get()), int(separation_z_size_label_entry.get())], locations_args_entry.get()))
  generate_button.grid(row=4, column=0, padx=10, sticky="ew")
  generate_window.mainloop()

def generate_lc(diameter, separation, add_args_str) :

  image_filenames = get_filenames(images_path_label_entry.get())
  _, image_extension = os.path.splitext(image_filenames[0])
  image_postfix = images_name_postfix_label_entry.get().lstrip('_')
  frames = pims.ImageSequenceND(images_path_label_entry.get()+'*'+ image_extension, axes_identifiers = [image_postfix, images_name_prefix_label_entry.get()])
  
  expr = ast.parse(f"dict({add_args_str}\n)", mode="eval")
  add_args = {kw.arg: ast.literal_eval(kw.value) for kw in expr.body.keywords}
  print("add_args", add_args)
  print("frames shape:", frames.shape[0])
  pt_loc = pd.DataFrame()
  for tt in range(0, frames.shape[0]):
    print('Locating stack: ', tt)
    pt_loc_temp = tp.locate(frames[tt], diameter=diameter, separation=separation, **add_args)
    pt_loc_temp['frame'] = tt
    pt_loc = pd.concat([pt_loc_temp, pt_loc], ignore_index=True)
  
  pt_loc.to_pickle(locations_path_entry.get() + ''.join(str(p_dim) for p_dim in diameter) + r'_' + ''.join(str(p_sep) for p_sep in separation) + r'.pkl')


def run_pt():
  pt_save_settings()

  # Calibration
  scale_x = float(scale_x_size_label_entry.get()) # um / pixel
  scale_y = float(scale_y_size_label_entry.get()) # um / pixel
  # if not 2D
  if not (buttons_vals[2]) :
    scale_z = float(scale_z_size_label_entry.get()) # um / pixel
    px_to_um = np.array([scale_x, scale_y, scale_z])
  else :
    px_to_um = np.array([scale_x, scale_y])

  # TODO add support of 2D PT

  image_filenames = get_filenames(images_path_label_entry.get())
  _, image_extension = os.path.splitext(image_filenames[0])
  image_postfix = images_name_postfix_label_entry.get().lstrip('_')
  frames = pims.ImageSequenceND(images_path_label_entry.get()+'*'+ image_extension, axes_identifiers = [image_postfix, images_name_prefix_label_entry.get()])
  
  if (buttons_vals[2]) :
    frames.bundle_axes = ['x', 'y']
  else :
    frames.bundle_axes = ['x', 'y', 'z']
  frames.iter_axes = 'w'

  time_size = frames.shape[0]

  # pt locations
  locations_filenames = get_filenames(locations_path_entry.get())
  if (len(locations_filenames) == 0) :
    generate_locations_window()
    return

  locations_filenames = sort_by_postfix(locations_filenames)
  _, image_extension = os.path.splitext(locations_filenames[0])

  pt_loc = pd.read_pickle(locations_path_entry.get() + locations_filenames[0])

  if (buttons_vals[2]) :
    dim = 2
    pt_loc['xum'] = pt_loc['x'] * scale_x
    pt_loc['yum'] = pt_loc['y'] * scale_y
    pt_pos_columns = ['xum', 'yum']
  else :
    dim = 3
    pt_loc['xum'] = pt_loc['x'] * scale_x
    pt_loc['yum'] = pt_loc['y'] * scale_y
    pt_loc['zum'] = pt_loc['z'] * scale_z
    pt_pos_columns = ['xum', 'yum', 'zum']

  # pt_locations_list = []
  # if (image_extension == '.csv') :
  #   for location_f in locations_filenames :
  #     df = pd.read_csv(locations_path_entry.get() + location_f)
  #     columns_names = df.columns.values.tolist()

  #     if (buttons_vals[2]) :
  #       dim = 2
  #       df = df.rename(columns={columns_names[0]: 'x', columns_names[1]: 'y'})
  #       df['xum'] = df.x * scale_x
  #       df['yum'] = df.y * scale_y
  #       pt_pos_columns = ['xum', 'yum']
  #     else :
  #       dim = 3
  #       df = df.rename(columns={columns_names[0]: 'x', columns_names[1]: 'y', columns_names[2]: 'z'})
  #       df['xum'] = df['x'] * scale_x
  #       df['yum'] = df['y'] * scale_y
  #       df['zum'] = df['z'] * scale_z
  #       pt_pos_columns = ['xum', 'yum', 'zum']

  #     pt_locations_list.append(df)
  # else :
  #   show_error("Particles' locations files must have .csv extension!")
  #   return

  # pt_loc = pd.concat(pt_locations_list)
  # pt_loc.iter_axes = 'frame'
  # print(pt_loc[pt_loc['frame']==1].to_string())
  pt_search_range = int(pt_search_range_label_entry.get())
  disp_interps = []

  print("Running 3D PT...")

  dic_results_filenames = get_filenames(dic_results_path_entry.get())
  dic_results_filenames = sort_by_postfix(dic_results_filenames)
  for tt0 in np.arange(time_size-1):
    tt1 = tt0 + 1
    grid_0 = load_array_from_file(dic_results_path_entry.get() + 'ref_'+str(tt0)+r'_'+str(tt1)+'.txt')
    grid_1 = load_array_from_file(dic_results_path_entry.get() + 'def_'+str(tt0)+r'_'+str(tt1)+'.txt')
    corr = load_array_from_file(dic_results_path_entry.get() + 'coefs_'+str(tt0)+r'_'+str(tt1)+'.txt')
    grid_0, grid_1, corr = RemoveWrongMatch(grid_1, grid_0, corr, 10, 0.25, 0.75, corr_th=30, max_ux=500)
    grid_0 = grid_0[:, :dim]
    grid_1 = grid_1[:, :dim]
    grid_0 = grid_0 * px_to_um
    grid_1 = grid_1 * px_to_um
    disp_01 = grid_1 - grid_0
    disp_interps.append(LinearNDInterpolatorExt(grid_0, disp_01))

  pt_link_path = pt_results_path_entry.get() + str(pt_search_range) + r'um_corr.pkl'

  @tp.predict.predictor
  def pred01(t1, particle):
    pos_pred_t1 = particle.pos + disp_interps[particle.t](particle.pos)[0]
    return pos_pred_t1
  
  # For data structure
  fTuple = ()
  for fid in np.arange(np.max(pt_loc.frame)+1):
        fTuple = fTuple + (pt_loc[pt_loc['frame']==fid],)

  # Link particles
  pt_link = pd.concat(tp.link_df_iter(fTuple, search_range = pt_search_range, pos_columns = pt_pos_columns,
                                     adaptive_stop = 0.01, adaptive_step = 0.95, predictor=pred01))
    
  # Save linked
  pt_link.to_pickle(pt_link_path)
  print("PT finished! Check results in ", pt_link_path)


def show_disps_dic_or_pt() :
  if (buttons_vals[4]) :
    show_disps_pt()
  else :
    show_disps_dic()

    

# Create frames 
text_frame = Frame(root)
text_frame.grid(row=0, column=0, sticky="ew")

case_frame = Frame(root)
case_frame.grid(row=0, column=4, sticky="ew")

images_name_settings_frame = Frame(root)
images_name_settings_frame.grid(row=1, column=0, pady=20, sticky="ew")

ROI_frame = Frame(root)
ROI_frame.grid(row=2, column=0, padx=30, sticky="ew")
dic_settings_frame = Frame(root)
dic_settings_frame.grid(row=3, column=0, sticky="ew")

pt_settings_name_frame = Frame(root)
pt_settings_name_frame.grid(row=3, column=1, sticky="ew")

buttons_frame = Frame(root)
buttons_frame.grid(row=5, column=0, sticky="ew")

run_frame = Frame(root)
run_frame.grid(row=2, column=2, sticky="ew")

show_frame = Frame(root)
show_frame.grid(row=4, column=1, sticky="ew")

# Create labels and text entry fields in the text frame
images_path_label = Label(text_frame, text="Images path:")
images_path_label.grid(row=0, column=0, sticky="ew")
images_path_label_entry = Entry(text_frame, width=100)
images_path_label_entry.grid(row=0, column=1, sticky="ew")
images_path_browse_button = Button(text_frame, text="browse", command=lambda: browse_directory(images_path_label_entry))
images_path_browse_button.grid(row=0, column=2, sticky="ew")

dic_results_path = Label(text_frame, text="DIC results path:")
dic_results_path.grid(row=1, column=0, sticky="ew")
dic_results_path_entry = Entry(text_frame, width=100)
dic_results_path_entry.grid(row=1, column=1, sticky="ew")
dic_results_path_browse_button = Button(text_frame, text="browse", command=lambda: browse_directory(dic_results_path_entry))
dic_results_path_browse_button.grid(row=1, column=2, sticky="ew")

dic_settings_path = Label(text_frame, text="DIC settings path:")
dic_settings_path.grid(row=2, column=0, sticky="ew")
dic_settings_path_entry = Entry(text_frame, width=100)
dic_settings_path_entry.grid(row=2, column=1, sticky="ew")
dic_settings_path_browse_button = Button(text_frame, text="browse", command=lambda: browse_directory(dic_settings_path_entry))
dic_settings_path_browse_button.grid(row=2, column=2, sticky="ew")
dic_settings_name = Label(text_frame, text="DIC settings filename:")
dic_settings_name.grid(row=2, column=3, sticky="ew")
dic_settings_name_entry = Entry(text_frame, width=30)
dic_settings_name_entry.grid(row=2, column=4, sticky="ew")
dic_settings_name_extension = Label(text_frame, text=".json")
dic_settings_name_extension.grid(row=2, column=5, sticky="ew")
settings_save_button = Button(text_frame, text="Save", command=dic_save_settings)
settings_save_button.grid(row=2, column=6, sticky="ew")
load_dic_settings_button = Button(text_frame, text="Load DIC settings", command=load_dic_settings)  # Define browse_pt_image function later
load_dic_settings_button.grid(row=2, column=7, sticky="ew")

locations_path = Label(text_frame, text="Particles' locations path:")
locations_path.grid(row=3, column=0, sticky="ew")
locations_path_entry = Entry(text_frame, width=100)
locations_path_entry.grid(row=3, column=1, sticky="ew")
locations_path_browse_button = Button(text_frame, text="browse", command=lambda: browse_directory(locations_path_entry))
locations_path_browse_button.grid(row=3, column=2, sticky="ew")

pt_results_path = Label(text_frame, text="PT results path:")
pt_results_path.grid(row=4, column=0, sticky="ew")
pt_results_path_entry = Entry(text_frame, width=100)
pt_results_path_entry.grid(row=4, column=1, sticky="ew")
pt_results_path_browse_button = Button(text_frame, text="browse", command=lambda: browse_directory(pt_results_path_entry))
pt_results_path_browse_button.grid(row=4, column=2, sticky="ew")


pt_settings_path = Label(text_frame, text="PT settings path:")
pt_settings_path.grid(row=5, column=0, sticky="ew")
pt_settings_path_entry = Entry(text_frame, width=100)
pt_settings_path_entry.grid(row=5, column=1, sticky="ew")
pt_settings_path_browse_button = Button(text_frame, text="browse", command=lambda: browse_directory(pt_settings_path_entry))
pt_settings_path_browse_button.grid(row=5, column=2, sticky="ew")
pt_settings_name = Label(text_frame, text="DIC settings filename:")
pt_settings_name = Label(text_frame, text="PT settings filename:")
pt_settings_name.grid(row=5, column=3, sticky="ew")
pt_settings_name_entry = Entry(text_frame, width=30)
pt_settings_name_entry.grid(row=5, column=4, sticky="ew")
pt_settings_name_extension = Label(text_frame, text=".json")
pt_settings_name_extension.grid(row=5, column=5, sticky="ew")
settings_save_button = Button(text_frame, text="Save", command=pt_save_settings)
settings_save_button.grid(row=5, column=6, sticky="ew")
load_settings_button = Button(text_frame, text="Load PT settings", command=load_pt_settings)  # Define browse_pt_image function later
load_settings_button.grid(row=5, column=7, sticky="ew")






#DIC
subset_size_label = Label(dic_settings_frame, text="Subset size:")
subset_size_label.grid(row=0, column=0, sticky="ew")
subset_size_label_entry = Entry(dic_settings_frame, width=5)
subset_size_label_entry.grid(row=0, column=1, sticky="ew")
subset_size_end_label = Label(dic_settings_frame, text="px")
subset_size_end_label.grid(row=0, column=2, sticky="ew")

subset_offset_label = Label(dic_settings_frame, text="Subset offset:")
subset_offset_label.grid(row=1, column=0, sticky="ew")
subset_offset_label_entry = Entry(dic_settings_frame, width=5)
subset_offset_label_entry.grid(row=1, column=1, sticky="ew")
subset_offset_end_label = Label(dic_settings_frame, text="px")
subset_offset_end_label.grid(row=1, column=2, sticky="ew")

z_search_label = Label(dic_settings_frame, text="z-search radius:")
z_search_label.grid(row=2, column=0, sticky="ew")
z_search_label_entry = Entry(dic_settings_frame, width=5)
z_search_label_entry.grid(row=2, column=1, sticky="ew")
z_search_end_label = Label(dic_settings_frame, text="px")
z_search_end_label.grid(row=2, column=2, sticky="ew")

z_spacing_label = Label(dic_settings_frame, text="z-spacing:")
z_spacing_label.grid(row=3, column=0, sticky="ew")
z_spacing_label_entry = Entry(dic_settings_frame, width=5)
z_spacing_label_entry.grid(row=3, column=1, sticky="ew")
z_spacing_end_label = Label(dic_settings_frame, text="px")
z_spacing_end_label.grid(row=3, column=2, sticky="ew")

downsampling_label = Label(dic_settings_frame, text="downsampling:")
downsampling_label.grid(row=4, column=0, sticky="ew")
downsampling_label_entry = Entry(dic_settings_frame, width=5)
downsampling_label_entry.grid(row=4, column=1, sticky="ew")
downsampling_end_label = Label(dic_settings_frame, text="px")
downsampling_end_label.grid(row=4, column=2, sticky="ew")

calculate_layers1_label = Label(dic_settings_frame, text="calculate layers from:")
calculate_layers1_label.grid(row=5, column=0, sticky="ew")
calc_layers_from_entry = Entry(dic_settings_frame, width=5)
calc_layers_from_entry.grid(row=5, column=1, sticky="ew")
calculate_layers2_label = Label(dic_settings_frame, text="calculate layers to:")
calculate_layers2_label.grid(row=5, column=2, sticky="ew")
calc_layers_to_entry = Entry(dic_settings_frame, width=5)
calc_layers_to_entry.grid(row=5, column=3, sticky="ew")

threshold_label = Label(dic_settings_frame, text="coef threshold:")
threshold_label.grid(row=6, column=0, sticky="ew")
threshold_label_entry = Entry(dic_settings_frame, width=5)
threshold_label_entry.insert(END, '0.15')
threshold_label_entry.grid(row=6, column=1, sticky="ew")

delta_threshold_label = Label(dic_settings_frame, text="delta coef threshold:")
delta_threshold_label.grid(row=7, column=0, sticky="ew")
delta_threshold_label_entry = Entry(dic_settings_frame, width=5)
delta_threshold_label_entry.insert(END, '0.015')
delta_threshold_label_entry.grid(row=7, column=1, sticky="ew")

crack_gradient_label = Label(dic_settings_frame, text="crack gradient threshold:")
crack_gradient_label.grid(row=8, column=0, sticky="ew")
crack_gradient_label_entry = Entry(dic_settings_frame, width=5)
crack_gradient_label_entry.insert(END, '0.01')
crack_gradient_label_entry.grid(row=8, column=1, sticky="ew")

# RUN DIC
run_label_button = Button(dic_settings_frame, text="Run DIC", command=run_dic)  # Define browse_pt_image function later
run_label_button.grid(row=9, column=0, padx=10, sticky="ew")


backward_value_label = Label(buttons_frame, text="Calculation:")
backward_value_label.grid(row=1, column=0, sticky="ew")
backward_value_button = Button(buttons_frame, text="Forward", command=change_backward)
backward_value_button.grid(row=1, column=1, sticky="ew")

# PT
scale_size_label = Label(pt_settings_name_frame, text="Scale ")
scale_size_label.grid(row=0, column=0, sticky="ew")
scale_x_size_label = Label(pt_settings_name_frame, text="x:")
scale_x_size_label.grid(row=0, column=1, sticky="ew")
scale_x_size_label_entry = Entry(pt_settings_name_frame, width=5)
scale_x_size_label_entry.grid(row=0, column=2, sticky="ew")
scale_y_size_label = Label(pt_settings_name_frame, text="y:")
scale_y_size_label.grid(row=0, column=3, sticky="ew")
scale_y_size_label_entry = Entry(pt_settings_name_frame, width=5)
scale_y_size_label_entry.grid(row=0, column=4, sticky="ew")
scale_z_size_label = Label(pt_settings_name_frame, text="z:")
scale_z_size_label.grid(row=0, column=5, sticky="ew")
scale_z_size_label_entry = Entry(pt_settings_name_frame, width=5)
scale_z_size_label_entry.grid(row=0, column=6, sticky="ew")
scale_size_label_end_label = Label(pt_settings_name_frame, text="um/px")
scale_size_label_end_label.grid(row=0, column=7, sticky="ew")

pt_search_range_label = Label(pt_settings_name_frame, text="PT search range:")
pt_search_range_label.grid(row=1, column=0, sticky="ew")
pt_search_range_label_entry = Entry(pt_settings_name_frame, width=5)
pt_search_range_label_entry.grid(row=1, column=2, sticky="ew")

run_label_button = Button(pt_settings_name_frame, text="Run PT", command=run_pt)  # Define browse_pt_image function later
run_label_button.grid(row=2, column=0, sticky="ew")



case_value_label = Label(case_frame, text="Case:")
case_value_label.grid(row=0, column=0, sticky="ew")
case_value_button = Button(case_frame, text="3D", command=change_case)
case_value_button.grid(row=0, column=1, sticky="ew")

# SHOW
show_text_label = Label(show_frame, text="Show displacements")
show_text_label.grid(row=0, column=0, sticky="ew")
show_dic_or_pt_button = Button(show_frame, text="DIC", command=show_dic_or_pt)
show_dic_or_pt_button.grid(row=0, column=1, sticky="ew")

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

show_button = Button(show_frame, text="Show", command=show_disps_dic_or_pt)
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



# Run the main loop
root.mainloop()
