import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D, axes3d
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter


import os
import inspect
import sys
if sys.version_info[0] < 3:
    import Tkinter as tk
    from Tkinter import filedialog
    from Tkinter import ttk

else:
    import tkinter as tk
    from tkinter import filedialog
    from tkinter import ttk

import numpy as np
import trimesh as trimesh
import pandas as pd
from scipy import interpolate as interp
import functools


class vectorize(np.vectorize):
    def __get__(self, obj, objtype):
        return functools.partial(self.__call__, obj)

class Composites():
    def __init__(self):
        return
    
    def UD_formulas(self):   
        def E_UD1(self, E_f1, phi, E_m):
            return 0.975 * (E_f1 * phi + E_m * ( 1 - phi) )
        
        def E_UD2(self, E_m, E_f2, phi, zeta=1):
            nu_E = ((E_f2/E_m) - 1 )/((E_f2/E_m) + zeta)
            return E_m * (1 + (zeta*nu_E*phi))/(1-(nu_E*phi))                              
    
    def Laminate_formulas(self):
        def Youngs_modulus(weave, E_UD1, E_UD2, G_BD=None, poisson_BD=None):
            if weave == 'CSM':
                return (3/8)*E_UD1 + (5/8)*E_UD2
            elif weave == 'BD' or weave == '0/90' or weave == '+':
                return 0.5*(E_UD1 + E_UD2)
            elif weave == 'DB' or weave == '+-45' or weave == 'x':
                if G_BD is None or poisson_BD is None:
                    raise ValueError("Need to define both G_BD and poisson_BD for this calculation")
                E_BD = 0.5*(E_UD1 + E_UD2)
                return ( 4 * E_BD ) / ((E_BD / G_BD) + 2*(1-poisson_BD))

class Hull(Composites):
    def __init__(self):
        global root 
        root = tk.Tk()
         
        
        self.properties = {}
        self.chine_beam_label = tk.StringVar()
        self.chine_beam_label.set("Chine beam   :   Undefined")
        self.materials_label = tk.StringVar()
        self.materials_label.set("Undefined")
        
        
        def on_closing():
            if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
                root.quit()
                root.destroy()
        root.protocol("WM_DELETE_WINDOW", on_closing)
        
        tk.Button(root, text='Load', command = lambda : self.load()).grid(row=1, column = 1, sticky='N')

        
        tk.Label(root, text='Material selection')\
            .grid(row = 1, column = 2, sticky='N')
           
        tk.Button(root, text='Show Hull',
                                command = lambda : self.mesh.show())\
            .grid(row=2, column = 1, sticky='N')
        
        self.material_input_frame(row = 2, column = 2)
        
        def design_pressures_button():
            try:
                self.ISO_12215(b=self.b, l=self.l, l_u=self.l_u, s=self.s, walking_area=self.walking_area)
            except AttributeError:
                self.ISO_12215()
                
        
        tk.Button(root, text='Calculate Design Pressures',
                            command = design_pressures_button )\
                                .grid(row = 1, column = 4)
        
        self.panel_input_frame(row = 2, column = 4)
        
        def get_class_properties():
            class_properties = []
            for k, d in self.__dict__.items():
                try:
                    float(d)
                except:
                    continue
                if k.startswith('_') is False:
                    class_properties.append((k, float(d)))
                else:
                    continue
            return '\n'.join('{}'.format(k) for k, d in class_properties), '\n'.join('{:.6g}'.format(d) for k, d in class_properties)
        
        properties_frame = tk.LabelFrame(root, text='Vessel Data')
        properties_frame.grid(row = 1, column = 5, rowspan=8)   
        properties_key = tk.StringVar()
        properties_value = tk.StringVar()
        properties_key.set(get_class_properties()[0])
        properties_value.set(get_class_properties()[1])
        tk.Label(properties_frame, textvariable=properties_key, anchor = 'e', justify='left')\
            .grid(row = 1, column = 1)
        tk.Label(properties_frame, textvariable=properties_value, anchor = 'e', justify='right')\
            .grid(row = 1, column = 2)
        
        def update():
            for attribute, value in self.__dict__.items():
                self.properties[attribute] = value
                properties_key.set(get_class_properties()[0])
                properties_value.set(get_class_properties()[1])
            root.update_idletasks()
            root.after(1000, update)
            
        update()
        
        self._deck_kl_figures = tk.Frame(root)
        self._deck_kl_figures.grid(row = 8, column = 1, columnspan=4)
        
        self._design_pressures_frame1 = tk.Frame(root)
        self._design_pressures_frame1.grid(row=9, column = 1, columnspan=6)
        
        root.mainloop()
    
        
    
    def load(self):
        try:
            del(self.B_c)
        except AttributeError:
            pass
        file_dialog = tk.Tk()
        file_dialog.withdraw()
        self.filepath = filedialog.askdirectory(title='Choose the folder containing the desired hull .stl and Maxsurf data')
        file_dialog.destroy()
        number_stls = 0
        try:
            for file in os.listdir(self.filepath):
                if file.endswith('stl'):
                    number_stls += 1
                    self.hull_path = os.path.join(self.filepath, file)
        except FileNotFoundError:
            print("You must select a file.")
            exit()
        
        if number_stls > 1:
            self.hull_path = filedialog.askopenfilename(title ='More than one .stl in folder, please specify')
            print(self.hull_path)

        self.hull_stl = os.path.abspath(self.hull_path)
        self.mesh = trimesh.load(self.hull_stl)
        # self.mesh.show()
        
        self.design_cat = 'A'
        self.speed = 25.0 # knots. Maximum previewed speed with apparent wind between 60 and 90 degrees in m_LDC condition, with B_c = b_WL

        self.LOA = max(self.mesh.vertices[:, 0]) - min(self.mesh.vertices[:, 0])
        self.aft_location = min(self.mesh.vertices[:, 0])

        self.equillibrium = pd.read_csv(os.path.abspath(os.path.join(self.filepath, 'Equilibrium.txt')), sep='\t', lineterminator='\n', index_col=0, skiprows=0).transpose().dropna()
        self.LWL = self.equillibrium['WL Length m'].to_numpy()[0]
        self.displacement = self.equillibrium['Displacement t'].to_numpy()[0] * 1000
        self.B_waterline = self.equillibrium['Beam max extents on WL m'].to_numpy()[0]
        self.m_LDC = self.displacement + (2 * 90)
        self.find_draught()
        
        self.stability = pd.read_csv(os.path.abspath(os.path.join(self.filepath, 'Stability.txt')), sep='\t', lineterminator='\n', index_col=0, skiprows=0).transpose().dropna()
        self.max_GZ = self.max_GZ_under_60_deg_heel()

        self.z_Qs = np.linspace(0, self.canoe_body_draught - min(self.mesh.vertices[:, 2]), num=25)
        self.Z_SDT = lambda x : ( (0.0286 * self.LWL + 0.0115) * (x / self.LWL) ) + (0.0571 * self.LWL) + 0.229
        
        self.find_hull_depth()
    
    def material_input_frame(self, row, column):
        self.B1_materials = pd.read_excel(os.path.join(os.getcwd(), "ISO Materials.xlsx"), sheet_name="Table B1", index_col=0)
        df = self.B1_materials.T
        categories = [cat.replace('\'', '').replace(']', '').replace('[', '') for cat in np.unique(self.B1_materials.loc['Category'])]
        
        self.material_frame = tk.Frame(root)
        
        material_category_box= ttk.Combobox(self.material_frame, width = 37, value=(categories))
        material_category_box.grid(row=3, column=1, columnspan=2, padx=10, pady=2, sticky='w')
        
        # filter materials to only those matching the chosen category
        def category_callback(eventObject):
            _ = eventObject.widget.get()
            self.material_category = material_category_box.get()
            self.filtered_materials = df.loc[df['Category'] == self.material_category]
            types = [t.replace('\'', '').replace(']', '').replace('[', '') for t in np.unique(self.filtered_materials['Type'])]
            material_type_box.config(values=types)
        
        # Make a dropdown boc to choose Type (Steel, Al, Fibre, Matrix, Core)
        material_type_box = ttk.Combobox(self.material_frame, width=37)
        material_type_box.grid(row=4, column=1, columnspan=2, padx=10, pady=2, sticky='w')
        material_type_box.bind('<Button-1>', category_callback)  
        
        # filter materials to only those matching the chosen type
        def type_callback(eventObject):
            _ = eventObject.widget.get()
            self.material_type = material_type_box.get()
            self.filtered_materials = self.filtered_materials.loc[self.filtered_materials['Type'] == self.material_type]
            materials = [t.replace('\'', '').replace(']', '').replace('[', '') for t in np.unique(self.filtered_materials.T.columns)]
            material_box.config(values=materials)
        
        # Make a dropdown box to choose Material 
        material_box = ttk.Combobox(self.material_frame, width=37)
        material_box.grid(row=5, column=1, columnspan=2, padx=10, pady=2, sticky='w')
        material_box.bind('<Button-1>', type_callback)

        selected_mat = tk.StringVar()
        # Update the data box
        def update_data(eventObject=None):
            self.material = self.B1_materials[material_box.get()].dropna()
            selected_mat.set(self.material)
        
        material_box.bind('<<ComboboxSelected>>', update_data)
        tk.Label(self.material_frame, textvariable=selected_mat).grid(row=6, column=1, sticky='w')
       
        # define function that checks if the selected material is a composite, and if so calls the self.SandwichInput() function
        def check_for_sandwich():
            if self.material_category == 'FRP':
                
                self.sandwich_input_window()
            
            return
        
        # Button for final input
        tk.Button(self.material_frame, text="Update", width=37, height=2, command=check_for_sandwich).grid(row=7, column=1, sticky='w')
        self.material_frame.grid(row=row, column=column, columnspan=2, rowspan=4)
        
    def sandwich_input_window(self):
        sandwich_window = tk.Toplevel(root)
        
        self.SInput = tk.Frame(sandwich_window)
        self._num_top_rows = 2
        number_table_rows = 1
        self._rows = number_table_rows + self._num_top_rows
                
        tk.Button(sandwich_window, text="Add row", command=lambda : self._add_row())\
            .grid(row = 3, column = 2)
        
        tk.Button(sandwich_window, text="Submit", command=lambda : [self.get(), root.update(), sandwich_window.destroy()])\
            .grid(row = 4, column = 2)
        
        def fibre_volume_density_reference():
            reference_window = tk.Toplevel(sandwich_window)
            df = pd.read_excel(os.path.join(os.getcwd(), "ISO Materials.xlsx"), sheet_name="Table C7", index_col=0)
            cols = list(df.columns)
            tree = ttk.Treeview(reference_window)
            tree.pack()
            tree["columns"] = cols
            for i in cols:
                tree.column(i, anchor="w")
                tree.heading(i, text=i, anchor='w')

            for index, row in df.iterrows():
                tree.insert("",0,text=index,values=list(row))
            reference_window.mainloop()

        # re-fetch the table of materials
        df = self.B1_materials.T
        
        # List of available materials for the Sandwich (fibre or core materials)
        self.filtered_materials = df.loc[df['Type'] == 'Fibre'].T.columns.append(df.loc[df['Type'] == 'Core'].T.columns)

        # List of materials used in the sandwich
        self.material_list = []
        
        # List of available fibre weaves / orientations
        self.weaves = ['CSM', 'UD', 'BD+', 'DBx', 'Qx']
        
        # List of weaves used in the sandwich
        self.weave_list = []

        # List of thicknesses used in the sandwich
        self.thickness_list = []
        # List of fibre volume densities used in the sandwich
        self.fibre_volume_density_list = []
        
        # List used to store the table of OptionMenues
        self.option_menus = []
        
        tk.Label(self.SInput, text='Lamination Process:').grid(row = 1, column = 1, sticky='N', padx = 30, pady= 7)
        
        # Selection of layup process
        self.lamination_process = tk.StringVar()
        processes = ['Hand layup, simple surface', 'Hand layup, complex surface', 'RTM ECO', 'Infusion', 'Prepreg void', 'Prepreg autoclave']
        tk.OptionMenu(self.SInput, self.lamination_process, *processes).grid(row = 1, column = 2, sticky='N', padx = 30, pady= 7)
        
        tk.Button(self.SInput, text='Reference Fibre volume densities', command=fibre_volume_density_reference).grid(row = 1, column = 4, sticky='E', padx = 30, pady= 7)
        
        tk.Label(self.SInput, text='Material', font='Helvetica 18 bold').grid(row = 2, column = 1, sticky='N', padx = 30, pady= 7)
        tk.Label(self.SInput, text='Weave / Orientation', font='Helvetica 18 bold').grid(row = 2, column = 2, sticky='N', padx = 30, pady= 7)
        tk.Label(self.SInput, text='Thickness (mm)', font='Helvetica 18 bold').grid(row = 2, column = 3, sticky='N', padx = 30, pady= 7)
        tk.Label(self.SInput, text='Fibre Volume Density', font='Helvetica 18 bold').grid(row = 2, column = 4, sticky='N', padx = 30, pady= 7)
        
        for row in range(self._num_top_rows, self._rows, 1):
            self._create_row(row)

        self.SInput.grid_rowconfigure(self._rows, weight=1)
        self.SInput.grid(row = 2, column = 1, columnspan = 4)
        sandwich_window.mainloop()

    def _create_row(self, row):
            table_row = row - self._num_top_rows
            
            vcmd = (self.SInput.register(self._validate),
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
            self.option_menus.append([])
            
            material_selected = tk.StringVar()
            self.option_menus[table_row].append(tk.OptionMenu(self.SInput, material_selected, *self.filtered_materials)\
                .grid(row = row + 1, column = 1, sticky='N', padx = 30, pady= 7))
            self.material_list.append(material_selected)
            
            weave_selected = tk.StringVar()
            self.option_menus[table_row].append(tk.OptionMenu(self.SInput, weave_selected, *self.weaves)\
                .grid(row = row + 1, column = 2, sticky='N', padx = 30, pady= 7))
            self.weave_list.append(weave_selected)
            
            thickness_selected = tk.StringVar()
            self.option_menus[table_row].append(tk.Entry(self.SInput, textvariable = thickness_selected, validate='key', validatecommand=vcmd)\
                .grid(row = row + 1, column = 3, sticky='N', padx = 30, pady= 7))
            self.thickness_list.append(thickness_selected)
            
            fibre_volume_density_selected = tk.StringVar()
            self.option_menus[table_row].append(tk.Entry(self.SInput, textvariable = fibre_volume_density_selected, validate='key', validatecommand=vcmd)\
                .grid(row = row + 1, column = 4, sticky='N', padx = 30, pady= 7))
            self.fibre_volume_density_list.append(fibre_volume_density_selected)
            
    def get(self):
        '''Return a list of lists, containing the data in the table'''
        self.material = pd.DataFrame(np.asarray([[i.get() for i in self.material_list],
                                        [i.get() for i in self.weave_list],
                                        [i.get() for i in self.thickness_list],
                                        [i.get() for i in self.fibre_volume_density_list]]).transpose(),
                            columns = ['Material', 'Weave', 'Thickness', 'Fibre Volume Density'])
        
        tk.Label(root, text='Sandwich information').grid(row = 6, column = 2)
        self.material_tree = ttk.Treeview(root)
        self.material_tree.grid(row = 7, column = 2)
        self.material_tree['columns'] = self.material.columns.tolist()
        
        for col in self.material.columns.tolist():
            self.material_tree.column(col, anchor='w')
            self.material_tree.heading(col, text=col, anchor='w')
        
        for index, row in self.material.iterrows():
            self.material_tree.insert("", "end", text=index, values=list(row))
        
    def _add_row(self):
        row = self._rows
        self._create_row(row)
        self._rows += 1

    def _validate(self, action, index, value_if_allowed,
                    prior_value, text, validation_type, trigger_type, widget_name):
        if value_if_allowed:
            try:
                float(value_if_allowed)
                return True
            except ValueError:
                return False
        else:
            return False
    
    def panel_input_frame(self, row, column):
        panel_input_frame = tk.LabelFrame(root, text="Structural Dimensions")
        panel_frame = tk.LabelFrame(panel_input_frame, text="Panels")
        stiffener_frame = tk.LabelFrame(panel_input_frame, text="Stiffeners")
        b = tk.StringVar()
        l = tk.StringVar()
        l_u = tk.StringVar()
        s = tk.StringVar()
        walking_area = tk.StringVar()
        walking_area.set("False")
        
        def walking_button_toggle():
            if str(walking_area.get()) == "False":
                walking_area.set("True")
                btn.config(text="True")
            elif str(walking_area.get()) == "True":
                walking_area.set("False")
                btn.config(text="False")
            else:
                raise(ValueError("walking area is not \"True\" or \"False\""))
        
        def update():
            self.b = float(b.get())
            self.l = float(l.get())
            self.l_u = float(l_u.get())
            self.s = float(s.get())
            self.walking_area = bool(walking_area.get())
        
        vcmd = (panel_frame.register(self._validate),
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        
        tk.Label(panel_frame, text="Short unsupported dimension (mm):").grid(row = 0, column = 0)
        tk.Entry(panel_frame, textvariable = b, validate = 'key', validatecommand=vcmd).grid(row = 1, column = 0)
        tk.Label(panel_frame, text="Long unsupported dimension (mm):").grid(row = 2, column = 0)
        tk.Entry(panel_frame, textvariable = l, validate = 'key', validatecommand=vcmd).grid(row = 3, column = 0)
        tk.Label(panel_frame, text="Walking Area:").grid(row = 4, column = 0)
        btn = tk.Button(panel_frame, text="False", command = lambda : walking_button_toggle())
        btn.grid(row = 5, column = 0)
        
        tk.Label(stiffener_frame, text="Small dimension (spacing) of a stiffener between axis (mm):").grid(row = 0, column = 0)
        tk.Entry(stiffener_frame, textvariable = s, validate = 'key', validatecommand=vcmd).grid(row = 1, column = 0)
        tk.Label(stiffener_frame, text="Large dimension (span) of a stiffener between axis (mm):").grid(row = 2, column = 0)
        tk.Entry(stiffener_frame, textvariable = l_u, validate = 'key', validatecommand=vcmd).grid(row = 3, column = 0)
        
        panel_frame.grid(row = 0, column = 0)
        stiffener_frame.grid(row = 1, column = 0)
        tk.Button(panel_input_frame, text='Update Structural sizes', command=lambda : update()).grid(row =2, column= 0)
        
        panel_input_frame.grid(row = row, column = column)
    
    def find_hull_depth(self):
        zs = self.mesh.vertices[:, 2]
        self.depth = max(zs) - min(zs)
        return self.depth

    def plot_hull(self):
        xs = self.mesh.vertices[:, 0]
        ys = self.mesh.vertices[:, 1]
        zs = self.mesh.vertices[:, 2]
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(xs, ys, zs)
        return fig, ax

    def find_draught(self):
        T_midship = self.equillibrium['Draft Amidships m']
        T_FP = self.equillibrium['Draft at FP m']
        T_AP = self.equillibrium['Draft at AP m']
        T_LCF = self.equillibrium['Draft at LCF m']
        draughts = np.array([T_midship, T_AP, T_FP, T_LCF])
        self.canoe_body_draught = max(draughts)
        return self.canoe_body_draught
    
    def calc_z_SDT(self):
        self.z_SDT = lambda x: (0.0286 * self.LWL + 0.115) * (x / self.LWL) + (0.0571 * self.LWL)+ 0.229
    
        frame = tk.LabelFrame(self._deck_kl_figures, text="Theoretical Deck Height")
        frame.grid(row = 1, column = 1)
        f = Figure(figsize=(5,5), dpi=100)
        ax = f.add_subplot(111)
        xs = np.arange(0, self.LWL, 0.1)
        ax.plot(xs, self.z_SDT(xs))
        ax.set_xlabel('Distance from aft (m)')
        ax.set_ylabel(r'$z_{SDT}$')
        ax.grid()
        
        canvas = FigureCanvasTkAgg(f, frame)
        canvas.draw()

        toolbar = NavigationToolbar2Tk(canvas, frame)
        toolbar.update()
        canvas.get_tk_widget().pack()
        
        return self.z_SDT

    def max_GZ_under_60_deg_heel(self):
        indices = self.stability.index.to_numpy()
        max_gz = -np.infty
        for index in indices:
            try:
                if float(index) <= 60.0:
                    data = self.stability
                    gz = data.at[index, 'GZ m']
                    if max_gz < gz: max_gz = gz
            except ValueError:
                pass
        self.max_GZ = max_gz
        return self.max_GZ

    def calc_k_SLS(self):
        self.max_GZ_under_60_deg_heel()
        if self.m_LDC > 5 * self.LWL**3:
            self.k_SLS = 1
        else:
            a = 10 * self.max_GZ * (self.LWL**0.5)
            b = self.m_LDC**0.33
            sls = (a/b)**0.5
            if sls < 1:
                self.k_SLS = 1
            else:
                self.k_SLS = sls
            
        return self.k_SLS

    def get_section(self, origin, plane=[1, 0, 0], plot=False):
        intersections = trimesh.intersections.mesh_plane(self.mesh, plane, origin)
        fig = plt.figure()
        for line in intersections:
            start, end = line
            _, y_start, z_start = start
            _, y_end, z_end = end

            if plot:
                fig.plot([y_start, y_end], [z_start, z_end])
        if plot:
            fig.xlabel('y : beam (m)')
            fig.ylabel('z : draught (m)')
            fig.show()

        return intersections

    def calc_BC(self, order=3, s=0.01):
        
        x = self.aft_location + 0.4 * self.LWL
        grad_50 = -np.tan(np.radians(50))  # calculate a 50 deg angle as a gradient

        section = self.get_section([x, 0, 0])
        section = np.reshape(section, (-1, 3)) # flatten the array so it's no longer 3D
        section = section[section[:, 1].argsort()]  # sort the coordinated by y value
        section = np.unique([tuple(row) for row in section], axis=0) # remove duplicate points
        
        y = section[:, 1]
        z = section[:, 2]
        tck = interp.splrep(y, z, s=s, k=order)  # pass a cubic spline through the section points
        y_pos = interp.splev(grad_50, tck, der=1)
        z_pos = interp.splev(y_pos, tck)

        z2 = interp.splev(y, tck)  # line for plotting
        
        plt.ion()
        fig, ax1 = plt.subplots()
        fig.canvas.set_window_title("Hull section at 0.4Lwl ({:.2f} m) fwd, gradient and location of B_c".format(x))
        ax1.set_title('Please click on the location of the chine beam. An estimate is shown in green. Location should be each at the chine, or where the hull gradient is 50 degrees.')
        ax1.set_xlabel('y : Beam (m)')
        ax1.grid()
        ax1.set_ylabel('Surface angle to horizontal (deg)')
        ax1.plot(y, -180*np.arctan(interp.splev(y, tck, der=1))/np.pi, label='Gradient')

        ax2 = ax1.twinx()
        ax2.set_aspect('equal')
        ax2.set_ylabel('z : Draught (m)')
        ax2.grid()
        ax2.plot(y, -1*z2, '-', label='Fitted spline')
        ax2.plot(y, -1*z, 'x-', label='Hull data points' )
        ax2.plot(y_pos, -1*z_pos, 'o', label='Estimated chine beam position')
        fig.tight_layout()
        # cid = fig.canvas.mpl_connect('button_press_event', onclick)
        fig.legend()
        coords = plt.ginput(timeout = 0, n=1)
        plt.close()

        self.B_c = 2 * abs(coords[-1][0])
        self.B_c_pos = (coords[-1][0], (coords[-1][1] - min(z)))
        return self.B_c, self.B_c_pos
    
    def calc_Beta_0d4(self):
        y, z = self.B_c_pos
        deadrise = np.tan((z / y))
        self.Beta_0d4 = np.degrees(deadrise)
        return self.Beta_0d4

    def calc_k_DYN1(self):
        if self.Beta_0d4<10.0:
            beta = 10.0
        elif self.Beta_0d4 > 30.0:
            beta = 30.0
        else:
            beta = self.Beta_0d4
        
        a = (self.LWL / (10 * self.B_c) ) + 0.084
        b = 50 - beta
        c = ( (self.speed ** 2) * (self.B_c ** 2) ) / self.m_LDC
        self.k_DYN1 = 0.32 * a * b * c
        return self.k_DYN1

    def calc_k_DYN2(self):
        if self.speed > 50.0:
            raise Warning("For For recreational and charter craft, the maximum speed shall not be taken >50 knots, but for “Heavy duty” work boats this speed may be greater (see Annex J).")
        k_DYN2 = (0.5 * self.speed) / (self.m_LDC ** 0.17)
        if k_DYN2 > 6: k_DYN2 = 6
        if k_DYN2 < 3: k_DYN2 = 3
        self.k_DYN2 = k_DYN2
        return self.k_DYN2

    def calc_k_DYN(self):
        try:
            self.B_c += 0
        except AttributeError:
            self.calc_BC()
        self.calc_Beta_0d4()
        DYN1 = self.calc_k_DYN1()
        DYN2 = self.calc_k_DYN2()
        if DYN1 < DYN2:
            self.k_DYN = DYN1
        else:
            self.k_DYN = DYN2
        
        self.k_DYN = 3.0
        return self.k_DYN
    
    def calc_k_DC(self):
        if self.design_cat == 'A':
            self.k_DC = 1.0
        else:
            raise Warning("Design category not A, this program is not set up for this")
        return self.k_DC

    def calc_k_L(self):
        """x in (m) is the longitudinal position of the centre of the panel or middle of stiffener forward of aft end of LWL in mLDC,
        x/LWL = 0 and 1 are respectively the aft end and fore end of LWL.

        Raises:
            Exception: k_SLS doesn't fall within the expected range

        Returns:
            lambda: a function to calculate k_L at a given position x along the hull
        """
        if self.k_SLS > 1.0:
            k_DYN = max([self.k_DYN1, self.k_DYN2])
            if self.B_waterline * k_DYN < 3: 
                print("B_c = B_WL * k_DYN shall not be taken <3")
        elif self.k_SLS == 1.0:
            k_DYN = 3
        else:
            raise Exception("self.k_SLS is less than 1")

        self.k_L = lambda x: ( (1.667 - ( 0.222 * k_DYN ) ) * (x / self.LWL) ) + (0.133 * k_DYN)\
                            if ( (1.667 - 0.222 * k_DYN) * (x / self.LWL) ) + (0.133 * k_DYN) < 1\
                            else 1

        xs = np.linspace(0, self.LWL, num=100)
        k_Ls = [self.k_L(x) for x in xs]
        
        frame = tk.LabelFrame(self._deck_kl_figures, text="Longitudinal Pressure Coefficient")
        frame.grid(row = 1, column = 2)
        f = Figure(figsize=(5,5), dpi=100)
        ax = f.add_subplot(111)
        # fig.canvas.set_window_title("Longitudinal pressure distribution factor")
        ax.set_xlabel(r"Longitudinal distance from aft $\frac{x}{LWL}$")
        ax.set_ylabel(r"$k_L$")
        ax.plot(xs / self.LWL, k_Ls, label='k_DYN = {:.4f}'.format(k_DYN))
        ax.set_ylim(0, 1)
        ax.set_xticks(np.arange(0, 1.01, 0.1))
        ax.set_yticks(np.arange(0, 1.01, 0.1))
        ax.set_xlim(0, 1)
        f.legend()
        ax.grid()
        # fig.show()
        
        canvas = FigureCanvasTkAgg(f, frame)
        canvas.draw()

        toolbar = NavigationToolbar2Tk(canvas, frame)
        toolbar.update()
        canvas.get_tk_widget().pack()

        return self.k_L

    def calc_k_AR(self, b, l, l_u, s):
        try:
            b = float(b)
            l = float(l)
            l_u = float(l_u)
            s = float(s)
        except:
            raise(TypeError("Inputs could not be converted to floats"))
        A_D_PLT = b * l * 10 ** (-6)
        A_D_STF = l_u * s * 10 ** (-6) if l_u * s * 10 ** (-6) > 0.33 * l_u ** 2 * 10 ** (-6) else 0.33 * l_u ** 2 * 10 ** (-6)

        k_R_PLT = 1.5 - (3 * 10 ** (-4) * b)
        k_R_STF = 1 - (2 * 10 ** (-4) * l_u)

        self.k_AR_PLT = (k_R_PLT * 0.1 * self.m_LDC ** 0.15) / (A_D_PLT ** 0.3)
        self.k_AR_STF = (k_R_STF * 0.1 * self.m_LDC ** 0.15) / (A_D_STF ** 0.3)
        self.k_AR = self.k_AR_PLT, self.k_AR_STF
        return self.k_AR
        
    def calc_P_BS(self):

        self.P_BS_BASE = (2 * (self.m_LDC ** 0.33) + 18) * self.k_SLS

        self.P_BS_MIN_PLT = max([
                            (0.3 * self.m_LDC ** 0.33 + 0.66 * self.LWL * self.k_DC),
                            10 * self.canoe_body_draught,
                            7
                            ])
        self.P_BS_MIN_STF = max(
                            [0.85 * self.P_BS_MIN_PLT,
                            7
                            ])

        self.P_BS_PLT = lambda x: (self.P_BS_BASE * self.k_AR_PLT * self.k_DC * self.k_L(x)) if (self.P_BS_BASE * self.k_AR_PLT * self.k_DC * self.k_L(x)) > self.P_BS_MIN_PLT else self.P_BS_MIN_PLT
        self.P_BS_STF = lambda x: (self.P_BS_BASE * self.k_AR_STF * self.k_DC * self.k_L(x)) if (self.P_BS_BASE * self.k_AR_STF * self.k_DC * self.k_L(x)) > self.P_BS_MIN_STF else self.P_BS_MIN_STF
        
        return self.P_BS_PLT, self.P_BS_STF

    @vectorize
    def calc_P_DS(self, x, z_Qs, walking_area=False):
        """Sailing craft deck and cockpit bottom design pressure

        Args:
            x (float): Longitudinal distance from the stern (m)
            z_Qs (float): Local height of a point Q, centre of a panel or stiffener above WL

        Returns:
            float, float : Design pressure for plate, stiffener
        """
        try:
            walking_area = bool(walking_area)
        except:
            raise(TypeError("Walking area could not be converted to Boolean"))
        P_MIN = 5.0 if walking_area else 3.5
        self.P_DS_BASE = (0.5 * (self.m_LDC ** 0.33)) + 12
        # self.P_DS = lambda x, z_Qs: (self.P_BS_BASE - (self.P_BS_BASE - P_DS_BASE) * min([(z_Qs / self.z_SDT), 1])) * self.k_AR * self.k_DC * self.k_L(x)
        P_DS_PLT = (self.P_BS_BASE - (self.P_BS_BASE - self.P_DS_BASE) * min([(z_Qs / self.z_SDT(x)), 1])) * self.k_AR_PLT * self.k_DC * self.k_L(x)
        P_DS_STF = (self.P_BS_BASE - (self.P_BS_BASE - self.P_DS_BASE) * min([(z_Qs / self.z_SDT(x)), 1])) * self.k_AR_STF * self.k_DC * self.k_L(x)
        if P_DS_PLT < P_MIN:
            P_DS_PLT = P_MIN
        if P_DS_STF < P_MIN:
            P_DS_STF = P_MIN
        return P_DS_PLT, P_DS_STF

    @vectorize
    def calc_P_SS(self, x, z_Qs):
        P_SS_MIN_PLT = max([self.P_BS_MIN_PLT - (z_Qs / self.z_SDT(x)) * (self.P_BS_MIN_PLT - 5), 5])
        P_SS_MIN_STF = max([0.85*P_SS_MIN_PLT, 5])

        P_SS_PLT = (self.P_BS_BASE - (self.P_BS_BASE - self.P_DS_BASE) * min([(z_Qs / self.z_SDT(x)), 1])) * self.k_AR_PLT * self.k_DC * self.k_L(x)
        P_SS_STF = (self.P_BS_BASE - (self.P_BS_BASE - self.P_DS_BASE) * min([(z_Qs / self.z_SDT(x)), 1])) * self.k_AR_STF * self.k_DC * self.k_L(x)
        
        if P_SS_PLT < P_SS_MIN_PLT:
            P_SS_PLT = P_SS_MIN_PLT
        if P_SS_STF < P_SS_MIN_STF:
            P_SS_STF = P_SS_MIN_STF
        return P_SS_PLT, P_SS_STF

    def calc_Design_Pressures(self, walking_area):
        
        ### bottom area pressure ###
        frame1 = tk.LabelFrame(self._design_pressures_frame1, text="Bottom")
        frame1.grid(row=1, column = 1)
        f1 = Figure(figsize=(5,5), dpi=100)
        canvas1 = FigureCanvasTkAgg(f1, frame1)
        canvas1.draw()
        canvas1.get_tk_widget().pack()

        toolbar = NavigationToolbar2Tk(canvas1, frame1)
        toolbar.update()
        canvas1._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        ax1 = f1.add_subplot(111)
        bottom_plate, bottom_stiffener = self.calc_P_BS()
        xs = np.linspace(0, self.LWL, num=50)
        self.P_BS_PLTs = [bottom_plate(x) for x in xs]
        self.P_BS_STFs = [bottom_stiffener(x) for x in xs]
        ax1.set_xlabel("Longitudinal position from aft (x/LWL)")
        ax1.set_ylabel("Design Pressure (kN/m2) ")
        ax1.plot( xs / self.LWL, self.P_BS_PLTs, label='Bottom Plate')
        ax1.plot( xs / self.LWL, self.P_BS_STFs, label='Bottom Stiffener')
        ax1.set_xlim(0, 1)
        ax1.set_xticks(np.arange(0, 1.01, 0.1))
        f1.legend()
        ax1.grid()
        
        ### Deck and cockpit bottom design pressure ###
        frame2 = tk.LabelFrame(self._design_pressures_frame1, text="Deck and Cockpit")
        frame2.grid(row=1, column = 2)
        f2 = Figure(figsize=(5,5), dpi=100)
        canvas2 = FigureCanvasTkAgg(f2, frame2)
        canvas2.draw()
        canvas2.get_tk_widget().pack()

        toolbar2 = NavigationToolbar2Tk(canvas2, frame2)
        toolbar2.update()
        canvas2._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        ax2 = f2.add_subplot(111, projection='3d')
        xx, zz = np.meshgrid(xs, self.z_Qs)
        self.P_DS_PLT, self.P_DS_STF = self.calc_P_DS(x = xx, z_Qs = zz, walking_area=walking_area)
        surf_plate = ax2.plot_surface(xx / self.LWL, self.P_DS_PLT, zz, label=r"Deck and cockpit bottom pressure- Plate")
        surf_stiff = ax2.plot_surface(xx / self.LWL, self.P_DS_STF, zz, label=r"Deck and cockpit bottom pressure- Stiffener")
        surf_plate._facecolors2d=surf_plate._facecolors3d
        surf_stiff._facecolors2d=surf_stiff._facecolors3d
        ax2.set_xlabel("Longitudinal position from aft (x/LWL)")
        ax2.set_ylabel("Design Pressure (kN/m2)")
        ax2.set_zlabel("Height of z_Qs, panel or stiffener center, above the waterline (m)")
        # f2.legend()

        ### Side design pressure ###
        frame3 = tk.LabelFrame(self._design_pressures_frame1, text="Side")
        frame3.grid(row=1, column = 3)
        f3 = Figure(figsize=(5,5), dpi=100)
        canvas3 = FigureCanvasTkAgg(f3, frame3)
        canvas3.draw()
        canvas3.get_tk_widget().pack()

        toolbar3 = NavigationToolbar2Tk(canvas3, frame3)
        toolbar3.update()
        canvas3._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        ax3 = f3.add_subplot(111, projection='3d')
        xx, zz = np.meshgrid(xs, self.z_Qs)
        self.P_SS_PLT, self.P_SS_STF = self.calc_P_SS(x = xx, z_Qs = zz)
        surf_plate = ax3.plot_surface(xx / self.LWL, self.P_SS_PLT, zz, label=r"Side pressure - Plate")
        surf_stiff = ax3.plot_surface(xx / self.LWL, self.P_SS_STF, zz, label=r"Side pressure - Stiffener")
        surf_plate._facecolors2d=surf_plate._facecolors3d
        surf_stiff._facecolors2d=surf_stiff._facecolors3d
        ax3.set_xlabel("ax3")
        # ax3.set_xlabel("Longitudinal position from aft (x/LWL)")
        ax3.set_ylabel("Design Pressure (kN/m2)")
        ax3.set_zlabel("Height of z_Qs, panel or stiffener center, above the waterline (m)")
        f3.legend()

        return (xs, self.P_BS_PLTs, self.P_BS_STFs), (xx, zz, self.P_DS_PLT, self.P_DS_STF), (xx, zz, self.P_SS_PLT, self.P_SS_STF)

    def calc_k_BB(self):
        if self.build_quality is not None and self.build_method is not None:
            method_in_list = True
            if self.build_quality == 'TESTED':
                if self.build_method == 'Hand laid':
                    self.k_BB = 1
                elif self.build_method == 'Infusion, pre-preg or similar':
                    self.k_BB = 1
                else:
                    method_in_list = False
            elif self.build_quality == 'HIGH':
                if self.build_method == 'Hand laid':
                    self.k_BB = 0.95
                elif self.build_method == 'Infusion, pre-preg or similar':
                    self.k_BB = 1
                else:
                    method_in_list = False
            elif self.build_quality == 'LOW':
                if self.build_method == 'Hand laid':
                    self.k_BB = 0.75
                elif self.build_method == 'Infusion, pre-preg or similar':
                    self.k_BB = 0.80
                else:
                    method_in_list = False
            else:
                raise Exception("Build quality is not correct (none of ['TESTED', 'HIGH', 'LOW']).")
            if method_in_list == False:
                raise Exception("Build method is not correct (none of ['Hand laid', 'Infusion, pre-preg or similar']).")
        else:
            self.k_BB = 1
            raise Warning("No build quality or method defined. This is only applicable to non-FRP constructions")
        return self.k_BB
    
    def calc_k_AM(self):
        if self.assessment_method is not None and self.material_type is not None:
            material_in_list = True
            if self.assessment_method == 1:
                if self.material_type == 'FRP' or self.material_type == 'wood':
                    self.k_AM = 0.9
                elif self.material_type == 'metal':
                    self.k_AM = 1
                else:
                    material_in_list = False
            elif self.assessment_method == 2:
                if self.material_type == 'FRP' or self.material_type == 'wood':
                    self.k_AM = 0.95
                elif self.material_type == 'metal':
                    self.k_AM = 1
                else:
                    material_in_list = False
            elif self.assessment_method == 3:
                if self.material_type == 'FRP' or self.material_type == 'wood':
                    self.k_AM = 1
                elif self.material_type == 'metal':
                    self.k_AM = 1
                else:
                    material_in_list = False
            elif self.assessment_method == 4:
                if self.material_type == 'FRP' or self.material_type == 'wood':
                    self.k_AM = 1
                elif self.material_type == 'metal':
                    self.k_AM = 1
                else:
                    material_in_list = False
            else:
                raise Exception("This assessment method is incorrect (none of [1, 2, 3, 4]).")
            if material_in_list == False:
                    raise Exception("This material type is incorrect (none of ['FRP', 'wood', 'metal']).")
        else:
            raise Exception("You must define self.assessment_method and self.material_type")
        return self.k_AM

    def ISO_12215(self, b=None, l=None, l_u=None, s=None, walking_area=None):

        # 12215-5, 6.2 Height of the theoretical hull deck limit
        self.calc_z_SDT()

        # 12215-5, 8.7 Pressure correcting factor for slamming
        self.calc_k_SLS()
        
        # 12215-5, 8.3 Dynamic Load factor
        self.calc_k_DYN()

        # 12215-5, 8.2 Design Category factor
        self.calc_k_DC()

        # 12215-5, 8.4 Longitudinal pressure distribution
        self.calc_k_L()

        # 12215-5, 8.5 Area pressure reduction factor
        if np.all([b != None, l != None, l_u != None, s != None]):
            self.calc_k_AR(b, l, l_u, s)

            # 12215-5, 9.2 Design pressure for a sailing craft
            self.calc_Design_Pressures(walking_area)
        
        # 12215-5, 10.1 Boat building quality factor
        # self.calc_k_BB()

        # 12215-5, 10.2 Assessment method factor
        # self.calc_k_AM()

    def print_HullScant_data(self):
        print('Loaded Displacement: {}'.format(self.m_LDC))
        print('Maximum GZ under 60 deg: {:.4f}'.format(self.max_GZ))
        print('Length of hull: {}'.format(self.LOA))
        print('Waterline Length: {}'.format(self.LWL))
        print('Waterline beam: {}'.format(self.B_waterline))
        print('Canoe Body draft: {}'.format(self.canoe_body_draught))
        print('Hull Overall Depth to deck: {}'.format(self.find_hull_depth()))
        print('Baseline Position below the loaded waterline: {}'.format(0))
        print('Deadrise Angle: {:.4f}'.format(self.Beta_0d4))
        print('Max speed: {} '.format(self.speed))

    

if __name__ == "__main__":
    Hull()
