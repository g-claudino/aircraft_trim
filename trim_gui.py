import tkinter as tk
import matplotlib
matplotlib.use('TkAgg')

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)

import pandas as pd
import os
import numpy as np


class App:
    def __init__(self):
        path = os.getcwd()
        self.foil_w = pd.read_table(path + '/Foils/23012.dat',skiprows=1,names=["X", "Y"], sep='\s+')
        self.foil_t = pd.read_table(path + '/Foils/0012.dat',skiprows=1,names=["X", "Y"], sep='\s+')
        self.polar_w = pd.read_table(path + '/Polars/23012_re_1e6.dat',sep='\s+')
        self.polar_t = pd.read_table(path + '/Polars/0012_re_1e6.dat', sep='\s+')

    # The main tkinter window 
        self.root = tk.Tk() 

    # setting the title and 
        self.root.title('Trimmer') 

    # setting the dimensions of 
    # the main window 
        self.root.geometry("1200x1200") 

    # button that would displays the plot 
        self.plot_button = tk.Button(master = self.root, 
                        height = 2, 
                        width = 10, 
                        text = "Trim",
                        command= self.ac_trim) 
    # place the button 
    # into the window 
        self.plot_button.pack() 

    # Wing Sliders
        self.w_mac = tk.DoubleVar()
        self.chord_w = tk.Scale(master=self.root, from_=0.10, to=10.0, orient='horizontal', label='W MAC', variable=self.w_mac, digits = 2, resolution = 0.1, command=self.scaler_event)
        self.chord_w.set(1.5)
        self.chord_w.pack()
        self.w_aoa = tk.DoubleVar()
        self.aoa_w = tk.Scale(master=self.root, from_=-10.0, to=10.0, orient='horizontal', label='W AOA', variable=self.w_aoa, digits = 2, resolution = 0.1, command=self.scaler_event)
        self.aoa_w.set(2.5)
        self.aoa_w.pack()

        self.w_x = tk.DoubleVar()
        self.x_w = tk.Scale(master=self.root, from_=-10.0, to=10.0, orient='horizontal', label='W X', variable=self.w_x, digits = 2, resolution = 0.1, command=self.scaler_event)
        self.x_w.set(0.0)
        self.x_w.pack()

        self.w_y = tk.DoubleVar()
        self.y_w = tk.Scale(master=self.root, from_=-10.0, to=10.0, orient='horizontal', label='W Y', variable=self.w_y, digits = 2, resolution = 0.1, command=self.scaler_event)
        self.y_w.set(0.5)
        self.y_w.pack()

    # Tail Sliders
        self.t_mac = tk.DoubleVar()
        self.chord_t = tk.Scale(master=self.root, from_=0.10, to=5.0, orient='horizontal', label='T MAC', variable=self.t_mac, digits = 2, resolution = 0.1, command=self.scaler_event)
        self.chord_t.set(0.8)
        self.chord_t.pack()
        
        self.t_aoa = 0

        self.t_x = tk.DoubleVar()
        self.x_t = tk.Scale(master=self.root, from_=-10.0, to=10.0, orient='horizontal', label='T X', variable=self.t_x, digits = 2, resolution = 0.1, command=self.scaler_event)
        self.x_t.set(8.0)
        self.x_t.pack()

        self.t_y = tk.DoubleVar()
        self.y_t = tk.Scale(master=self.root, from_=-10.0, to=10.0, orient='horizontal', label='T Y', variable=self.t_y, digits = 2, resolution = 0.1, command=self.scaler_event)
        self.y_t.set(0.0)
        self.y_t.pack()

    # CG Sliders
        self.cg_x = tk.DoubleVar()
        self.x_cg = tk.Scale(master=self.root, from_=-5.0, to=5.0, orient='horizontal', label='CG X', variable=self.cg_x, digits = 2, resolution = 0.1, command=self.scaler_event)
        self.x_cg.set(2.5)
        self.x_cg.pack()

        self.cg_y = tk.DoubleVar()
        self.y_cg = tk.Scale(master=self.root, from_=-5.0, to=5.0, orient='horizontal', label='CG Y', variable=self.cg_y, digits = 2, resolution = 0.1, command=self.scaler_event)
        self.y_cg.set(0.0)
        self.y_cg.pack()



        # create a figure
        figure = Figure(figsize=(6, 4), dpi=100)

        # create FigureCanvasTkAgg object
        self.figure_canvas = FigureCanvasTkAgg(figure, self.root)

        # create the toolbar
        NavigationToolbar2Tk(self.figure_canvas, self.root)

        # create axes
        self.axes = figure.add_subplot()

        
        self.axes.plot(self.foil_w['X']*self.w_mac.get()+self.w_x.get(),self.foil_w['Y']*self.w_mac.get()+self.w_y.get())

        self.figure_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
        self.root.mainloop()

    def plot_update(self,event):
        x_wing = self.foil_w['X'].values*self.w_mac.get()
        y_wing = self.foil_w['Y'].values*self.w_mac.get()
        w_aoa = self.w_aoa.get()*np.pi/180
        self.x_wing_r = np.multiply(x_wing,np.cos(w_aoa))-np.multiply(y_wing,np.sin(w_aoa))+self.w_x.get()
        self.y_wing_r = np.multiply(x_wing,np.sin(w_aoa))+np.multiply(y_wing,np.cos(w_aoa))+self.w_y.get()
        
        x_tail = self.foil_t['X'].values*self.t_mac.get()
        y_tail = self.foil_t['Y'].values*self.t_mac.get()
        t_aoa = self.t_aoa*np.pi/180
        self.x_tail_r = np.multiply(x_tail,np.cos(t_aoa-w_aoa))-np.multiply(y_tail,np.sin(t_aoa-w_aoa))+self.t_x.get()
        self.y_tail_r = np.multiply(x_tail,np.sin(t_aoa-w_aoa))+np.multiply(y_tail,np.cos(t_aoa-w_aoa))+self.t_y.get()
        print(self.y_tail_r)
        
        self.axes.plot(self.x_wing_r,self.y_wing_r,'b'),
        self.axes.plot(self.x_tail_r,self.y_tail_r,'r'),
        self.axes.plot(self.cg_x.get(),self.cg_y.get(),'ko')
        self.axes.legend(['Wing','Tail','CG'])
        self.figure_canvas.draw()
        self.axes.clear()

    def ac_trim(self):
        w_CL = self.polar_w['CL'].loc[self.polar_w['alpha'] == self.w_aoa.get()].values
        w_CD = self.polar_w['CD'].loc[self.polar_w['alpha'] == self.w_aoa.get()].values
        w_CM = self.polar_w['CM'].loc[self.polar_w['alpha'] == self.w_aoa.get()].values
        w_x = 0.25*(np.min(self.x_wing_r)+np.max(self.x_wing_r))
        w_y = 0.25*(np.min(self.y_wing_r)+np.max(self.y_wing_r))
        w_part = w_CM+w_CL*(w_x-self.cg_x.get())+w_CD*(w_y-self.cg_y.get())
        w_mac = self.w_mac.get()

        t_CL = self.polar_t['CL'].values
        t_CD = self.polar_t['CD'].values
        t_x = 0.25*(np.min(self.x_tail_r)+np.max(self.x_tail_r))
        t_y = 0.25*(np.min(self.y_tail_r)+np.max(self.y_tail_r))
        t_part = t_CL*(t_x-self.cg_x.get())+t_CD*(t_y-self.cg_y.get())
        t_mac = self.t_mac.get()
        a_ht = np.interp(-w_part*w_mac,t_part*t_mac/w_mac,self.polar_t['alpha'].values)
        self.t_aoa = a_ht
        print(str(a_ht))
        self.plot_update(self)
        
    def scaler_event(self, event):
        self.plot_update(self)
    # run the gui 
app=App()

