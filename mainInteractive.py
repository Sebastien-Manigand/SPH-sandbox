#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 10:23:25 2020

@author: sebastien
"""

import tkinter as tk
from tkinter import ttk
import random as rd
import matplotlib.pyplot as plt
exec(open('GridSPH.py').read())

COLOR_BUTTON_BG_PLAY = '#88AA88'
COLOR_BUTTON_BG_STOP = '#AA8888'
MOVIE_DT = 25
DISPLAY_INFO_PARTICLE_NEIGHBORS = True

MODEL_OUTPUT = True
MODEL_OUTPUT_FILENAME = 'simulation_output.dat'


class AppliCanevas(tk.Tk):
    
    def __init__(self):
        self.running = False
#        self.infoPanel_opened = False

        self.stepNumber = 0
        tk.Tk.__init__(self)
        self.width = 840
        self.height = 480
        self.model = GridSPH()
#        self.model.testScene()
        self.model.uniformFlowScene()
        
        self.create_widgets()
        
        
        
        
        

    def create_widgets(self):
        # création canevas
        self.canv = tk.Canvas(self, bg="#EEEEFF", height=self.height, width=self.width)
        self.canv.pack(side=tk.LEFT)
        
        # boutons
        self.bouton_init = tk.Button(self, text="  Initialize  ", command=self.initialize, height=2, width=25)
        self.bouton_init.pack(side=tk.TOP)
        
        self.bouton_run = tk.Button(self, text="  Play  ", state='disabled', command=self.startModel, height=2, width=25, background=COLOR_BUTTON_BG_PLAY)
        self.bouton_run.pack(side=tk.TOP)
        
        self.bouton_quitter = tk.Button(self, text="  Quitter  ", command=self.quit, height=2, width=25)
        self.bouton_quitter.pack(side=tk.BOTTOM)
        
        self.displayFrame = tk.LabelFrame(self, text="Display options", padx=2, pady=2)
        self.var_displayCollision = tk.IntVar()
        self.var_displayCollision.set(self.model.DISPLAY_OBS_WALLCELL)
        self.var_displayPressure = tk.IntVar()
        self.var_displayPressure.set(self.model.DISPLAY_OBS_PRESSURECELL)
        self.var_displayNormals = tk.IntVar()
        self.var_displayNormals.set(self.model.DISPLAY_NORMVECTORS)
        self.button_displayCollision = tk.Checkbutton(self.displayFrame, text="Collision cells", variable=self.var_displayCollision, onvalue=1, offvalue=0)
        self.button_displayPressure = tk.Checkbutton(self.displayFrame, text="Pressure cells", variable=self.var_displayPressure, onvalue=1, offvalue=0)
        self.button_displayNormals = tk.Checkbutton(self.displayFrame, text="Surfaces vectors", variable=self.var_displayNormals, onvalue=1, offvalue=0)
        self.list_colorParticle = ttk.Combobox(self.displayFrame, values=[  "Particle state",
                                                                            "Density",
                                                                            "Pressure"])
        self.list_colorParticle.current(0)
        self.list_colorParticle.bind("<<ComboboxSelected>>", self.setColorParticle)
        self.button_displayCollision.pack(side=tk.TOP, fill=tk.X)
        self.button_displayPressure.pack(side=tk.TOP, fill=tk.X)
        self.button_displayNormals.pack(side=tk.TOP, fill=tk.X)
        self.list_colorParticle.pack()
        
        self.displayFrame.pack(side=tk.TOP)
        
        self.infoFrame = tk.LabelFrame(self, text="Simulation  information", padx=2, pady=2)
        self.particlesLabel = tk.Text(self.infoFrame)
        self.scrbarInfoLabel = tk.Scrollbar(self.infoFrame, orient=tk.VERTICAL, command=self.particlesLabel.yview)
        self.particlesLabel.config(fg='gray', width=25, height=8, yscrollcommand=self.scrbarInfoLabel.set)
        self.scrbarInfoLabel.pack(side=tk.RIGHT)
        self.particlesLabel.pack(side=tk.LEFT)
        self.infoFrame.pack(side=tk.TOP)

        
        self.infoLabel = tk.Label(self, text = "Not initialized", fg='gray', width=25)
        self.infoLabel.pack(side=tk.BOTTOM)
        
        
        

    def initialize(self):
        self.initialized = True
        self.stepNumber = 0
        self.bouton_run.configure(state='normal')
        self.updateInfo()
        self.infoLabel.configure(fg='black')
        self.drawScene()
        
        if MODEL_OUTPUT:
            self.model.createOutputFile(MODEL_OUTPUT_FILENAME)
            
        
    def startModel(self):
        self.running = True
        self.bouton_run.configure(text="  Stop  ", command=self.stopModel)
        self.updateInfo()
        self.infoLabel.configure(fg='black')
        self.runModel()
        
        
    def runModel(self):
        self.stepNumber += 1
        self.model.moveParticles()
        self.drawScene()
        self.updateInfo()
        self.infoLabel.configure(fg='black')
        if self.running:
            self.after(MOVIE_DT, self.runModel)
        
        
    def stopModel(self):
        self.running = False
        self.bouton_run.configure(text="  Play  ", command=self.startModel)
        self.updateInfo()
        self.infoLabel.configure(fg='gray')
        
        if MODEL_OUTPUT:
            f = open(MODEL_OUTPUT_FILENAME, 'r')
            data = []
            for l in f.readlines():
                if l[0] != '%':
                    data.append(l.split(' '))
                    while '' in data[-1]:
                        data[-1].remove('')
            data = np.asarray(data, dtype=float).T
            
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.plot(data[0], data[1], color=COLOR_PARTICLE_DEFAULT, label='Fluid')
            ax.plot(data[0], data[2], color=COLOR_PARTICLE_WALL_CORRECTED, label='Colliding')
            ax.plot(data[0], data[3], color=COLOR_PARTICLE_INLET, label='Inlet')
            ax.plot(data[0], data[4], color=COLOR_PARTICLE_OUTLET, label='Outlet')
            ax.plot(data[0], data[5], color=COLOR_PARTICLE_SOLID, label='Solid')
            ax.plot(data[0], data[6], color='#FF8888', label='Destroyed')
            ax.set_xlim(data[0][0], data[0][-1])
            ax.legend(loc='upper left', fontsize=12)
            fig.tight_layout()
            fig.savefig('simulation_output_Nparticles.pdf')
            
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.plot(data[0], data[7], color='blue', label='mean Vx')
            ax.plot(data[0], data[8], color='red', label='mean Vy')
            ax.set_xlim(data[0][0], data[0][-1])
            ax.legend(loc='upper left', fontsize=12)
            fig.tight_layout()
            fig.savefig('simulation_output_Velocity.pdf')


    
    def drawScene(self):
        self.model.DISPLAY_OBS_PRESSURECELL = self.var_displayPressure.get()
        self.model.DISPLAY_OBS_WALLCELL = self.var_displayCollision.get()
        self.model.DISPLAY_BOX_WALLCELL = self.var_displayCollision.get()
        self.model.DISPLAY_NORMVECTORS = self.var_displayNormals.get()
        self.model.plotTK(self.canv, self.width, self.height)
    
    
#    def openInfoPanel(self):
#        if not self.infoPanel_opened:
#            self.infoPanel_opened = True
#            self.infoPanel = tk.Toplevel(self.master)
#            self.appInfoPanel = InfoPanel(self.infoPanel)
#        else:
#            self.infoPanel_opened = False
#            self.appInfoPanel.close_window()

    def updateInfo(self):
        self.model.updateStat()
                
        self.infoLabel.configure(text = 't = {0:.6f} s'.format(self.model.simu_t))
        s = '\nMean Velocity:'
        s += '\n  Vx = {0:.3f} m.s-1'.format(self.model.mean_velocity[0])
        s += '\n  Vy = {0:.3f} m.s-1'.format(self.model.mean_velocity[1])
        s += '\nMean Pressure:'
        s += '\n  P = {0:.3f} a.p.'.format(self.model.mean_pressure)
        s += '\nMean Density:'
        s += '\n  P = {0:.3f} a.p.'.format(self.model.mean_density)
        s += '\n'
        s += 'Particles:'
        s += '\n  TOTAL  = {0}'.format(self.model.stat_Npar_fluid + self.model.stat_Npar_colliding + self.model.stat_Npar_inlet + self.model.stat_Npar_outlet + self.model.stat_Npar_solid)
        s += '\n  FLUID  = {0}'.format(self.model.stat_Npar_fluid)
        s += '\n  COLLID = {0}'.format(self.model.stat_Npar_colliding)
        s += '\n  INLET  = {0}'.format(self.model.stat_Npar_inlet)
        s += '\n  OUTLET = {0}'.format(self.model.stat_Npar_outlet)
        s += '\n  SOLID  = {0}'.format(self.model.stat_Npar_solid)
        s += '\n  DESTRO = {0}'.format(self.model.stat_Npar_destroyed)
        s += '\n'
        s += '\nAdaptative stepping:'
        s += '\n  Vmax = {0:.3f} m.s-1'.format(self.model.simu_dt_adapt_vmax)
        s += '\n  dt = {0:.4f} s'.format(self.model.simu_dt)
#        self.particlesLabel.configure(text = s)
        scrbarA, scrbarB = self.scrbarInfoLabel.get()
        self.particlesLabel.delete(1.0, tk.END)
        self.particlesLabel.insert(tk.INSERT, s)
        self.scrbarInfoLabel.set(scrbarA, scrbarB)
        
        if MODEL_OUTPUT:
            self.model.updateOutputFile(MODEL_OUTPUT_FILENAME)

    def setColorParticle(self, event):
        select = self.list_colorParticle.get()
        if select == 'Particle state':
            self.model.DISPLAY_COLORPARTICLE = 0
        elif select == 'Density':
            self.model.DISPLAY_COLORPARTICLE = 1
        elif select == 'Pressure':
            self.model.DISPLAY_COLORPARTICLE = 2
            
            
    

#class InfoPanel:
#    def __init__(self, master, model):
#        self.master = master
#        self.model = model
#        self.create_widgets()
#
#    def create_widgets(self):
#        self.frame = tk.Frame(self.master)
#        self.particlesLabel = tk.Label(self.frame, text = "info panel", fg='gray', width=50, height=19)
#        self.particlesLabel.pack(side=tk.TOP)
#        self.frame.pack()
#        
#    def updateInfo(self):
#        self.model.updateStat()
#                
#        self.infoLabel.configure(text = 't = {0:.6f} s'.format(self.model.simu_t))
#        s = 'Particles:'
#        s += '\n  TOTAL  = {0}'.format(self.model.stat_Npar_fluid + self.model.stat_Npar_colliding + self.model.stat_Npar_inlet + self.model.stat_Npar_outlet + self.model.stat_Npar_solid)
#        s += '\n  FLUID  = {0}'.format(self.model.stat_Npar_fluid)
#        s += '\n  COLLID = {0}'.format(self.model.stat_Npar_colliding)
#        s += '\n  INLET  = {0}'.format(self.model.stat_Npar_inlet)
#        s += '\n  OUTLET = {0}'.format(self.model.stat_Npar_outlet)
#        s += '\n  SOLID  = {0}'.format(self.model.stat_Npar_solid)
#        s += '\n  DESTRO = {0}'.format(self.model.stat_Npar_destroyed)
#        s += '\n'
#        s += '\nMean Velocity:'
#        s += '\n  Vx = {0:.3f} m.s-1'.format(self.model.mean_velocity[0])
#        s += '\n  Vy = {0:.3f} m.s-1'.format(self.model.mean_velocity[1])
#        s += '\nMean Pressure:'
#        s += '\n  P = {0:.3f} a.p.'.format(self.model.mean_pressure)
#        s += '\n'
#        s += '\nAdaptative stepping:'
#        s += '\n  Vmax = {0:.3f} m.s-1'.format(self.model.simu_dt_adapt_vmax)
#        s += '\n  dt = {0:.4f} s'.format(self.model.simu_dt)
#        self.particlesLabel.configure(text = s)
#    
#    def close_windows(self):
#        self.destroy()
        
        
        

if __name__ == "__main__":
    app = AppliCanevas()
    app.title("Gridded SPH simulation")
    app.mainloop()
    
    
    