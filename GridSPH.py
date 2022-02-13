#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:15:00 2020

@author: sebastien
"""

import time
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, CirclePolygon, Circle
import tkinter as tk
mp.rcParams['font.family'] = 'serif'
rand = np.random
rand.seed(42)


CELL_FLUID = 0 # fluid cell
CELL_WALL = 1 # fluid cell with half plan movement not permitted (typically, along a surface)
CELL_SOLID = 2 # solid cell
CELL_WALL_THICKNESS = 0.1
CELL_PRESSURE_THICKNESS = 0.3

PARTICLE_SIZE = 3
PARTICLE_STATE_NORMAL = 0
PARTICLE_STATE_WALL_CORRECTED = 1
PARTICLE_STATE_SOLID = 2
PARTICLE_STATE_INLET = 3
PARTICLE_STATE_OUTLET = 4
PARTICLE_STATE_TODESTROY = 5

CALCUL_INTERACTIONS = True
MODEL_OUTPUT_FILENAME = 'simulation_output.dat'


COLOR_INLET = '#CCFFCC'
COLOR_OUTLET = '#FFCCCC'
COLOR_CELL_FLUID = '#CCCCFF'
COLOR_CELL_WALL = '#FFEEDD'
COLOR_CELL_SOLID = '#DDDDDD'
COLOR_PARTICLE_DEFAULT = '#0000BB'
COLOR_PARTICLE_WALL_CORRECTED = '#BBBB00'
COLOR_PARTICLE_INLET = '#00BB00'
COLOR_PARTICLE_OUTLET = '#BB0000'
COLOR_PARTICLE_SOLID = '#888888'
COLOR_CELL_REF = ['#DDDDDD', '#CCCCFF', '#CCFFCC', '#FFEEDD', '#FFCCCC', '#FF8888']
                  
#DISPLAY_NORMVECTORS = True
#DISPLAY_NORMVECTOR_SIZE = 0.5
#DISPLAY_NEIGHBORREGION = False
#DISPLAY_OBS_PRESSURECELL = False
#DISPLAY_OBS_WALLCELL = True
#DISPLAY_BOX_WALLCELL = True
#DISPLAY_INLETOUTLET = True



class GridSPH:
    
    def __init__(self):
        self.DISPLAY_NORMVECTORS = 1
        self.DISPLAY_NORMVECTOR_SIZE = 0.5
        self.DISPLAY_NEIGHBORREGION = False
        self.DISPLAY_OBS_PRESSURECELL = 0
        self.DISPLAY_OBS_WALLCELL = 1
        self.DISPLAY_BOX_WALLCELL = 1
        self.DISPLAY_INLETOUTLET = True
        self.DISPLAY_COLORPARTICLE = 0
        self.domain_xmin = -1
        self.domain_xmax = 1
        self.domain_ymin = -1
        self.domain_ymax = 1
        self.domainPath = None
        self.simu_t = 0
        self.simu_dt_adaptative = True
        self.simu_dt_adapt_amax = 0
        self.simu_dt_adapt_vmax = 0
        self.simu_dt = 0.01
        self.obstacle = {'edges': [],
                         'faces': [],
                         'path': None,
                         'wallWidth': 0,
                         'wallEdges': [],
                         'wall': [],
                         'pressureWidth': 0,
                         'pressureEdges': [],
                         'pressureWall': [],
                         }
        self.inlet = None
        self.outlet  = None
        self.box = []
        self.particles = []
        self.stat_Npar_destroyed = 0
        self.stat_Npar_fluid = 0
        self.stat_Npar_colliding = 0
        self.stat_Npar_inlet = 0
        self.stat_Npar_outlet = 0
        self.stat_Npar_solid = 0
        self.mean_velocity = [0, 0]
        return None
    
    
    def addBox(self):
        self.box.append({'edges': [],
                         'faces': [],
                         'path': None,
                         'wallWidth': 0,
                         'wallEdges': [],
                         'wall': [],
                         'pressureWidth': 0,
                         'pressureEdges': [],
                         'pressureWall': [],
                         })
    
    def addParticle(self, x=0, y=0, mass=0, pressure = 0, neighborR = 0.5):
        self.particles.append({'pos': [x, y],
                               'velocity': [0, 0],
                               'acceleration': [0, 0],
                               'pressure': pressure,
                               'pressureStatic': pressure,
                               'pressureGrad': [0, 0],
                               'mass': mass,
                               'density': 0,
                               'neighborR': neighborR,
                               'neighbors': [],
                               'gaussGradI': 1.0, # amplitude of the gaussian used to calculate the spatial gradient
                               'state': PARTICLE_STATE_NORMAL,
                               })
            
    def duplicateParticle(self, old, pos):
        self.particles.append({'pos': pos,
                               'velocity': old['velocity'],
                               'acceleration': old['acceleration'],
                               'pressure': old['pressure'],
                               'pressureStatic': old['pressureStatic'],
                               'pressureGrad': old['pressureGrad'],
                               'mass': old['mass'],
                               'density': old['density'],
                               'neighborR': old['neighborR'],
                               'neighbors': old['neighbors'],
                               'gaussGradI': old['gaussGradI'],
                               'state': old['state'],
                               })
            
    def createInletOutlet(self, inletDomain, outletDomain, direction = 'toRight', velocity = 1.0):
        vx, vy = 0, 0
        if direction == 'toRight':
            vx = velocity
        elif direction == 'toLeft':
            vx = -velocity
        elif direction == 'toUp':
            vy = velocity
        elif direction == 'toDown':
            vy = -velocity
        self.inlet = {'xmin': inletDomain[0],
                      'xmax': inletDomain[1],
                      'ymin': inletDomain[2],
                      'ymax': inletDomain[3],
                      'edges': [[inletDomain[0], inletDomain[2]], [inletDomain[1], inletDomain[2]], 
                                [inletDomain[1], inletDomain[3]], [inletDomain[0], inletDomain[3]]],
                      'path': Path([[inletDomain[0], inletDomain[2]], [inletDomain[1], inletDomain[2]],
                                    [inletDomain[1], inletDomain[3]], [inletDomain[0], inletDomain[3]],
                                    [inletDomain[0], inletDomain[2]]], closed=True),
                      'direction': direction,
                      'velocity': [vx, vy],
                      }
        self.outlet = {'xmin': outletDomain[0],
                      'xmax': outletDomain[1],
                      'ymin': outletDomain[2],
                      'ymax': outletDomain[3],
                      'edges': [[outletDomain[0], outletDomain[2]], [outletDomain[1], outletDomain[2]],
                                [outletDomain[1], outletDomain[3]], [outletDomain[0], outletDomain[3]]],
                      'path': Path([[outletDomain[0], outletDomain[2]], [outletDomain[1], outletDomain[2]],
                                    [outletDomain[1], outletDomain[3]], [outletDomain[0], outletDomain[3]],
                                    [outletDomain[0], outletDomain[2]]], closed=True),
                      'direction': direction,
                      'velocity': [vx, vy],
                      }
        
        
    def addInletParticle(self, outPar):
        if self.inlet['direction'] == 'toRight':
            x, y = self.inlet['xmin'], outPar['pos'][1]
        elif self.inlet['direction'] == 'toLeft':
            x, y = self.inlet['xmax'], outPar['pos'][1]
        elif self.inlet['direction'] == 'toUp':
            x, y = outPar['pos'][0], self.inlet['ymin']
        elif self.inlet['direction'] == 'toDown':
            x, y = outPar['pos'][0], self.inlet['ymax']
        self.duplicateParticle(outPar, [x, y])
        self.particles[-1]['state'] = PARTICLE_STATE_INLET
#        self.Stat_Npar_inlet += 1
    
    
    def testObstacle(self, cx, cy, r, N):
        # definition of the edges
        for i in range(N):
            self.obstacle['edges'].append([cx + r*np.cos(-2*np.pi * i / N), cy + r*np.sin(-2*np.pi * i / N)])
            self.obstacle['faces'].append([i, (i+1)%N])
        self.obstacle['path'] = Path(self.obstacle['edges'] + [self.obstacle['edges'][0]], closed=True)
        self.obstacle['wallWidth'] = CELL_WALL_THICKNESS
        self.obstacle['pressureWidth'] = CELL_PRESSURE_THICKNESS
        self.setCellsFromObstacle()
        

    def obstacleNACA4(self, naca, N = 100, center = [0, 0], length = 1.0, rotate=0.0):
        # definition of the edges
        M = int(naca[0])
        P = int(naca[1])
        XX = int(naca[2:])
        t = XX/100.0
        m = M/100.0
        p = P/10.0
        lX = np.linspace(0, 1, num=int(N/2))
        lX = lX**2 * length
        lY = t*length/0.2*(0.2969*np.sqrt(lX/length) - 0.1260*(lX/length) - 0.3516*(lX/length)**2 + 0.2843*(lX/length)**3 - 0.1015*(lX/length)**4)
        
        def cambrure(x, l):
            if x < p*l:
                yc = l * m*x/l/p**2 * (2*p - x/l)
            else:
                yc = m*(l-x)/(1-p)**2 * (1 + x/l - 2*p)
            return yc
        cambrure = np.vectorize(cambrure)
        cY = cambrure(lX, length)
        
        def dcambrure(x, l):
            if x < p*l:
                yc = l * 2*m/p**2 * (p - x/l)
            else:
                yc = l * 2*m/(1-p)**2 * (p - x/l)
            return yc
        dcambrure = np.vectorize(dcambrure)
        dcY = dcambrure(lX, length)
        theta = np.arctan(dcY)
        
        
        lX = lX - length/2*np.ones_like(lX)
        
        
        for i in range(len(lX)-1): # extrados
            self.obstacle['edges'].append([lX[i] - lY[i]*np.sin(theta[i]), cY[i] + lY[i]*np.cos(theta[i])])
        for i in range(len(lX)-1): # intrados
            self.obstacle['edges'].append([lX[len(lX)-1 - i] + lY[len(lX)-1 - i]*np.sin(theta[len(lX)-1 - i]), cY[len(lX)-1 - i] - lY[len(lX)-1 - i]*np.cos(theta[len(lX)-1 - i])])
#            self.obstacle['edges'].append([lX[len(lX)-1 - i], -lY[len(lX)-1 - i]])
        
#        for i in range(len(lX)-1):
#            self.obstacle['edges'].append([lX[i], lY[i]])
#        for i in range(len(lX)-1):
#            self.obstacle['edges'].append([lX[len(lX)-1 - i], -lY[len(lX)-1 - i]])
            
        for i in range(len(self.obstacle['edges'])):
            x_ = (self.obstacle['edges'][i][0])*np.cos(rotate*np.pi/180) - (self.obstacle['edges'][i][1])*np.sin(rotate*np.pi/180) + center[0]
            y_ = (self.obstacle['edges'][i][0])*np.sin(rotate*np.pi/180) + (self.obstacle['edges'][i][1])*np.cos(rotate*np.pi/180) + center[1]
            self.obstacle['edges'][i] = [x_, y_]
        N = len(self.obstacle['edges'])
        for i in range(N):
            self.obstacle['faces'].append([i, (i+1)%N])
        self.obstacle['path'] = Path(self.obstacle['edges'] + [self.obstacle['edges'][0]], closed=True)
        self.obstacle['wallWidth'] = CELL_WALL_THICKNESS
        self.obstacle['pressureWidth'] = CELL_PRESSURE_THICKNESS
        self.setCellsFromObstacle()
        
        
        

    def setCellsFromObstacle(self):
        # estimation of the normal vectors, and the shift vector for defining the wall and pressure cells
        N = len(self.obstacle['edges'])
        vectEdges = [[0, 0] for x in range(N)]
        for i in range(N):
            ux = self.obstacle['edges'][self.obstacle['faces'][i][1]][0] - self.obstacle['edges'][self.obstacle['faces'][i][0]][0]
            uy = self.obstacle['edges'][self.obstacle['faces'][i][1]][1] - self.obstacle['edges'][self.obstacle['faces'][i][0]][1]
            normd = np.sqrt(ux**2 + uy**2)
            normx = -uy / normd
            normy = ux / normd
            self.obstacle['faces'][i].append([normx, normy])
            vectEdges[self.obstacle['faces'][i][0]][0] += normx
            vectEdges[self.obstacle['faces'][i][0]][1] += normy
            vectEdges[self.obstacle['faces'][i][1]][0] += normx
            vectEdges[self.obstacle['faces'][i][1]][1] += normy
            
        # defining the wall and pressure cell edges
        for i in range(N):
            ox = self.obstacle['edges'][i][0]
            oy = self.obstacle['edges'][i][1]
            l = np.sqrt(vectEdges[i][0]**2 + vectEdges[i][1]**2)
            self.obstacle['wallEdges'].append([ox + vectEdges[i][0]/l * self.obstacle['wallWidth'], oy + vectEdges[i][1]/l * self.obstacle['wallWidth']])
            self.obstacle['pressureEdges'].append([ox + vectEdges[i][0]/l * self.obstacle['pressureWidth'], oy + vectEdges[i][1]/l * self.obstacle['pressureWidth']])
        
        # defining the wall and pressure cells
        N = len(self.obstacle['faces'])
        for i in range(N):
            self.obstacle['wall'].append(Path([self.obstacle['edges'][i],
                                               self.obstacle['edges'][(i+1)%N],
                                               self.obstacle['wallEdges'][(i+1)%N],
                                               self.obstacle['wallEdges'][i],
                                               self.obstacle['edges'][i]], closed=True))
            self.obstacle['pressureWall'].append(Path([self.obstacle['edges'][i],
                                               self.obstacle['edges'][(i+1)%N],
                                               self.obstacle['pressureEdges'][(i+1)%N],
                                               self.obstacle['pressureEdges'][i],
                                               self.obstacle['edges'][i]], closed=True))
    
    
    
    def closedBox(self):
        self.box = []
        self.addBox()
        self.box[0]['edges'].append([self.domain_xmin, self.domain_ymin])
        self.box[0]['edges'].append([self.domain_xmax, self.domain_ymin])
        self.box[0]['edges'].append([self.domain_xmax, self.domain_ymax])
        self.box[0]['edges'].append([self.domain_xmin, self.domain_ymax])
        for i in range(4):
            self.box[0]['faces'].append([i, (i+1)%4])
        self.box[0]['path'] = Path(self.box[0]['edges'] + [self.box[0]['edges'][0]], closed=True)
        self.box[0]['wallWidth'] = 0.1
        self.box[0]['pressureWidth'] = 0.5
        self.setCellsFromBoxes()
    
    
    
    def horizontalBox(self):
        self.box = []
        self.addBox()
        self.box[0]['edges'].append([self.domain_xmin, self.domain_ymin])
        self.box[0]['edges'].append([self.domain_xmax, self.domain_ymin])
        self.box[0]['faces'].append([0, 1])
        self.box[0]['path'] = Path(self.box[0]['edges'] , closed=False)
        self.box[0]['wallWidth'] = 0.1
        self.addBox()
        self.box[1]['edges'].append([self.domain_xmax, self.domain_ymax])
        self.box[1]['edges'].append([self.domain_xmin, self.domain_ymax])
        self.box[1]['faces'].append([0, 1])
        self.box[1]['path'] = Path(self.box[1]['edges'] , closed=False)
        self.box[1]['wallWidth'] = 0.1
        self.setCellsFromBoxes()
    
    
    
    def setCellsFromBoxes(self):
        for b in range(len(self.box)):
            # estimation of the normal vectors, and the shift vector for defining the wall and pressure cells
            N = len(self.box[b]['edges'])
            vectEdges = [[0, 0] for x in range(N)]
            N = len(self.box[b]['faces'])
            for i in range(N):
                ux = self.box[b]['edges'][self.box[b]['faces'][i][1]][0] - self.box[b]['edges'][self.box[b]['faces'][i][0]][0]
                uy = self.box[b]['edges'][self.box[b]['faces'][i][1]][1] - self.box[b]['edges'][self.box[b]['faces'][i][0]][1]
                normd = np.sqrt(ux**2 + uy**2)
                normx = -uy / normd
                normy = ux / normd
                self.box[b]['faces'][i].append([normx, normy])
                vectEdges[self.box[b]['faces'][i][0]][0] += normx
                vectEdges[self.box[b]['faces'][i][0]][1] += normy
                vectEdges[self.box[b]['faces'][i][1]][0] += normx
                vectEdges[self.box[b]['faces'][i][1]][1] += normy
                
            # defining the wall cell edges
            N = len(self.box[b]['edges'])
            for i in range(N):
                ox = self.box[b]['edges'][i][0]
                oy = self.box[b]['edges'][i][1]
                l = np.sqrt(vectEdges[i][0]**2 + vectEdges[i][1]**2)
                self.box[b]['wallEdges'].append([ox + vectEdges[i][0]/l * self.box[b]['wallWidth'], oy + vectEdges[i][1]/l * self.box[b]['wallWidth']])
                
            # defining the wall cells
#            N = len(self.box[b]['edges'])
            for i in range(len(self.box[b]['faces'])):
                self.box[b]['wall'].append(Path([self.box[b]['edges'][i],
                                                   self.box[b]['edges'][(i+1)%N],
                                                   self.box[b]['wallEdges'][(i+1)%N],
                                                   self.box[b]['wallEdges'][i],
                                                   self.box[b]['edges'][i]], closed=True))
            
    
    
    def setFluidParticles(self, Nx, Ny = 1, domain = [-5, -3, -5, 5], p = 1.0, m = 1.0, arangement = 'random'):
        self.particles = []
        nR = np.sqrt((self.domain_xmax-self.domain_xmin)*(self.domain_ymax-self.domain_ymin)*4/(Nx*Ny))
        if arangement == 'random':
            for i in range(Nx*Ny):
                self.addParticle(x = np.random.uniform(domain[0], domain[1]), 
                                 y = np.random.uniform(domain[2], domain[3]), 
                                 pressure = p,
                                 mass = m, 
                                 neighborR = nR)
        elif arangement == 'uniform':
            gridX, gridY = np.meshgrid(np.linspace(domain[0], domain[1], num=Nx), np.linspace(domain[2], domain[3], num=Ny))
            for i in range(Nx):
                for j in range(Ny):
                    self.addParticle(x = gridX[j][i], 
                                     y = gridY[j][i], 
                                     pressure = p,
                                     mass = m, 
                                     neighborR = nR)
        elif arangement == 'uniformDisturbed':
            gridX, gridY = np.meshgrid(np.linspace(domain[0], domain[1], num=Nx), np.linspace(domain[2], domain[3], num=Ny))
            dX = abs(domain[0] - domain[1])/Nx
            dY = abs(domain[2] - domain[3])/Ny
            for i in range(Nx):
                for j in range(Ny):
                    self.addParticle(x = gridX[j][i] + np.random.normal(scale=dX/5), 
                                     y = gridY[j][i] + np.random.normal(scale=dY/5), 
                                     pressure = p,
                                     mass = m, 
                                     neighborR = nR)
    
    
    
    def updateNeighbors(self):
        for i in range(len(self.particles)):
            self.particles[i]['neighbors'] = []
        for i in range(len(self.particles)-2):
            neighborRegion = Path.circle(center=self.particles[i]['pos'], radius=self.particles[i]['neighborR'])
            for j in range(i+1,len(self.particles)):
                if j != i:
                    if neighborRegion.contains_point(self.particles[j]['pos']):
                        if not Path([self.particles[i]['pos'], self.particles[j]['pos'], self.particles[i]['pos']], closed=True).intersects_path(self.obstacle['path']):
                            self.particles[i]['neighbors'].append(j)
                            self.particles[j]['neighbors'].append(i)
             
                
                        
    def testRepulsion(self, pos1, pos2, C = -0.3):
        d2 = ((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
        d = np.sqrt(d2)
        a = C / d2
        return a * (pos2[0] - pos1[0])/d, a * (pos2[1] - pos1[1])/d
        
        

    def testScene(self):
        self.domain_xmin = -5.0
        self.domain_xmax = 5.0
        self.domain_ymin = -5.0
        self.domain_ymax = 5.0
        self.domainPath = Path([[self.domain_xmin,self.domain_ymin],
                                [self.domain_xmax,self.domain_ymin],
                                [self.domain_xmax,self.domain_ymax],
                                [self.domain_xmin,self.domain_ymax],
                                [self.domain_xmin,self.domain_ymin],
                                ], closed=True)
#        self.testObstacle(0, 0, 2.0, 50)
        self.obstacleNACA4('0015', 50, length=4)
        self.closedBox()
        
#        self.setFluidParticles(500, domain = [-4.9, -2, -3, 3], arangement = 'random')
        self.setFluidParticles(200, domain = [-4.9, 4.9, -4.9, 4.9], arangement = 'random')
        for i in range(len(self.particles)):
            self.particles[i]['velocity'] = [rand.uniform(-0.1, 0.1), rand.uniform(-0.1, 0.1)]
        i = 0
        while(i < len(self.particles)):
            if self.obstacle['path'].contains_point(self.particles[i]['pos']):
                self.particles.pop(i)
            else:
                i += 1
#            self.particles[i]['velocity'] = [1.0, 0.0]
                
                
                
    
    def uniformFlowScene(self):
        self.domain_xmin = -7.0
        self.domain_xmax = 7.0
        inletWidth = 2.0
        self.domain_ymin = -4.0
        self.domain_ymax = 4.0
        self.init_U = 1.0
        self.init_P = 1e0
        self.domainPath = Path([[self.domain_xmin,self.domain_ymin],
                                [self.domain_xmax,self.domain_ymin],
                                [self.domain_xmax,self.domain_ymax],
                                [self.domain_xmin,self.domain_ymax],
                                [self.domain_xmin,self.domain_ymin],
                                ], closed=True)
        self.testObstacle(0, 0, 1.0, 20)
#        self.obstacleNACA4('4415', 30, length=4, rotate=0.0)
        self.horizontalBox()
        
#        self.setFluidParticles(42, Ny=24, domain = [self.domain_xmin+CELL_WALL_THICKNESS/2, self.domain_xmax-CELL_WALL_THICKNESS/2, self.domain_ymin+CELL_WALL_THICKNESS/2, self.domain_ymax-CELL_WALL_THICKNESS/2], 
        self.setFluidParticles(22, Ny=12, domain = [self.domain_xmin+CELL_WALL_THICKNESS/2, self.domain_xmax-CELL_WALL_THICKNESS/2, self.domain_ymin+CELL_WALL_THICKNESS/2, self.domain_ymax-CELL_WALL_THICKNESS/2], 
                               arangement = 'uniform', # 'uniform', 'random', 'uniformDisturbed'
                               m = 1.0,
                               p = self.init_P) 
        for i in range(len(self.particles)):
            self.particles[i]['velocity'] = [self.init_U, 0]
        i = 0
        while(i < len(self.particles)):
            if self.obstacle['path'].contains_point(self.particles[i]['pos']):
                self.particles.pop(i)
            else:
                i += 1
        
        if CALCUL_INTERACTIONS:
            self.updateNeighbors()
        
            I = np.zeros(len(self.particles))
            for i in range(len(self.particles)):
                I[i] += 1.0#self.particles[i]['pressure']
                for j in range(len(self.particles[i]['neighbors'])):
                    dist = np.sqrt((self.particles[i]['pos'][0] - self.particles[self.particles[i]['neighbors'][j]]['pos'][0])**2 + (self.particles[i]['pos'][1] - self.particles[self.particles[i]['neighbors'][j]]['pos'][1])**2)
#                    I[i] += self.particles[self.particles[i]['neighbors'][j]]['pressure'] * np.exp(-dist**2 / (2*(self.particles[i]['neighborR']/3)**2))
                    I[i] += np.exp(-dist**2 / (2*(self.particles[i]['neighborR']/3)**2))
            Imean = np.average(I)
            print("Imean = {0}".format(Imean))
            for i in range(len(self.particles)):
                if Imean != 0:
                    self.particles[i]['gaussGradI'] = 1.0/Imean
                else:
                    self.particles[i]['gaussGradI'] = 0.0
            
            
        self.createInletOutlet([self.domain_xmin, self.domain_xmin+inletWidth, self.domain_ymin, self.domain_ymax], [self.domain_xmax-inletWidth, self.domain_xmax, self.domain_ymin, self.domain_ymax], direction = 'toRight', velocity = self.init_U)
        for i in range(len(self.particles)):
            if self.inlet['path'].contains_point(self.particles[i]['pos']):
                self.particles[i]['state'] = PARTICLE_STATE_INLET
                self.particles[i]['velocity'] = self.inlet['velocity']
            elif self.outlet['path'].contains_point(self.particles[i]['pos']):
                self.particles[i]['state'] = PARTICLE_STATE_OUTLET
                self.particles[i]['velocity'] = self.outlet['velocity']
        
        
    def gauss2d(self, I, s, d):
        return I*np.exp(-d**2/(2*s**2))
    
    def gradP(self):
        for i in range(len(self.particles)):
            gP = [0, 0]
            for j in range(len(self.particles[i]['neighbors'])):
                dist = np.sqrt( (self.particles[self.particles[i]['neighbors'][j]]['pos'][0] - self.particles[i]['pos'][0])**2 + (self.particles[self.particles[i]['neighbors'][j]]['pos'][1] - self.particles[i]['pos'][1])**2 )
                dist = max(dist, CELL_WALL_THICKNESS)
                dP = (self.particles[self.particles[i]['neighbors'][j]]['pressure'] - self.particles[i]['pressure']) / dist
                gP[0] += dP * (self.particles[self.particles[i]['neighbors'][j]]['pos'][0] - self.particles[i]['pos'][0])
                gP[1] += dP * (self.particles[self.particles[i]['neighbors'][j]]['pos'][1] - self.particles[i]['pos'][1])
            self.particles[i]['pressureGrad'] = [gP[0], gP[1]]
    
    def calcP(self):
        P = np.zeros(len(self.particles))
        for i in range(len(self.particles)):
            figauss = self.gauss2d(self.particles[i]['gaussGradI'], self.particles[i]['neighborR']/3, 0)
            P[i] = self.particles[i]['pressureStatic'] * figauss
            for j in range(len(self.particles[i]['neighbors'])):
                dist = np.sqrt( (self.particles[self.particles[i]['neighbors'][j]]['pos'][0] - self.particles[i]['pos'][0])**2 + (self.particles[self.particles[i]['neighbors'][j]]['pos'][1] - self.particles[i]['pos'][1])**2 )
                fjgauss = self.gauss2d(self.particles[i]['gaussGradI'], self.particles[i]['neighborR']/3, dist)
                P[i] += self.particles[self.particles[i]['neighbors'][j]]['pressureStatic'] * fjgauss
        for i in range(len(self.particles)):
            self.particles[i]['pressure'] = P[i]
        
    def calcDensity(self):
        m = np.zeros(len(self.particles))
        for i in range(len(self.particles)):
            figauss = self.gauss2d(self.particles[i]['gaussGradI'], self.particles[i]['neighborR']/3, 0)
            m[i] = self.particles[i]['mass'] * figauss
            for j in range(len(self.particles[i]['neighbors'])):
                dist = np.sqrt( (self.particles[self.particles[i]['neighbors'][j]]['pos'][0] - self.particles[i]['pos'][0])**2 + (self.particles[self.particles[i]['neighbors'][j]]['pos'][1] - self.particles[i]['pos'][1])**2 )
                fjgauss = self.gauss2d(self.particles[i]['gaussGradI'], self.particles[i]['neighborR']/3, dist)
                m[i] += self.particles[self.particles[i]['neighbors'][j]]['mass'] * fjgauss
        for i in range(len(self.particles)):
            self.particles[i]['density'] = m[i] * np.pi*(self.particles[i]['neighborR']/3)**2
        
    def moveParticles(self):
        t0 = time.time()*1000.0
        
        dvx = np.zeros(len(self.particles))
        dvy = np.zeros(len(self.particles))
        v = np.zeros(len(self.particles))
        a = np.zeros(len(self.particles))
        
        # ACCELERATION -> VELOCITY
        if CALCUL_INTERACTIONS:
            self.updateNeighbors()
            self.calcDensity()
            self.calcP()
            self.gradP()
            
        for i in range(len(self.particles)):
            self.particles[i]['acceleration'] = [-self.particles[i]['pressureGrad'][0], -self.particles[i]['pressureGrad'][1]]
#            self.particles[i]['acceleration'] = [0, 0]
#            if CALCUL_INTERACTIONS:
#                for j in range(len(self.particles[i]['neighbors'])):
#                    ax, ay = self.testRepulsion(self.particles[i]['pos'], self.particles[self.particles[i]['neighbors'][j]]['pos'])
#                    self.particles[i]['acceleration'][0] += ax
#                    self.particles[i]['acceleration'][1] += ay
            
            
            if self.particles[i]['state'] != PARTICLE_STATE_INLET and self.particles[i]['state'] != PARTICLE_STATE_OUTLET:
                a[i] = np.sqrt(self.particles[i]['acceleration'][0]**2 + self.particles[i]['acceleration'][1]**2)
                v[i] = np.sqrt(self.particles[i]['velocity'][0]**2 + self.particles[i]['velocity'][1]**2)
                dvx[i] = self.particles[i]['acceleration'][0]*self.simu_dt
                dvy[i] = self.particles[i]['acceleration'][1]*self.simu_dt
            
        if self.simu_dt_adaptative:
            self.simu_dt_adapt_vmax = max(v)
            self.simu_dt_adapt_amax = max(a)
            D = self.simu_dt_adapt_vmax**2 + 4*self.simu_dt_adapt_amax*CELL_WALL_THICKNESS
            if self.simu_dt_adapt_amax == 0:
                new_dt = self.simu_dt
            else:
                new_dt = (self.simu_dt_adapt_vmax - np.sqrt(D)) / (-2*self.simu_dt_adapt_amax)
            dvx = dvx/self.simu_dt*new_dt
            dvy = dvy/self.simu_dt*new_dt
            self.simu_dt = new_dt
#            print("new_dt = {0}".format(self.simu_dt_adapt_amax))
            
        for i in range(len(self.particles)):
            if self.particles[i]['state'] != PARTICLE_STATE_INLET and self.particles[i]['state'] != PARTICLE_STATE_OUTLET:
                self.particles[i]['velocity'] = [
                        self.particles[i]['velocity'][0] + dvx[i],#self.particles[i]['acceleration'][0]*self.simu_dt,
                        self.particles[i]['velocity'][1] + dvy[i]#self.particles[i]['acceleration'][1]*self.simu_dt
                        ]
                    
        # VELOCITY -> POSITION
        for i in range(len(self.particles)):
            if self.particles[i]['state'] != PARTICLE_STATE_SOLID:
                if (self.obstacle['path'].contains_point(self.particles[i]['pos']) or
                    not(self.domainPath.contains_point(self.particles[i]['pos']))):
#                    print("particle {0} frozen!".format(i))
                    self.particles[i]['velocity'] = [0, 0]
                    self.particles[i]['state'] = PARTICLE_STATE_SOLID
                    
                elif self.particles[i]['state'] != PARTICLE_STATE_INLET and self.particles[i]['state'] != PARTICLE_STATE_OUTLET:
                    self.particles[i]['state'] = PARTICLE_STATE_NORMAL
                    for j in range(len(self.obstacle['wall'])):
                        if self.obstacle['wall'][j].contains_point(self.particles[i]['pos']):
#                            print("particle {0} in obstacle detected!".format(i))
                            ps = (self.particles[i]['velocity'][0]*self.obstacle['faces'][j][2][0] + self.particles[i]['velocity'][1]*self.obstacle['faces'][j][2][1]) / (self.obstacle['faces'][j][2][0]**2 + self.obstacle['faces'][j][2][1]**2)
                            if ps < 0:
                                self.particles[i]['velocity'][0] -= ps*self.obstacle['faces'][j][2][0]
                                self.particles[i]['velocity'][1] -= ps*self.obstacle['faces'][j][2][1]
                                self.particles[i]['state'] = PARTICLE_STATE_WALL_CORRECTED
                            
                    
                    for b in range(len(self.box)):
                        for j in range(len(self.box[b]['wall'])):
                            if self.box[b]['wall'][j].contains_point(self.particles[i]['pos']):
#                                print("particle {0} in box {1} detected!".format(i, b))
                                ps = (self.particles[i]['velocity'][0]*self.box[b]['faces'][j][2][0] + self.particles[i]['velocity'][1]*self.box[b]['faces'][j][2][1]) / (self.box[b]['faces'][j][2][0]**2 + self.box[b]['faces'][j][2][1]**2)
                                if ps < 0:
                                    self.particles[i]['velocity'][0] -= ps*self.box[b]['faces'][j][2][0]
                                    self.particles[i]['velocity'][1] -= ps*self.box[b]['faces'][j][2][1]
                                    self.particles[i]['state'] = PARTICLE_STATE_WALL_CORRECTED
                            
                self.particles[i]['pos'] = [
                        self.particles[i]['pos'][0] + self.particles[i]['velocity'][0]*self.simu_dt,
                        self.particles[i]['pos'][1] + self.particles[i]['velocity'][1]*self.simu_dt
                        ]
                # changing particle state for inlet-outlet
                if not self.inlet is None:
                    if self.particles[i]['state'] == PARTICLE_STATE_INLET:
                        if not self.inlet['path'].contains_point(self.particles[i]['pos']):
                            self.particles[i]['state'] = PARTICLE_STATE_NORMAL
                            self.addInletParticle(self.particles[i])
                            
                    elif self.particles[i]['state'] == PARTICLE_STATE_OUTLET:
                        if not self.outlet['path'].contains_point(self.particles[i]['pos']):
#                            print('Destroying call for particle {0}'.format(i))
                            self.particles[i]['state'] = PARTICLE_STATE_TODESTROY
                    else:  #if self.particles[i]['state'] != PARTICLE_STATE_INLET or self.particles[i]['state'] != PARTICLE_STATE_OUTLET:
                        if self.inlet['path'].contains_point(self.particles[i]['pos']):
#                            self.particles[i]['state'] = PARTICLE_STATE_INLET
                            self.particles[i]['velocity'] = self.inlet['velocity']
                        elif self.outlet['path'].contains_point(self.particles[i]['pos']):
                            self.particles[i]['state'] = PARTICLE_STATE_OUTLET
                            self.particles[i]['velocity'] = self.outlet['velocity']
                    
        i = 0
        while(i < len(self.particles)):
            if self.particles[i]['state'] == PARTICLE_STATE_TODESTROY:
                self.particles.pop(i)
                self.stat_Npar_destroyed += 1
            else:
                i += 1
                
        self.simu_t += self.simu_dt
                
        dt = int(time.time()*1000.0 - t0)
        print("Time for moving the particles: {0} ms".format(dt))
    
    
    
    def updateStat(self):
        self.stat_Npar_fluid = 0
        self.stat_Npar_colliding = 0
        self.stat_Npar_inlet = 0
        self.stat_Npar_outlet = 0
        self.stat_Npar_solid = 0
        
        self.mean_velocity = [0, 0]
        self.mean_pressure = 0.0
        self.mean_density = 0.0
        
        for i in range(len(self.particles)):
            if self.particles[i]['state'] == PARTICLE_STATE_NORMAL: 
                self.stat_Npar_fluid += 1
                self.mean_velocity[0] += self.particles[i]['velocity'][0]
                self.mean_velocity[1] += self.particles[i]['velocity'][1]
                self.mean_pressure += self.particles[i]['pressure']
                self.mean_density += self.particles[i]['density']
            if self.particles[i]['state'] == PARTICLE_STATE_WALL_CORRECTED: 
                self.stat_Npar_colliding += 1
                self.mean_velocity[0] += self.particles[i]['velocity'][0]
                self.mean_velocity[1] += self.particles[i]['velocity'][1]
                self.mean_pressure += self.particles[i]['pressure']
                self.mean_density += self.particles[i]['density']
            if self.particles[i]['state'] == PARTICLE_STATE_INLET: 
                self.stat_Npar_inlet += 1
            if self.particles[i]['state'] == PARTICLE_STATE_OUTLET: 
                self.stat_Npar_outlet += 1
            if self.particles[i]['state'] == PARTICLE_STATE_SOLID: 
                self.stat_Npar_solid += 1
        self.mean_velocity[0] /= (self.stat_Npar_fluid+self.stat_Npar_colliding)
        self.mean_velocity[1] /= (self.stat_Npar_fluid+self.stat_Npar_colliding)
        self.mean_pressure /= (self.stat_Npar_fluid+self.stat_Npar_colliding)
        self.mean_density /= (self.stat_Npar_fluid+self.stat_Npar_colliding)
        
        
        
        
    def createOutputFile(self, filename = MODEL_OUTPUT_FILENAME):
        fout = open(filename, 'w')
        fout.write('% Simulation output data\n')
        fout.write('% ---------------------------------------------------------------------------------------------\n')
        fout.write('% Time[s]   Nfluid   Ncoll   Ninlet   Noutlet   Nsolid   Ndestr   MeanVx[m.s-1]   MeanVy[m.s-1]\n')
        fout.write('% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
        fout.close()
        
    def updateOutputFile(self, filename = MODEL_OUTPUT_FILENAME):
        fout = open(filename, 'a')
        fout.write('{0:10.6f}  {1:5}  {2:5}  {3:5}  {4:5}  {5:5}  {6:5}  {7:6.3f}  {8:6.3f}\n'.format(
                self.simu_t,
                self.stat_Npar_fluid,
                self.stat_Npar_colliding,
                self.stat_Npar_inlet,
                self.stat_Npar_outlet,
                self.stat_Npar_solid,
                self.stat_Npar_destroyed,
                self.mean_velocity[0], self.mean_velocity[1]))
        fout.close()
        
    
    
    def plotMPL(self):
        fig, ax = plt.subplots(figsize=(8, 8))
        plotDX = self.domain_xmax-self.domain_xmin
        plotDY = self.domain_ymax-self.domain_ymin
        ax.set_xlim(self.domain_xmin - plotDX/20, self.domain_xmax + plotDX/20)
        ax.set_ylim(self.domain_ymin - plotDY/20, self.domain_ymax + plotDY/20)
                
        for i in range(len(self.obstacle['pressureWall'])):
            ax.add_patch(PathPatch(self.obstacle['pressureWall'][i], edgecolor='black', facecolor=COLOR_CELL_FLUID, fill=True, linewidth=0.2))
        for i in range(len(self.obstacle['wall'])):
            ax.add_patch(PathPatch(self.obstacle['wall'][i], edgecolor='black', facecolor=COLOR_CELL_WALL, fill=True, linewidth=0.2))   
        ax.add_patch(PathPatch(self.obstacle['path'], edgecolor='black', hatch='//', fill=False, linewidth=1.0))
    
        for b in range(len(self.box)):
            for i in range(len(self.box[b]['wall'])):
                ax.add_patch(PathPatch(self.box[b]['wall'][i], edgecolor='black', facecolor=COLOR_CELL_WALL, fill=True, linewidth=0.2))   
            ax.add_patch(PathPatch(self.box[b]['path'], edgecolor='black', fill=False, linewidth=1.0))
    
        for i in range(len(self.particles)):
            ax.add_patch(Circle(self.particles[i]['pos'], 0.025, color='blue'))
        
        fig.tight_layout()
        fig.savefig('plotScene.pdf')

    
    def colorScale(self, xlist):
    #print " building color scale from : %s"%str(xlist)
        r, g, b = 0,0,0
        color = []
        emin = min(xlist)
        emax = max(xlist)
        D = float(emax-emin)
    #    print "emin = %.2f"%emin
    #    print "emax = %.2f"%emax
    #    print "D    = %.2f"%D
        for l in range(len(xlist)):
            e = xlist[l]-emin
            if e < 0: 
                r=0
                g=0
                b=255
            elif e >= 0 and e < D/4.0:
                r=0
                g=int(4*e/D*255)
                b=255
            elif e >= D/4.0 and e < 2*D/4.0:
                r=0
                g=255
                b=int((1-(4*e-D)/D)*255)
            elif e >= 2*D/4.0 and e < 3*D/4.0:
                r=int((4*e-2*D)/D*255)
                g=255
                b=0
            elif e >= 3*D/4.0 and e < D:
                r=255
                g=int((1-(4*e-3*D)/D)*255)
                b=0
            else: 
                r=255
                g=0
                b=0
            color.append("#" + ("{:0>2x}{:0>2x}{:0>2x}".format(r, g, b)))
        
        #print "color scale is : %s"%str(color)
        return color
    
    
    def plotTK(self, canvas, width=400, height=400):
#        t0 = time.time()*1000.0
        canvas.delete("all")

        def coordX(x_):
            return 1+(x_ - self.domain_xmin) * width / (self.domain_xmax - self.domain_xmin)
        def coordY(y_):
            return height - (y_ - self.domain_ymin) * height / (self.domain_ymax - self.domain_ymin)
        def coordEdges(edg):
            coords = np.zeros_like(edg, dtype=float).flatten()
            for i in range(len(edg)):
                coords[2*i] = (coordX(edg[i][0]))
                coords[2*i+1] = (coordY(edg[i][1]))
            return coords.tolist()
        
        if self.DISPLAY_INLETOUTLET:
            if not self.inlet is None:
                canvas.create_polygon(coordEdges(self.inlet['edges']), fill=COLOR_INLET, outline=COLOR_INLET, width=1)
            if not self.outlet is None:
                canvas.create_polygon(coordEdges(self.outlet['edges']), fill=COLOR_OUTLET, outline=COLOR_OUTLET, width=1)
            
        
        canvas.create_polygon(coordEdges(self.obstacle['edges']), fill=COLOR_CELL_SOLID, outline="black", width=1)
        N = len(self.obstacle['faces'])
        if self.DISPLAY_OBS_PRESSURECELL:
            for i in range(len(self.obstacle['pressureWall'])):
                canvas.create_polygon(coordEdges([self.obstacle['edges'][i],
                                                   self.obstacle['edges'][(i+1)%N],
                                                   self.obstacle['pressureEdges'][(i+1)%N],
                                                   self.obstacle['pressureEdges'][i]]), fill=COLOR_CELL_FLUID, outline="black", width=0.2)
        if self.DISPLAY_OBS_WALLCELL:
            for i in range(len(self.obstacle['wall'])):
                canvas.create_polygon(coordEdges([self.obstacle['edges'][i],
                                                   self.obstacle['edges'][(i+1)%N],
                                                   self.obstacle['wallEdges'][(i+1)%N],
                                                   self.obstacle['wallEdges'][i]]), fill=COLOR_CELL_WALL, outline="black", width=0.2)
        if self.DISPLAY_NORMVECTORS:
            for i in range(N):
                canvas.create_line(coordX((self.obstacle['edges'][i][0]+self.obstacle['edges'][(i+1)%N][0])/2), coordY((self.obstacle['edges'][i][1]+self.obstacle['edges'][(i+1)%N][1])/2), 
                                   coordX((self.obstacle['edges'][i][0]+self.obstacle['edges'][(i+1)%N][0])/2 + self.DISPLAY_NORMVECTOR_SIZE*self.obstacle['faces'][i][2][0]), coordY((self.obstacle['edges'][i][1]+self.obstacle['edges'][(i+1)%N][1])/2 + self.DISPLAY_NORMVECTOR_SIZE*self.obstacle['faces'][i][2][1]), 
                                   arrow=tk.LAST)
        
        for b in range(len(self.box)):
            N = len(self.box[b]['edges'])
            for i in range(len(self.box[b]['wall'])):
                if self.DISPLAY_BOX_WALLCELL:
                    canvas.create_polygon(coordEdges([self.box[b]['edges'][i],
                                                      self.box[b]['edges'][(i+1)%N],
                                                      self.box[b]['wallEdges'][(i+1)%N],
                                                      self.box[b]['wallEdges'][i]]), fill=COLOR_CELL_WALL, outline="black", width=0.2)
                canvas.create_line(coordEdges([self.box[b]['edges'][i],
                                                  self.box[b]['edges'][(i+1)%N]]), fill='black')
                if self.DISPLAY_NORMVECTORS:
                    canvas.create_line(coordX((self.box[b]['edges'][i][0]+self.box[b]['edges'][(i+1)%N][0])/2), coordY((self.box[b]['edges'][i][1]+self.box[b]['edges'][(i+1)%N][1])/2), 
                                       coordX((self.box[b]['edges'][i][0]+self.box[b]['edges'][(i+1)%N][0])/2 + self.DISPLAY_NORMVECTOR_SIZE*self.box[b]['faces'][i][2][0]), coordY((self.box[b]['edges'][i][1]+self.box[b]['edges'][(i+1)%N][1])/2 + self.DISPLAY_NORMVECTOR_SIZE*self.box[b]['faces'][i][2][1]), 
                                       arrow=tk.LAST)
    
        if self.DISPLAY_COLORPARTICLE == 1:
            Xlist = []
            for i in range(len(self.particles)):
                if not(self.particles[i]['state'] == PARTICLE_STATE_SOLID or
                   self.particles[i]['state'] == PARTICLE_STATE_INLET or
                   self.particles[i]['state'] == PARTICLE_STATE_OUTLET):
                    Xlist.append(self.particles[i]['density'])
            colors = self.colorScale(Xlist)
        
        if self.DISPLAY_COLORPARTICLE == 2:
            Xlist = []
            for i in range(len(self.particles)):
                if not(self.particles[i]['state'] == PARTICLE_STATE_SOLID or
                   self.particles[i]['state'] == PARTICLE_STATE_INLET or
                   self.particles[i]['state'] == PARTICLE_STATE_OUTLET):
                    Xlist.append(self.particles[i]['pressure'])
            colors = self.colorScale(Xlist)
        
        j = 0
        for i in range(len(self.particles)):
            if self.DISPLAY_COLORPARTICLE == 0:
                color = COLOR_PARTICLE_DEFAULT
                if self.particles[i]['state'] == PARTICLE_STATE_WALL_CORRECTED:
                    color = COLOR_PARTICLE_WALL_CORRECTED
                elif self.particles[i]['state'] == PARTICLE_STATE_SOLID:
                    color = COLOR_PARTICLE_SOLID
                elif self.particles[i]['state'] == PARTICLE_STATE_INLET:
                    color = COLOR_PARTICLE_INLET
                elif self.particles[i]['state'] == PARTICLE_STATE_OUTLET:
                    color = COLOR_PARTICLE_OUTLET
            elif self.DISPLAY_COLORPARTICLE == 1 or self.DISPLAY_COLORPARTICLE == 2:
                if(self.particles[i]['state'] == PARTICLE_STATE_SOLID or
                   self.particles[i]['state'] == PARTICLE_STATE_INLET or
                   self.particles[i]['state'] == PARTICLE_STATE_OUTLET):
                    color = COLOR_PARTICLE_SOLID
                else:
                    color = colors[j]
                    j += 1
                
            
                
            canvas.create_oval(coordX(self.particles[i]['pos'][0])-PARTICLE_SIZE, coordY(self.particles[i]['pos'][1])-PARTICLE_SIZE, coordX(self.particles[i]['pos'][0])+PARTICLE_SIZE, coordY(self.particles[i]['pos'][1])+PARTICLE_SIZE, outline=color, fill=color)
            
            if self.DISPLAY_NEIGHBORREGION:
                canvas.create_oval(coordX(self.particles[i]['pos'][0]-self.particles[i]['neighborR']), coordY(self.particles[i]['pos'][1]-self.particles[i]['neighborR']), coordX(self.particles[i]['pos'][0]+self.particles[i]['neighborR']), coordY(self.particles[i]['pos'][1]+self.particles[i]['neighborR']), dash=(2,2), outline='black', fill='')

        
#        dt = int(time.time()*1000.0 - t0)
#        print("Time for updating the scene: {0} ms".format(dt))
    
    


