# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 20:54:11 2024

@author: liuxi
"""

import gudhi
import numpy as np


class Acupuncture():
    def __init__(self,filename):
        self.point = []  # points of the data
        self.short_chain_point = []  # chain points up to the pt3dclear
        self.long_chain_point = []   # chain points up to the file
        self.short_chain_number = 0
        self.long_chain_number = 0
        self.neighbor_distance = []
        self.prepare(filename)
        self.beta = 0   # index for branching loop
        self.alpha = 0  # index for branching angle
        self.c = 1      # index for curvature
        self.e = 0      # index for entanglement
        self.graph = {'V':[],'E':[]}
        
        
    
    @classmethod
    def distance_of_points(cls,p1,p2):
        temp = pow(p1[0]-p2[0],2) + pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2)
        return pow(temp,0.5)
    
    @classmethod
    def get_angle(cls,v1,v2):
        one = [v1[0],v1[1],v1[2]]
        two = [v2[0],v2[1],v2[2]]
        temp1 = one[0]*two[0] + one[1] * two[1] + one[2]*two[2]
        len1 = pow( pow(one[0],2) + pow(one[1],2) + pow(one[2],2),0.5)
        len2 = pow( pow(two[0],2) + pow(two[1],2) + pow(two[2],2),0.5)
        
        v = round(temp1/(len1*len2),6)
        if v<-1 or v>1:
            print(v,'error')
            return
        angle = np.arccos(v)
        t = angle*180.0/np.pi
        return t  
    
    def prepare(self,filename):
        f = open(filename)
        contents = f.readlines()
        f.close()
        
        for i in range(len(contents)):
            if contents[i][0:11]=='  pt3dclear':
                self.short_chain_number = self.short_chain_number + 1
            if contents[i][0:15]=='create filament':
                self.long_chain_number = self.long_chain_number + 1
        self.short_chain_point = [ [] for i in range(self.short_chain_number) ]
        self.long_chain_point = [ [] for i in range(self.long_chain_number) ]
        
        
        s_n = -1
        l_n = -1
        for i in range(len(contents)):
            if contents[i][0:9]=='  pt3dadd':
                line = contents[i][10::].split(',')
                self.short_chain_point[s_n].append([ float(line[0]),float(line[1]),float(line[2])])
                self.long_chain_point[l_n].append([ float(line[0]),float(line[1]),float(line[2])])
                if [ float(line[0]),float(line[1]),float(line[2])] not in self.point:
                    self.point.append([ float(line[0]),float(line[1]),float(line[2])])
            if contents[i][0:11]=='  pt3dclear':
                s_n = s_n + 1
            if contents[i][0:15]=='create filament':
                l_n = l_n + 1
        
        for i in range(self.short_chain_number):
            if len(self.short_chain_point[i])<2:
                continue
            else:
                for j in range(1,len(self.short_chain_point[i])):
                    dis = Acupuncture.distance_of_points(self.short_chain_point[i][j], self.short_chain_point[i][j-1])
                    self.neighbor_distance.append(dis)
        self.neighbor_distance = sorted(self.neighbor_distance)
                
    def compute_index_beta(self,min_length=6):
        # betti1 index for branching loop
        point = []
        for item in self.short_chain_point:
            if len(item)<min_length:
                continue
            else:
                for p in item:
                    if p not in point:
                        point.append(p)
        rips_complex = gudhi.RipsComplex(points=point,max_edge_length=self.neighbor_distance[-1])
        simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
        diag = simplex_tree.persistence()
        for item in diag:
            if item[0]==1 and item[1][1]==float('inf'):
                self.beta = self.beta + 1
    
    def get_graph(self):
        for chain in self.short_chain_point:
            if len(chain)<3:
                continue
            if chain[0] not in self.graph['V']:
                self.graph['V'].append(chain[0])
            if chain[-1] not in self.graph['V']:
                self.graph['V'].append(chain[-1])
            index1 = self.graph['V'].index(chain[0])
            index2 = self.graph['V'].index(chain[-1])
            m = min(index1,index2)
            M = max(index1,index2)
            if [m,M] not in self.graph['E']:
                self.graph['E'].append([m,M])
    
    def compute_index_alpha(self):
        # alpha index for branching angle
        self.get_graph()
        for i in range(len(self.graph['V'])):
            d = 0
            for item in self.graph['E']:
                if i in item:
                    d = d + 1
            if d>2:
                self.alpha = self.alpha + 1
        
    def get_curvature_hop_n(self,point,N):
        d = [ [] for i in range(len(point)) ]
        for i in range(len(point)):
            p = point[i]
            for j in range(N,len(p)-N,N):
                dis1 = Acupuncture.distance_of_points(p[j], p[j-N])
                dis2 = Acupuncture.distance_of_points(p[j], p[j+N])
                arc = dis1 + dis2
                v1 = [ p[j][0]-p[j-N][0], p[j][1]-p[j-N][1], p[j][2]-p[j-N][2]  ]
                v2 = [ p[j+N][0]-p[j][0], p[j+N][1]-p[j][1], p[j+N][2]-p[j][2]  ]
                angle = Acupuncture.get_angle(v1, v2)
                curvature = angle/arc
                d[i].append(curvature)
        res = []
        for d1 in d:
            if len(d1)>0:
                res.append(max(d1))
        return np.mean(res)
    
    def compute_index_c(self,k):
        def remove_repeat(point):
            res = [ [] for i in range(len(point)) ]
            for i in range(len(point)):
                p = point[i]
                for j in range(len(p)):
                    if p[j] in res[i]:
                        continue
                    else:
                        res[i].append( [p[j][0],p[j][1],p[j][2]] )
            return res
        # c = c1 * c2 * ... * ck
        point = remove_repeat(self.long_chain_point)
        for i in range(1,k+1):
            self.c = self.c * self.get_curvature_hop_n(point, i)
       
    @classmethod
    def compute_index_e(cls,vessel,nerve,r=100):
        # e index for entanglement between vessel and nerve/f-nerve
        dis_ls = []
        for i in range(len(vessel.point)):
            for j in range(len(nerve.point)):
                dis = Acupuncture.distance_of_points(vessel.point[i], nerve.point[j])
                dis_ls.append(dis)
        e = 0
        for item in dis_ls:
            if item<=r:
                e = e + 1
        return e
    
