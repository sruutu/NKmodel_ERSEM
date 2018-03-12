# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 12:18:16 2015

"""
__author__ =  'Sampsa Ruutu'
__version__=  'ERSEM final'

import random
from itertools import combinations, product
from math import exp, sqrt
import numpy as np
from bitstring import BitArray
import matplotlib.pyplot as plt
import scipy.stats
import time
import pickle
import pandas as pd
import matplotlib.cm as cm
import networkx as nx

class Unit:
    def __init__(self,org,parts,string,landscape,z,num,strategy):
        self.num = num # Number to identify a unit within the organization
        self.string = string # Determines the location on the landscape 
        self.landscape = landscape # Landscape instance
        self.z = z # Maximum search distance
        self.strategy = strategy # Determines how the unit performs search
        self.org = org # Organization that the unit is part of
        self.n = len(self.string) # Length of the search string
        self.localoptimumflag=False # Has the local optimum been reached?
        self.search_list = []
        self.search_list_start = 0
        
        if self.strategy=='partition':
            '''
            Assigns the unit to a part based on the unit's number. Low/high-
            values specify the indexes on the string to be searched.
            '''
            self.assigned_part = self.num%parts
            self.elements_per_part = int(self.n/parts)
            self.low = int(self.assigned_part*self.elements_per_part)
            self.high = self.low+self.elements_per_part
            
        else:
            self.low = 0
            self.high = self.n        
        
        self.updatePosition(string)

    def updatePosition(self,string):
        self.string =string
        self.fitness = self.fitnessAt(string)
                             
    def localsearch(self):  
        if self.localoptimumflag:
            return self.string, self.fitness
                
        # Create a list of possible combinations of string indexes to change
        if(self.search_list == []):
            for L in range(1, self.z+1):
                for subset in combinations(range(self.low,self.high), L):
                    self.search_list.append(subset)
            random.seed()
            random.shuffle(self.search_list)
        
        # Substrings based on partition
        str_begin = self.string[0:self.low]
        str_end = self.string[self.high:]
        offset = self.low    
        
        # Go through the items on the search list. 
        maxsearches = len(self.search_list)
        for i in range(maxsearches):            
            if i+self.search_list_start == len(self.search_list):
                self.localoptimumflag = True
                return self.string, self.fitness
            indexes=self.search_list[i+self.search_list_start]            
            localsearchstring = list(self.string[self.low:self.high])
            for j in range(len(indexes)):
                if localsearchstring[indexes[j]-offset] == '0':
                    localsearchstring[indexes[j]-offset] = '1'
                else:
                    localsearchstring[indexes[j]-offset] = '0'
            localsearchstring = ''.join(localsearchstring)
            newstring = str_begin + localsearchstring + str_end

            # Return the first string and fitness that improves performance        
            if self.fitnessAt(self.string) < self.fitnessAt(newstring):
                self.updatePosition(newstring)
                self.search_list = []
                self.search_list_start = 0
                return newstring, self.fitnessAt(newstring)
                
        self.search_list_start = self.search_list_start + i + 1
        return self.string, self.fitness
        
    def fitnessAt(self,string):
        return self.landscape.fitness(string)
        
class Landscape:
    def __init__(self,n=10,parts=2,changes=0,ls_num=0,r=0,
                 unit=0,symmetric=False):
        self.n = n   
        self.ls_num = ls_num # Landscape number
        self.r = r # Parameter to vary the correlation of landscapes
        self.unit = unit
        self.network = self.active_indexes_modular(parts,changes,symmetric)
        
    def active_indexes_modular(self,parts,changes=0,symmetric=False):  
        # Same seed for units with the same landscape number
        random.seed(self.ls_num)
        elements_per_part = int(self.n/parts)
        b=np.zeros([self.n,self.n])
        c=np.zeros([self.n,self.n])
        for row in range(self.n):
            for col in range(self.n):
                if int(row/elements_per_part) == int(col/elements_per_part):
                    b[row,col] = 1
                    c[row,col] = 1+int(row/elements_per_part)  
                    
        for j in range(changes):
            while True:
                row1=random.randint(0,self.n-1)
                col1=random.randint(0,self.n-1)     
                row2=random.randint(0,self.n-1)
                col2=random.randint(0,self.n-1) 

                if (row1==col1 or row2==col2) or (row1==row2 and col1==col2):
                    continue
                if b[row1,col1]==0 or b[row2,col2]==1:
                    continue
                else:
                    break
    
            b[row1,col1] = 0            
            b[row2,col2] = 1 

            if symmetric:
                b[col1,row1] = 0            
                b[col2,row2] = 1 

        b=b.astype(int)
     
        activeindexes= []
        for location in range(0,self.n):
            activeindex=[]
            for i in range(0,self.n):
                if b[location][i]==1:
                    activeindex.append(i)
            activeindexes.append(activeindex)
        
        return activeindexes        
   
    def fitness(self,string):
        fitness = 0
        fitness1 = 0
        c = (self.r)**2 + (1-self.r)**2
        
        for location in range(0,self.n): 
            activestring=''
            for j in self.network[location]:
                activestring += string[j]
            
            # Create a unique seed for random number generator
            seed1= exp(1)*location + \
                   exp(2)*self.ls_num + \
                   int(activestring)
            random.seed(seed1)
            randomi1 = random.gauss(0,1/sqrt(c))
            
            seed2=exp(1)*location + \
                  exp(2)*self.ls_num + \
                  int(activestring) + \
                  exp(3)*(1+self.unit)
                  
            random.seed(seed2)
            randomi2 = random.gauss(0,1/sqrt(c))
            
            fitness += randomi1
            
            fitness1 += randomi2
            
    
        #print(fitness/self.n)
        #print(fitness1/self.n)
        return ((1-self.r)*fitness + self.r*fitness1)/self.n
        
        
    def isLocalOpt(self,string):
        for i in range(self.n):
            testi = list(string)
            if testi[i] == '0':
                testi[i] = '1'
            else:
                testi[i] = '0'
            testi = ''.join(testi)
            if self.fitness(string)<self.fitness(testi):
                return False
        return True

    def localOpts(self):
        n = self.n
        count = 0
        fitnesses = []
        for i in range(0,2**n):
            string=BitArray(uint=i, length=n).bin
            if self.isLocalOpt(string):
                count = count + 1
                fitnesses.append(self.fitness(string))
        return count,np.sort(fitnesses),np.mean(fitnesses)


    def allFitnesses(self):
        n = self.n
        fitnesses = []
        strings = []
        for i in range(0,2**n):
            string=BitArray(uint=i, length=n).bin
            strings.append(string)
            fitnesses.append(self.fitness(string))
        return strings,fitnesses
        
                
class Organization:
    def __init__(self,n,parts,changes,r,org_size,ls_num,z,strategy,symmetric):
        self.units = []
        self.localoptimumflag=False
        self.strategy=strategy
        self.r = r
        self.parts = parts
        initstrings = []
        
        for i in range(0,org_size):
            initstring = ''
            for j in range(0,n):
                random.seed(exp(1)*i+exp(2)*j+1000)
                initstring += str(random.randint(0,1))
            initstrings.append(initstring)
            l = Landscape(n,parts,changes,ls_num,r,i,symmetric)
            self.units.append(Unit(self,parts,initstring,l,z,i,strategy))

        self.updateOrgStats()
        self.form = -1000

    def updateOrgStats(self):
        ss = []
        ff = []
        ff_sub =[]
        
        for unit in self.units:
            ss.append(unit.string)
            ff.append(unit.fitness)
            
            if unit.num>0:
                ff_sub.append(unit.fitness)                
            
        self.form = len(set(ss)) # Different org forms
        self.fit = np.mean(ff) 
        self.fit_sub = np.mean(ff_sub)
        
    def orgFits(self):
        ff = []
        for unit in self.units:
            ff.append(unit.fitness)
        return ff
                         
    def findBest(self):
        bestfitness=0
        bestnumber = 0
        for unit in self.units:
            if unit.fitness > bestfitness:
                bestfitness = unit.fitness
                bestnumber = unit.num
        return self.units[bestnumber].string, bestfitness

    def substring_combinations(self):
        substringss = []
        for part in range(self.parts):
            substrings = []
            for unit in self.units:
                if unit.assigned_part == part:
                    substrings.append(unit.string[unit.low:unit.high])
            substringss.append(substrings)
                
        lista = list(product(*substringss))
        lista2 = []
        for i in range(len(lista)):
            lista2.append(''.join(lista[i])) # possible substring combinations
        
        return lista2
        
    def bestAdoptionVoluntary(self):
        beststring, bestfitness= self.findBest()                  
        
        for unit in self.units:
            if unit.fitness < unit.fitnessAt(beststring): 
                unit.updatePosition(beststring)
                
    def partitionVoluntaryAdoption(self):
        lista = self.substring_combinations()
        
        for unit in self.units:
            beststring = unit.string
            bestfitness = unit.fitness
            for i in range(len(lista)):
                string = lista[i]
                fitness = unit.fitnessAt(string)
                if fitness > bestfitness:
                    beststring = string
                    bestfitness = fitness
            unit.updatePosition(beststring)
                        
    def bestAdoptionVoluntaryFromAll(self):
        allstrings=[]   
        for unit in self.units:
            allstrings.append(unit.string)    
        
        for unit in self.units:
            bestfitness = unit.fitness
            beststring = unit.string
            for i in range(len(allstrings)):
                if bestfitness< unit.fitnessAt(allstrings[i]):
                    beststring = allstrings[i]
                    bestfitness = unit.fitnessAt(allstrings[i])
            unit.updatePosition(beststring)  
            
    def unitData(self):
        s = []
        f = []
        for unit in self.units:
            if not(unit.string in s):
                s.append(unit.string)
                f.append(unit.fitness)
        return s,f
              
    def sim(self):
        if self.strategy=='C':
            while not self.units[0].localoptimumflag:
                self.units[0].localsearch()
            for unit in self.units:
                unit.updatePosition(self.units[0].string)
        else:
            for unit in self.units:
                while not unit.localoptimumflag:
                    unit.localsearch()

        if self.strategy=='partition':
            self.partitionVoluntaryAdoption()
        elif self.strategy=='BAV':
            self.bestAdoptionVoluntary()
        elif self.strategy=='BAVall':
            self.bestAdoptionVoluntaryFromAll()

        self.updateOrgStats()




class Experiment:
    def __init__(self,PARTS,CHANGES,STRATEGY,LANDSCAPES,ORG_SIZE,z3=False,
                 symmetric=True,N=12):
                     
        cor = np.asarray([0,0.2,0.4,0.6,0.8,1])
        rr=(cor+np.sqrt(-(cor-1)*cor)-1)/(2*cor-1)
        #rr[5]=0.5
        self.rrange = rr                   
                     
        t0=time.time()        
        alpha = 0.05        
        self.n = N
        self.symmetric = symmetric
        if STRATEGY=='C':
            self.z = 2  # <----------------------
        else:
            self.z = 1
        
        if z3:
            self.z=3
        
        #w=np.asarray(RRANGE)
        #self.correlation = 1-w**2 / (w**2 + (1-w)**2)
        self.correlation = cor        
        self.strategy = STRATEGY
        self.landscapes = LANDSCAPES
        self.org_size = ORG_SIZE
        self.parts = PARTS
        self.changes = CHANGES   

        print(self.strategy + ' - ch:' + str(self.changes))
        
        self.fit = []
        self.fitconf = []
        self.orgs_r = [] 
        self.form = []
        self.formconf = []
        self.allfits = []
        self.allforms = []
        
        for r in self.rrange:
            orgs = [];
            fits = [];  # list of fitnesses on different landscapes
            fits_all = []
            forms = []
            for ls_num in range(self.landscapes):
                orgs.append(Organization(self.n,self.parts,self.changes,r,
                                         self.org_size,ls_num,self.z,
                                         self.strategy,
                                         self.symmetric))        
                
                orgs[ls_num].sim()
                  
                if self.strategy=='C':
                    fits.append(orgs[ls_num].fit_sub)
                    fits_all.append(orgs[ls_num].fit)
                else:
                    fits.append(orgs[ls_num].fit)
                    
                forms.append(orgs[ls_num].form)
                                        
            self.form.append(np.mean(forms))                            
            self.fit.append(np.mean(fits))
            self.allfits.append(fits)
            self.allforms.append(forms)            
            
            ci=scipy.stats.t.ppf(1-alpha/2,self.landscapes-1)             
            # ddof=1 -> sample stdev (keskiarvon keskivirhe)            
            self.fitconf.append(ci*np.std(fits,ddof=1)/sqrt(self.landscapes))
            self.formconf.append(ci*np.std(forms,ddof=1)/sqrt(self.landscapes)) 
            self.orgs_r.append(orgs)
                        
        t1=time.time()
        print('Duration:'+str(round((t1-t0)/60,2)))
        
if __name__ == '__main__':
    
    def baseline_sims(N,PARTS,ch,expname):
        #ORG_SIZE = 3 * 2**12
        ORG_SIZE = 2**12        
        landscapes = 50
        t0=time.time()        
        #orgs = []
        form = []
        fitmean = []
        fitstdev = []
        fitmin = []
        fitmax = []
        fitrange = []
        experiment = []

        for ls_num in range(landscapes):
            print(ls_num)
            o = Organization(N,PARTS,ch,0,ORG_SIZE,ls_num,1,'D',True)
            o.sim()
            allfits = o.orgFits()
            form.append(o.form)
            fitmean.append(o.fit)
            fitstdev.append(np.std(allfits,ddof=1))
            fitmin.append(np.min(allfits))
            fitmax.append(np.max(allfits))
            fitrange.append(np.max(allfits) - np.min(allfits))
            #experiment.append('N='+str(N)+' parts='+str(PARTS)+' ch='+str(ch))
            experiment.append(expname)            
                  
        t1=time.time()
        print('Duration:'+str(round((t1-t0)/60,2)))
        return pd.DataFrame(data=np.transpose([form,fitmean,fitstdev,fitmin,
                                               fitmax,fitrange,experiment]),
                                               columns=['form','mean','stdev',
                                               'min','max','range','experiment'])
    
    
    
    def plot_figure1():
        num=random.randint(1,1000)
        num=595
        print(num)
        
        ll1=[]
        ll2=[]
        wrange=[0,0.5,1]
        
        n=6
        ch=1
        parts=2
        for w in wrange:
            ll1.append(Landscape(n,parts,ch,num,w,0))
            ll2.append(Landscape(n,parts,ch,num,w,1))
                 
        #vrt. Ganco & Agarwal      
        xstrings=['000','001','011','010','110','111','101','100']
        ystrings=['000','001','011','010','110','111','101','100'] 
                
        n=len(xstrings[0])+len(ystrings[0])
        
        zz1=[]
        zz2=[]
        for w in wrange:
            zz1.append(np.zeros((2**int(n/2),2**int(n/2))))
            zz2.append(np.zeros((2**int(n/2),2**int(n/2))))
        
        for i in range(len(xstrings)):
            for j in range(len(ystrings)):
                for index in range(len(wrange)):
                    zz1[index][i,j]=ll1[index].fitness(xstrings[i]+ystrings[j])
                    zz2[index][i,j]=ll2[index].fitness(xstrings[i]+ystrings[j])
            
        x=range(2**int(n/2))
        y=range(2**int(n/2))
        
        index=0
        X, Y = np.meshgrid(x, y)
   
        fig=plt.figure()
        index=0 # index of w
         
        for w in wrange:        
            cmap1=cm.gray
            cmap2=cm.Greys

            
            #levels = np.linspace(0,1.2)
            levels = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3] 
            
            ax=plt.subplot(231+index);
            #C=plt.contourf(X,Y,zz1[index],levels,cmap=cmap1)
            C=plt.imshow(zz1[index],interpolation='none',vmin=0,vmax=1.4,
                         cmap=cm.Greys)

            plt.title('w='+str(w))
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            if index==0:
                plt.ylabel('Unit 1')
            #cbar = plt.colorbar(C)    

            for i in range(len(X)):
                for j in range(len(X)):
                    val = zz1[index][j][i]
                    #if val>0:
                    if(round(val,1)>=1):
                        cc='1'
                    else:
                        cc='0'
                    plt.text(i,j,ystrings[j]+'\n'+xstrings[i],size=6,color=cc,
                             va='center',ha='center')
                    

            #plt.xticks(y,ystrings,size=6)
            #plt.yticks(x,xstrings,size=6)  
            #ax.tick_params(direction='out', length=1, width=1)

            ax=plt.subplot(234+index);
            #C=plt.contourf(X,Y,zz2[index],levels,cmap=cmap2)
            C=plt.imshow(zz2[index],interpolation='nearest',vmin=0,vmax=1.4,
                         cmap=cm.Greys)      
            
            
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([]) 
            #cbar = plt.colorbar(C) 
            if index==0:
                plt.ylabel('Unit 2')
                
            for i in range(len(X)):
                for j in range(len(X)):
                    val = zz2[index][j][i]
                    #if val>0:
                    if(round(val,1)>=1):
                        cc='1'
                    else:
                        cc='0'
                    plt.text(i,j,ystrings[j]+'\n'+xstrings[i],size=6,color=cc,
                             va='center',ha='center')    
            
            #plt.xticks(y,ystrings,size=6)
            #plt.yticks(x,xstrings,size=6)            
            #ax.tick_params(direction='out', length=1, width=1)      
            
            index+=1
             
        
        plt.subplots_adjust(left=None, bottom=0.1, right=None, top=0.9,
                            wspace=0.1, hspace=0.1)
        fig.set_size_inches(7,4);
        plt.savefig('figure1.jpg',dpi=250)
        return zz1
                      
    def run_experiments():
        t00=time.time() 
        LS = 500;
        eee=dict()  
        eee['correlation'] = np.asarray([0,0.2,0.4,0.6,0.8,1])
                       
        ### FIGURE 4 
        units = 8               
        parts = 1           
        eee['C2 p1 ch0 u8'] = Experiment(parts,0,'C',LS,units).fit
        eee['C3 p1 ch0 u8'] = Experiment(parts,0,'C',LS,units,True).fit
        eee['BAV p1 ch0 u8'] = Experiment(parts,0,'BAV',LS,units).fit
        
        parts = 2           
        eee['C2 p2 ch0 u8'] = Experiment(parts,0,'C',LS,units).fit
        eee['C3 p2 ch0 u8'] = Experiment(parts,0,'C',LS,units,True).fit
        eee['BAV p2 ch0 u8'] = Experiment(parts,0,'BAV',LS,units).fit 
        eee['BAV p2 ch15 u8'] = Experiment(parts,15,'BAV',LS,units).fit        
                
        parts = 4           
        eee['C2 p4 ch0 u8'] = Experiment(parts,0,'C',LS,units).fit
        eee['C3 p4 ch0 u8'] = Experiment(parts,0,'C',LS,units,True).fit
        eee['BAV p4 ch0 u8'] = Experiment(parts,0,'BAV',LS,units).fit
        eee['BAV p4 ch15 u8'] = Experiment(parts,15,'BAV',LS,units).fit
        
        ### FIGURE 5
        
        units = 4
        
        parts = 1           
        eee['BAV p1 ch0 u4'] = Experiment(parts,0,'BAV',LS,units).fit
        eee['P p1 ch0 u4'] = Experiment(parts,0,'partition',LS,units).fit
                              
        parts = 2           
        eee['BAV p2 ch0 u4'] = Experiment(parts,0,'BAV',LS,units).fit
        eee['P p2 ch0 u4'] = Experiment(parts,0,'partition',LS,units).fit
        eee['BAV p2 ch15 u4'] = Experiment(parts,15,'BAV',LS,units).fit
        eee['P p2 ch15 u4'] = Experiment(parts,15,'partition',LS,units).fit                    
                
        parts = 4           
        eee['BAV p4 ch0 u4'] = Experiment(parts,0,'BAV',LS,units).fit
        eee['P p4 ch0 u4'] = Experiment(parts,0,'partition',LS,units).fit        
        eee['BAV p4 ch15 u4'] = Experiment(parts,15,'BAV',LS,units).fit
        eee['P p4 ch15 u4'] = Experiment(parts,15,'partition',LS,units).fit 
        
        units = 12
        
        parts = 1           
        eee['BAV p1 ch0 u12'] = Experiment(parts,0,'BAV',LS,units).fit
        eee['P p1 ch0 u12'] = Experiment(parts,0,'partition',LS,units).fit
                              
        parts = 2           
        eee['BAV p2 ch0 u12'] = Experiment(parts,0,'BAV',LS,units).fit
        eee['P p2 ch0 u12'] = Experiment(parts,0,'partition',LS,units).fit
        eee['BAV p2 ch15 u12'] = Experiment(parts,15,'BAV',LS,units).fit
        eee['P p2 ch15 u12'] = Experiment(parts,15,'partition',LS,units).fit                    
                
        parts = 4           
        eee['BAV p4 ch0 u12'] = Experiment(parts,0,'BAV',LS,units).fit
        eee['P p4 ch0 u12'] = Experiment(parts,0,'partition',LS,units).fit        
        eee['BAV p4 ch15 u12'] = Experiment(parts,15,'BAV',LS,units).fit
        eee['P p4 ch15 u12'] = Experiment(parts,15,'partition',LS,units).fit 
                
        
        print('Total duration:'+str(round((time.time()-t00)/60,2)))

        return pd.DataFrame.from_dict(eee)     
        

    def create_network(n,modules,seed,changes=0):
        random.seed(seed)
        elements_per_module = int(n/modules)
        b=np.zeros([n,n])
        c=np.zeros([n,n])
        for row in range(n):
            for col in range(n):
                if int(row/elements_per_module) == int(col/elements_per_module):
                    b[row,col] = 1
                    c[row,col] = 1+int(row/elements_per_module)                
                    
        for j in range(changes):
            while True:
                row1=random.randint(0,n-1)
                col1=random.randint(0,n-1)     
                row2=random.randint(0,n-1)
                col2=random.randint(0,n-1) 
                
                if (row1==col1 or row2==col2) or (row1==row2 and col1==col2) or b[row1,col1]==0 or b[row2,col2]==1:
                    continue
                else:
                    break
    
            b[row1,col1] = 0            
            b[row2,col2] = 1
    
            b[col1,row1] = 0            
            b[col2,row2] = 1 
                     
        return b
        
    def plotnetwork(sp,parts,ch,seed,otsikko,colors):
        N=12  
        b = create_network(N,parts,1000,ch)
        g = nx.DiGraph(b)
        gg = g.to_undirected()
        plt.subplot(sp)
        nx.draw_spring(gg,node_color=colors,node_size=50)
        plt.title(otsikko,fontsize=8)
          
    def plot_figure4(sp,D0,C2,D15,C3,otsikko):
        
        ax=plt.subplot(sp)
        plt.plot(cor,D0,'k.--',label='Decentralisation, decomposable (0 changes)')
        plt.plot(cor,D15,'k.-',label='Decentralisation, non-decomposable (15 changes)')
        plt.plot(cor,C2,'k--',label='Centralisation (search distance=2)')
        plt.plot(cor,C3,'k-',label='Centralisation (search distance=3)') 
        
        ax.fill_between(cor,D0,D15, facecolor='black', interpolate=True, alpha=0.15)
        ax.fill_between(cor,C2,C3, facecolor='black', interpolate=True, alpha=0.15)    
        
        plt.ylim([ylim_low,ylim_high])
        plt.xlim([xlim_low,xlim_high])
        plt.title(otsikko,fontsize=8)
        ax.set_xticks([xlim_low,xlim_high])                                                                                                 
        #ax.set_yticks([0,1])  
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left') 
        #plt.xlabel('Correlation',fontsize=xlabelsize)
        #ax.minorticks_on()
        ax.grid(True)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        #ax.grid(which='minor', linestyle=':', linewidth='0.5')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0+ box.height * 0.2, box.width,
                         box.height*0.8])
        
        if sp==133:
            leg=ax.legend(loc='upper center', bbox_to_anchor=(-0.7,-0.2),
                          ncol=2,framealpha=0,fontsize=8)
        if sp==131:
            plt.ylabel('Mean fitness (end of simulation)\n8 units',fontsize=8)
        if sp==132:
            plt.xlabel('Correlation',fontsize=xlabelsize)
        if sp==131:
            ax.spines['left'].set_visible(True)

    
    def plot_figure5(sp,D0,P0,D15,P15,otsikko):
        
        ax=plt.subplot(sp)
        plt.plot(cor,D0,'k.--',label='Decentralisation (decomposable, 0 changes)')
        plt.plot(cor,D15,'k.-',label='Decentralisation (non-decomposable, 15 changes)')
        plt.plot(cor,P0,'k--',label='Partition (decomposable, 0 changes)')
        plt.plot(cor,P15,'k-',label='Partition (non-decomposable, 15 changes)') 
        
        ax.fill_between(cor,D0,D15, facecolor='black', interpolate=True, alpha=0.15)
        ax.fill_between(cor,P0,P15, facecolor='black', interpolate=True, alpha=0.15)    
        
        plt.ylim([ylim_low,ylim_high])
        plt.xlim([xlim_low,xlim_high])
        plt.title(otsikko,fontsize=8)
        ax.set_xticks([xlim_low,xlim_high])                                                                                                 
        #ax.set_yticks([0,1])  
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left') 
        #plt.xlabel('Correlation',fontsize=xlabelsize)
        #ax.minorticks_on()
        ax.grid(True)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        #ax.grid(which='minor', linestyle=':', linewidth='0.5')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0+ box.height * 0.2, box.width,
                         box.height*0.8])
        
        if sp==236:
            leg=ax.legend(loc='upper center', bbox_to_anchor=(-0.75,-0.3),
                          ncol=2,framealpha=0,fontsize=8)
        if sp==231:
            plt.ylabel('Mean fitness (end of simulation)\n4 units',fontsize=8)
        if sp==234:
            plt.ylabel('Mean fitness (end of simulation)\n12 units',fontsize=8)
        if sp==235:
            plt.xlabel('Correlation',fontsize=xlabelsize)
        if sp==231 or sp==234:
            ax.spines['left'].set_visible(True)
            
        
    ####################################################
    ####################################################

    
     
         
    ''' 
    # Figure 1 : landscapes
    zz1 = plot_figure1() 
    '''    
    
    
    '''
    # Figure 2 : networks    
    fig=plt.figure()          
    seed = random.seed()
    plotnetwork(231,2,0,seed,'Decomposable\n0 changes\nK=5 (2 parts)',['k','k','k','k','k','k','w','w','w','w','w','w'])
    plotnetwork(232,2,1,seed,'Parly decomposable\n1 change\nK=5 (2 parts)',['k','k','k','k','k','k','w','w','w','w','w','w'])
    plotnetwork(233,2,15,seed,'Non-decomposable\n15 changes\nK=5 (2 parts)',['k','k','k','k','k','k','w','w','w','w','w','w'])
    plotnetwork(234,3,0,seed,'K=3 (3 parts)',['k','k','k','k','w','w','w','w','0.5','0.5','0.5','0.5'])
    plotnetwork(235,3,1,seed,'K=3 (3 parts)',['k','k','k','k','w','w','w','w','0.5','0.5','0.5','0.5'])
    plotnetwork(236,3,15,seed,'K=3 (3 parts)',['k','k','k','k','w','w','w','w','0.5','0.5','0.5','0.5'])
    fig.set_size_inches(6,5);
    plt.savefig('figure2.jpg',dpi=250)            
    '''


    
    # Figure 3:     
    #p2ch0 = baseline_sims(12,2,0,1)
    #p2ch0.to_csv('p2_ch0.txt',sep='\t')
    #p2ch15 = baseline_sims(12,2,15,2)
    #p2ch15.to_csv('p2_ch15.txt',sep='\t')
    #p3ch0 = baseline_sims(12,3,0,3)
    #p3ch0.to_csv('p3_ch0.txt',sep='\t')
    #p3ch15 = baseline_sims(12,3,15,4)
    #p3ch15.to_csv('p3_ch15.txt',sep='\t')
    
    p2ch0 = pd.DataFrame.from_csv('p2_ch0.txt',sep='\t')
    p2ch15 = pd.DataFrame.from_csv('p2_ch15.txt',sep='\t')
    p3ch0 = pd.DataFrame.from_csv('p3_ch0.txt',sep='\t')
    p3ch15 = pd.DataFrame.from_csv('p3_ch15.txt',sep='\t')
    
    
    ax=plt.subplot(121)
    plt.plot(p3ch0['form'],p3ch0['mean'],'kx',label='K=3',markersize=10,alpha=0.7)
    plt.plot(p2ch0['form'],p2ch0['mean'],'ko',label='K=5',markersize=10,alpha=0.5)
    plt.title('Decomposable\n(0 changes)')    
    plt.ylabel('Mean fitness (end of simulation)')
    plt.ylim([0.3,1.1])
    plt.xlim([0,160])
    plt.xlabel('Organizational forms\n(end of simulation)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.xticks([0,160])
    #ax.legend(loc=4,ncol=2,framealpha=0.5)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+ box.height * 0.2, box.width, box.height*0.8])
   
    ax=plt.subplot(122)
    plt.plot(p3ch15['form'],p3ch15['mean'],'kx',label='K=3',markersize=10,alpha=0.7)
    plt.plot(p2ch15['form'],p2ch15['mean'],'ko',label='K=5',markersize=10,alpha=0.5)
    plt.title('Non-decomposable\n(15 changes)')
    plt.ylim([0.3,1.1])
    plt.xlim([0,160])
    plt.xlabel('Organizational forms\n(end of simulation)')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.xticks([0,160])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+ box.height * 0.2, box.width, box.height*0.8])
    leg=ax.legend(loc='upper center', bbox_to_anchor=(0,-0.3), ncol=2,framealpha=0,fontsize=8)
    plt.savefig('figure3.jpg',dpi=250)     
        
        

    # Run simulation experiments or read simulation results from file
    #df = run_experiments()
    #df.to_csv('results.txt',sep='\t')
    #df = pd.DataFrame.from_csv('results.txt',sep='\t')
    

    '''
    # Figure 4 : Simulation results I
    ylim_high = 1
    ylim_low = 0
    xlim_high = 1
    xlim_low = 0
    xlabelsize=8
    cor = np.asarray([0,0.2,0.4,0.6,0.8,1])            
    fig = plt.figure()
    plot_figure4(131,df['BAV p4 ch0 u8'],df['C2 p4 ch0 u8'],df['BAV p4 ch15 u8'],df['C3 p4 ch0 u8'],'Easy\nK=2 (4 parts)')
    plot_figure4(132,df['BAV p2 ch0 u8'],df['C2 p2 ch0 u8'],df['BAV p2 ch15 u8'],df['C3 p2 ch0 u8'],'Intermediate\nK=5 (2 parts)')
    plot_figure4(133,df['BAV p1 ch0 u8'],df['C2 p1 ch0 u8'],df['BAV p1 ch0 u8'],df['C3 p1 ch0 u8'],'Difficult\nK=11 (1 part)')
    #fig = plt.gcf()
    fig.set_size_inches(6,4)
    plt.savefig('figure4.jpg', dpi=250)
    '''
    
    '''
    # Figure 5 : Simulation results II
    fig = plt.figure()
    plot_figure5(231,df['BAV p4 ch0 u4'],df['P p4 ch0 u4'],df['BAV p4 ch15 u4'],df['P p4 ch15 u4'],'Easy\nK=2 (4 parts)')
    plot_figure5(232,df['BAV p2 ch0 u4'],df['P p2 ch0 u4'],df['BAV p2 ch15 u4'],df['P p2 ch15 u4'],'Intermediate\nK=5 (2 parts)')
    plot_figure5(233,df['BAV p1 ch0 u4'],df['P p1 ch0 u4'],df['BAV p1 ch0 u4'],df['P p1 ch0 u4'],'Difficult\nK=11 (1 part)')
    plot_figure5(234,df['BAV p4 ch0 u12'],df['P p4 ch0 u12'],df['BAV p4 ch15 u12'],df['P p4 ch15 u12'],'Easy\nK=2 (4 parts)')
    plot_figure5(235,df['BAV p2 ch0 u12'],df['P p2 ch0 u12'],df['BAV p2 ch15 u12'],df['P p2 ch15 u12'],'Intermediate\nK=5 (2 parts)')
    plot_figure5(236,df['BAV p1 ch0 u12'],df['P p1 ch0 u12'],df['BAV p1 ch0 u12'],df['P p1 ch0 u12'],'Difficult\nK=11 (1 part)')      
    #fig = plt.gcf()
    fig.set_size_inches(6,5)
    plt.savefig('figure5.jpg', dpi=250)              
    '''