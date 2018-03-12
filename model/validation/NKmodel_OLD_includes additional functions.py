# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 12:18:16 2015

"""
__author__ =  'Sampsa Ruutu'
__version__=  '0.425'

import random
from itertools import combinations, product
from math import exp, sqrt
import numpy as np
from bitstring import BitArray
import matplotlib.pyplot as plt
import scipy.stats
import time
#from scipy.stats.stats import pearsonr
import pickle
import pandas as pd

class Unit:
    '''
    An organizational unit. 
    '''
    def __init__(self,org,parts,string,landscape,z,num,strategy,maxsearches,mut=False):
        '''
        Initializes the organizationa unit. If the strategy is set to
        'partition', the unit is assigned to search a specific partition.
        '''
        self.num = num # Number to identify a unit within the organization
        self.string = string # Determines the location of the landscape 
        self.landscape = landscape # Landscape instance
        self.z = z # Maximum search distance
        self.strategy = strategy # Determines how the unit performs search
        self.org = org # Organization that the unit is part of
        self.n = len(self.string) # Length of the search string
        self.localoptimumflag=False # Has the local optimum been reached?
        self.globaloptimumflag=False # Has the global optimum been reached?
        self.mut = mut # Is mutation on?
        self.search_list = []
        self.search_list_start = 0
        self.maxsearches = maxsearches

        # Additional variables currently not in use        
        #self.strings = [] # List of past strings
        #self.fitnesses = []  # List of past fitnesses
        #self.incrementallyadoptedflag = False

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
        '''        
        Updates the string and fitness. 
        '''
        self.string =string
        self.fitness = self.fitnessAt(string)
        
        #self.strings.append(self.string)
        #self.fitnesses.append(self.fitness)
                
    def imitatePosition(self,imitatestring):
        '''
        Imitates the string of another unit. Used in 'Best Adoption'. If
        mutation is on, the imitation differs by one bit.
        '''        
        if self.mut:        
            random.seed()
            i=random.randint(0,len(imitatestring)-1) 
            mutated = list(imitatestring)
            if imitatestring[i]=='0':
                mutated[i]='1'
            else:
                mutated[i]='0'
            mutatedstring=''.join(mutated)
            if self.fitness < self.fitnessAt(mutatedstring):
                self.updatePosition(mutatedstring)
        else:
            self.updatePosition(imitatestring)
                
    def localsearch(self):  
        '''
        The number of bits to be changed in local search is determined by
        the search length z. Assumption that the unit is able to identify a
        better string in the neighbourhood but is unable to evaluate the
        maximum gradient.
        '''            
            
        # Search is not performed if the unit is already at a local optimum
        if self.localoptimumflag:
            return self.string, self.fitness
                
        # Create a list of possible combinations of string indexes to change
        
        if(self.search_list == []):
            for L in range(1, self.z+1):
                for subset in combinations(range(self.low,self.high), L):
                    self.search_list.append(subset)
                
            # Randomize the list
            random.seed()
            random.shuffle(self.search_list)
        
        # Substrings based on the
        str_begin = self.string[0:self.low]
        str_end = self.string[self.high:]
        offset = self.low    
        
        # Go through the items on the search list. 
        #maxsearches = len(self.search_list)
        
        for i in range(self.maxsearches):
            
            if i+self.search_list_start == len(self.search_list):
                #print('local optimum')
                self.localoptimumflag = True
                return self.string, self.fitness

            indexes=self.search_list[i+self.search_list_start]
            
            localsearchstring = list(self.string[self.low:self.high])
                        
            #print('start='+str(self.search_list_start)+', i='+str(i)+ ' => ' + str(i+self.search_list_start))            
            
            # Change the bits of the active indexes (1->0 and 0->1)                        
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
                #print('found better')
                return newstring, self.fitnessAt(newstring)
                
        self.search_list_start = self.search_list_start + i + 1
        
        #print('no change, i='+str(i))
        return self.string, self.fitness
        
    def fitnessAt(self,string):
        '''
        Returns the fitness at a given point on the landscape.
        '''
        return self.landscape.fitness(string)
        
    def toGlobalMax(self):
        if not(self.globaloptimumflag):
            beststring,bestfitness = self.landscape.globalmax(self.n)
            self.updatePosition(beststring)
            self.globaloptimumflag=True

    
class Landscape:
    ''' Landscape is the fitness landscape of an organizational unit. '''
    def __init__(self,network_type,n=10,k=4,parts=2,changes=0,ls_num=0,r=0,
                 unit=0,symmetric=False):
        ''' Landscape initialization '''
        self.n = n   
        self.k=k
        self.ls_num = ls_num # Landscape number
        self.r = r # Parameter to vary the correlation of landscapes
        self.unit = unit
        
        if network_type=='modular':
            self.network = self.active_indexes_modular(parts,changes,symmetric)
            
        elif network_type=='nearest': 
            self.network = self.active_indexes_nearest()
        else:
            self.network = self.active_indexes_file()            
        
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

        
    def active_indexes_modular(self,parts,changes=0,symmetric=False):
        ''' Creation of a network with  '''
        
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
        #print(b)
       
        activeindexes= []
        for location in range(0,self.n):
            activeindex=[]
            for i in range(0,self.n):
                if b[location][i]==1:
                    activeindex.append(i)
            activeindexes.append(activeindex)
        
        #print(activeindexes)
        return activeindexes        
        
   
    def active_indexes_nearest(self): 
        '''
        Active indexes -> locus itself and the following k locuses
        (with wrap around)
        '''
        activeindexes = []
        
        for location in range(0,self.n):
            # List of indexes that contribute to the locus fitness    
            activeindex=[] 
            for j in range(0,self.k+1):
                activeindex.append((location+j)%self.n)
            activeindexes.append(activeindex)
            
        return activeindexes
        
    def active_indexes_file(self): 
        activeindexes = []
        txt = open('lok.txt')
        str1=txt.read().split('\n')
        txt.close()
        for location in range(0,self.n):
            activeindex=[]
            for i in range(0,self.n):
                if str1[location][i]=='1':
                    activeindex.append(i)
            activeindexes.append(activeindex)
        
        return activeindexes
               
   
   
    def fitness(self,string):
        #n= len(string)
        fitness = 0
        fitness1 = 0
        
        c = (self.r)**2 + (1-self.r)**2
        
        # käydään läpi jokainen bitti stringissa
        for location in range(0,self.n): 
        
            # muodostetaan stringi perustuen aktiivisiin indekseihin   
            activestring=''
            #print(self.network[location])
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
    
        return ((1-self.r)*fitness + self.r*fitness1)/self.n
                
    def globalmax(self):
        n = self.n
        bestfitness = 0
        beststring = ''
        for i in range(0,2**n):
            string=BitArray(uint=i, length=n).bin
            fitness=self.fitness(string)
            if fitness>bestfitness:
                bestfitness=fitness
                beststring=string
        return beststring,bestfitness
        
    def globalmin(self):
        n = self.n
        bestfitness = 0
        beststring = ''
        for i in range(0,2**n):
            string=BitArray(uint=i, length=n).bin
            fitness=self.fitness(string)
            if fitness<bestfitness:
                bestfitness=fitness
                beststring=string
        return beststring,bestfitness
        
    def allFitnesses(self):
        n = self.n
        fitnesses = []
        strings = []
        for i in range(0,2**n):
            string=BitArray(uint=i, length=n).bin
            strings.append(string)
            fitnesses.append(self.fitness(string))
        return strings,fitnesses
        
    def fitnesses(self,n):
        fitnesses = []
        for i in range(0,2**n):
            string=BitArray(uint=i, length=n).bin
            fitnesses.append(self.fitness(string))
        return fitnesses
        
        
class Organization:
    def __init__(self,network_type,n,k,parts,changes,r,org_size,ls_num,z,
                 randomstart,strategy,mut,budget,symmetric,maxsearches):
        #print('new org')
        self.units = [] # Organization koostuu listasta yksikköjä
        self.localoptimumflag=False
        self.incrementaladoptionflag = False
        self.strategy=strategy
        self.r = r
        self.parts = parts
        self.budget = budget
        
        initstrings = []
        
        for i in range(0,org_size):
            initstring = ''
            for j in range(0,n):
                if randomstart: # satunnaisesti valittu aloituspiste
                    random.seed(exp(1)*i+exp(2)*j+1000)
                    initstring += str(random.randint(0,1))
                else: # aloitus pisteestä '0000..'
                    initstring += '0'
            initstrings.append(initstring)
            l = Landscape(network_type,n,k,parts,changes,ls_num,r,i,symmetric)
            self.units.append(Unit(self,parts,initstring,l,z,i,strategy,maxsearches,mut))
            
        #self.forms = []
        #self.fits = []
        #self.fitsmax = []
        #self.fits_sub = []
        
        self.updateOrgStats()
        
        self.form = -1000

    
    def updateOrgStats(self):
        ss = []
        ff = []
        
        #ss_sub =[]
        ff_sub =[]
        
        for unit in self.units:
            ss.append(unit.string)
            ff.append(unit.fitness)
            
            if unit.num>0:
                #ss_sub.append(unit.string)
                ff_sub.append(unit.fitness)                
            
        self.form = len(set(ss)) # Different org forms
        #self.forms.append(self.form) # history of org forms
        
        self.fit = np.mean(ff) 
        #self.fits.append(self.fit) # list of avg fitnesses
        #self.fitmax = np.max(ff)
        #self.fitsmax.append(self.fitmax)
        
        self.fit_sub = np.mean(ff_sub)
        #self.fits_sub.append(self.fit_sub)
                         
        
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
                unit.imitatePosition(beststring)
                
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
            unit.imitatePosition(beststring)  
            
        
    def bestAdoptionEnforced(self):
        beststring, bestfitness= self.findBest()                  
        
        for unit in self.units:
            if unit.fitness < bestfitness:
                unit.imitatePosition(beststring)
                
    def unitData(self):
        s = []
        f = []
        for unit in self.units:
            if not(unit.string in s):
                s.append(unit.string)
                f.append(unit.fitness)
        return s,f
              
    def sim_vanha(self): # Organisaation simulointi
        if self.strategy=='C':
            while self.budget>0:
                self.units[0].localsearch()
                self.budget-=1
            for unit in self.units:
                unit.updatePosition(self.units[0].string)
        else:
            while self.budget>=0:
                for unit in self.units:
                    if self.budget>=0:
                        unit.localsearch() 
                        self.budget -=1
                    else:
                        break 

        
        if self.strategy=='partition':
            self.partitionVoluntaryAdoption()
        elif self.strategy=='BAV':
            self.bestAdoptionVoluntary()
        elif self.strategy=='BAVall':
            self.bestAdoptionVoluntaryFromAll()
        elif self.strategy=='BAE':
            self.bestAdoptionEnforced()    
        

        self.updateOrgStats()
        
    def sim(self): # Organisaation simulointi
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
        elif self.strategy=='BAE':
            self.bestAdoptionEnforced()    
        

        self.updateOrgStats()

  


  
        
