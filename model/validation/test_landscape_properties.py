# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 18:28:35 2017

@author: srsampsa
"""




def landscape_properties(parts,changes):
    landscape_mean = []
    landscape_stdev = []
    peak_count = []
    peak_mean = []
    n = 12
    unit = 0
    for ls_num in range(100):
        print(ls_num)
        l = Landscape(n,parts,changes,ls_num,0,0,True)
        s,f  = l.allFitnesses()
        landscape_mean.append(np.mean(f))
        landscape_stdev.append(np.std(f,ddof=1))
        peakcount, peakfits, peakmean = l.localOpts()
        peak_count.append(peakcount)
        peak_mean.append(peakmean)
    
    df = pd.DataFrame(np.transpose([landscape_mean,landscape_stdev, peak_count,peak_mean]),columns=['ls mean','ls stdev','peak count','peak mean'])  
        
    return df

'''    
p1ch0 = landscape_properties(1,0)
p1ch0.to_csv('landscape_p1ch0.txt',sep='\t')

p2ch0  = landscape_properties(2,0)
p2ch0.to_csv('landscape_p2ch0.txt',sep='\t')
p2ch15  = landscape_properties(2,15)
p2ch15.to_csv('landscape_p2ch15.txt',sep='\t')
p3ch0  = landscape_properties(3,0)
p3ch0.to_csv('landscape_p3ch0.txt',sep='\t')
p3ch15  = landscape_properties(3,15)
p3ch15.to_csv('landscape_p3ch15.txt',sep='\t')    
'''

#plt.plot(p1ch0['ls mean'],p1ch0['ls mean'])
#sns.pairplot(p1ch0,vars=['ls mean','ls stdev','peak count','peak mean'],kind='scatter',diag_kind='hist')
#sns.pairplot(p2ch0,vars=['ls mean','ls stdev','peak count','peak mean'],kind='scatter',diag_kind='hist')
#sns.pairplot(p3ch0,vars=['ls mean','ls stdev','peak count','peak mean'],kind='scatter',diag_kind='hist')
#sns.pairplot(p2ch15,vars=['ls mean','ls stdev','peak count','peak mean'],kind='scatter',diag_kind='hist')
#sns.pairplot(p3ch15,vars=['ls mean','ls stdev','peak count','peak mean'],kind='scatter',diag_kind='hist')


p1ch0 = pd.DataFrame.from_csv('landscape_p1ch0.txt',sep='\t')
p2ch0 = pd.DataFrame.from_csv('landscape_p2ch0.txt',sep='\t')
p2ch15 = pd.DataFrame.from_csv('landscape_p2ch15.txt',sep='\t')
p3ch0 = pd.DataFrame.from_csv('landscape_p3ch0.txt',sep='\t')
p3ch15 = pd.DataFrame.from_csv('landscape_p3ch15.txt',sep='\t')

#df = pd.DataFrame.from_csv('p2p3ch0ch15.txt',sep='\t')
var1='peak count'
var2='ls stdev'

plt.plot(p1ch0[var1],p1ch0[var2],'bx')
plt.plot(p2ch0[var1],p2ch0[var2],'kx')
plt.plot(p2ch15[var1],p2ch15[var2],'ko')
plt.plot(p3ch0[var1],p3ch0[var2],'rx')
plt.plot(p3ch15[var1],p3ch15[var2],'ro')
plt.xlabel(var1)
plt.ylabel(var2)

#peak count -- peak mean
#peak count -- ls mean


'''
var =  'ls mean'
ls_mean = []
ls_mean.append(np.mean(p1ch0[var]))
ls_mean.append(np.mean(p2ch0[var]))
ls_mean.append(np.mean(p2ch15[var]))
ls_mean.append(np.mean(p3ch0[var]))
ls_mean.append(np.mean(p3ch15[var]))

ls_mean_conf = []
landscapes = 100
ls_mean_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p1ch0[var],ddof=1)/sqrt(landscapes))
ls_mean_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p2ch0[var],ddof=1)/sqrt(landscapes))
ls_mean_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p2ch15[var],ddof=1)/sqrt(landscapes))
ls_mean_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p3ch0[var],ddof=1)/sqrt(landscapes))
ls_mean_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p3ch15[var],ddof=1)/sqrt(landscapes))


var =  'ls stdev'
ls_stdev = []
ls_stdev.append(np.mean(p1ch0[var]))
ls_stdev.append(np.mean(p2ch0[var]))
ls_stdev.append(np.mean(p2ch15[var]))
ls_stdev.append(np.mean(p3ch0[var]))
ls_stdev.append(np.mean(p3ch15[var]))

ls_stdev_conf = []
landscapes = 100
ls_stdev_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p1ch0[var],ddof=1)/sqrt(landscapes))
ls_stdev_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p2ch0[var],ddof=1)/sqrt(landscapes))
ls_stdev_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p2ch15[var],ddof=1)/sqrt(landscapes))
ls_stdev_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p3ch0[var],ddof=1)/sqrt(landscapes))
ls_stdev_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p3ch15[var],ddof=1)/sqrt(landscapes))

var =  'peak count'
peak_count = []
peak_count.append(np.mean(p1ch0[var]))
peak_count.append(np.mean(p2ch0[var]))
peak_count.append(np.mean(p2ch15[var]))
peak_count.append(np.mean(p3ch0[var]))
peak_count.append(np.mean(p3ch15[var]))

peak_count_conf = []
landscapes = 100
peak_count_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p1ch0[var],ddof=1)/sqrt(landscapes))
peak_count_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p2ch0[var],ddof=1)/sqrt(landscapes))
peak_count_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p2ch15[var],ddof=1)/sqrt(landscapes))
peak_count_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p3ch0[var],ddof=1)/sqrt(landscapes))
peak_count_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p3ch15[var],ddof=1)/sqrt(landscapes))

var =  'peak mean'
peak_mean = []
peak_mean.append(np.mean(p1ch0[var]))
peak_mean.append(np.mean(p2ch0[var]))
peak_mean.append(np.mean(p2ch15[var]))
peak_mean.append(np.mean(p3ch0[var]))
peak_mean.append(np.mean(p3ch15[var]))

peak_mean_conf = []
landscapes = 100
peak_mean_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p1ch0[var],ddof=1)/sqrt(landscapes))
peak_mean_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p2ch0[var],ddof=1)/sqrt(landscapes))
peak_mean_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p2ch15[var],ddof=1)/sqrt(landscapes))
peak_mean_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p3ch0[var],ddof=1)/sqrt(landscapes))
peak_mean_conf.append(scipy.stats.t.ppf(1-0.05/2,landscapes-1)*np.std(p3ch15[var],ddof=1)/sqrt(landscapes))


plt.subplot(221)
plt.errorbar([0,1,2,3,4],ls_mean,yerr=ls_mean_conf,label='ls mean')
plt.legend()
plt.subplot(222)
plt.errorbar([0,1,2,3,4],ls_stdev,yerr=ls_stdev_conf,label='ls stdev')
plt.legend()
plt.subplot(223)
plt.errorbar([0,1,2,3,4],peak_count,yerr=peak_count_conf,label='peak count')
plt.legend()
plt.subplot(224)
plt.errorbar([0,1,2,3,4],peak_mean,yerr=peak_mean_conf,label='peak mean')
plt.legend()
'''

#Mean over landscapes and within landscapes = 0
#stdev within landscapes 


# high K -> more peaks
#ch 0 -> more peaks 