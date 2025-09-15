#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 18:53:37 2022

@author: lena
"""

### Librairies

import numpy as np
import matplotlib.pyplot as plt

### Parameters

dx = 1                       # Spatial interval 
T = 1000                     # Time at which simulation stops (it can stop earlier if one of the type disappear from the environment)
dt = 0.1 # np.round(m*dx**2/2,10)  # Time interval
conv_timing = "ger"         # Conversion timing : "ger" or "zyg"
r = 0.1                     # Intrasic growth rate
c = 0.9                     # Conversion rate
h = 0.4                     # Dominance
replicats = 100               # Number of runs for the simulation


# Carrying capacity for a spatial interval of size 1
K = 10**8
# K_range = [10**3, 10**5, 10**7]

# Disadvantage for drive
#s = 0.9                     
s_range = np.round(np.arange(0.4, 0.8, 0.1),3) #[0.55]
nb_l = len(s_range)

# Migration probability
#m = 0.2  
#sigma_range = np.arange(0.25, 4.25, 0.25)              
m_range = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2] #  np.round(sigma_range*((2*dt)/(dx**2)),3)  
#  [0.01, 0.02, 0.03, 0.04, 0.05]       #[0.1, 0.2, 0.3]
nb_c = len(m_range)

for s in s_range:
    for m in m_range : 
        print("\n----- s =", s, ", m =", m)
        title = f"s_{s}_m_{m}" # f"K_{int(np.log10(K))}_s_{s}"
        
        # Speed and space for the wave not to go outside the window    
        s_1 = c/(1-h*(1-c))   
        s_2 = c/(2*c*h + h*(1-c))
        lin = c*(1-2*s*h)-(1-c)*s*h  
        
        sigma = np.sqrt(np.round(m*((dx**2)/(2*dt)),3))
        v_cont = 2*sigma*np.sqrt(lin)
        nb_sites = (int(v_cont*T/dx*(1+1/5)//100)+2)*100
       
        for run in range(replicats):
            print(run)

            ### Initialization
        
            nD = np.zeros(nb_sites).astype(int); nD[:(nb_sites//5)] = K*dx
            nW = np.zeros(nb_sites).astype(int); nW[(nb_sites//5):] = K*dx
            fD = np.zeros(nb_sites)
            fW = np.zeros(nb_sites)
            chasing = 0
            nb_drive_around_last_wt = np.zeros(len(np.arange(300, T, dt)))
            i = 0
            last_wt = 5
            
            ### Evolution in time
            
            for t in np.arange(0, T, dt): 
                t = np.round(t,3)  
                
                # Mean number drive alleles around the last wild-type allele
                if t >= 300 : 
                    
                    last_wt = last_wt - 5 + np.where(nW[last_wt-5:]>0)[0][0]
                    nb_drive_around_last_wt[i] = nD[last_wt]
                    
                    if np.sum(nb_drive_around_last_wt[i-1:i+1]) == 0 : 
                        last_wt = last_wt + 20
                        last_wt = last_wt + np.where(nW[last_wt:]>0)[0][0]
                                       
                    i = i+1
                    
                    #plt.plot(nD, color = "crimson")
                    #plt.plot(nW, color = "cornflowerblue")
                    #plt.title(f"{nD[last_wt]}, t = {t}, vect= {nb_drive_around_last_wt[i-2:i+1]}")
                    #plt.xlim(last_wt-50, last_wt+50)
                    #plt.ylim(0, 150)
                    #plt.vlines(last_wt, 0, 150)
                    #plt.hlines(nD[last_wt], last_wt-50, last_wt+50)
                    #plt.show()
               
        
                ### Stop the simulation if the wave goes outside the window 
                
                if np.where(nD==max(nD))[0][0] > len(nD)-10 or t==T-dt : 
            
                    #print("t =",t)
                    #print ("chasing =", chasing)
                
                    # Save datas 
                    file = open(f"chasing/chasing_{title}.txt", "a") 
                    file.write(f"\n{chasing}")  
                    file.close() 
                    
                    
                    #print('\nMean drive number :', np.mean(nb_drive_around_last_wt))
                    
                    #if np.mean(nb_drive_around_last_wt) == 0 : 
                    #           print("!!!!!!!!!!!!!!!!PROBLEME!!!!!!!!!!!!!!")
                    #           nW_save = nW
                    #           nD_save = nD
                    
                    # Save datas 
                    file = open(f"chasing/nb_drive_{title}.txt", "a") 
                    file.write(f"\n{np.mean(nb_drive_around_last_wt)}")  
                    file.close() 

                    
                    #if run == 0 : 
                    #    fig, ax = plt.subplots()
                    #    bins = [x - 0.5 for x in range(0, nb_sites+1)]
                    #    plt.hist([np.arange(nb_sites), np.arange(nb_sites)], bins = bins, weights = [nD, nW], 
                    #                  histtype = 'barstacked', label = ['drive','wt'], color = ['crimson','cornflowerblue'])  
                    #    ax.set(xlabel='Space', ylabel='Number of alleles', ylim = [0,1.1*K*dx])
                    #    ax.set_title(f"Time {t}, K = 10^{int(np.log10(K))}, s = {s}")
                    #    plt.legend(loc="upper left")
                    #    plt.show()
                    
                    ### Plot graph
                    
                    
                    #plt.plot(nD_save, color = "crimson")
                    #plt.plot(nW_save, color = "cornflowerblue")
                    #plt.yscale('log')
                    #plt.xlim(1000, 1100)
                    #plt.show()
                    
        
                    break
            
                ### Chasing
                
                if np.where(nW>20)[0][0] < np.where(nD>0)[0][0]:
                    chasing = 1
    
     
                ### Birth and Death   
            
                # Index for empty and non empty sites
                extinct_index = np.where(nD+nW==0)[0]
                survive_index = np.delete(np.arange(nb_sites), extinct_index)
                # For non empty site, the fecundity is given by the following.
                sv_pop = nD[survive_index] + nW[survive_index]; sv_nD = nD[survive_index]; sv_nW = nW[survive_index]
                
                if conv_timing == "zyg" :
                    fD[survive_index] = ( 1 + r*(1-sv_pop/(K*dx)) ) * ( (1-s)*sv_nD + (1-s*h)*(1-c)*sv_nW +2*c*(1-s)*sv_nW ) /sv_pop
                    fW[survive_index] = ( 1 + r*(1-sv_pop/(K*dx)) ) * ( (1-c)*(1-s*h)*sv_nD + sv_nW ) /sv_pop           
                if conv_timing == "ger" : 
                    fD[survive_index] = ( 1 + r*(1-sv_pop/(K*dx)) ) * ( (1-s)*sv_nD + (1-s*h)*(1+c)*sv_nW ) /sv_pop
                    fW[survive_index] = ( 1 + r*(1-sv_pop/(K*dx)) ) * ( (1-c)*(1-s*h)*sv_nD + sv_nW ) /sv_pop
                    
                # For empty site, the fecundity is 0.
                fD[extinct_index] = 0
                fW[extinct_index] = 0
                # Check that all fecundity values are numbers.
                if len(np.argwhere(np.isnan(fD))) != 0 or len(np.argwhere(np.isnan(fW))) != 0 : 
                    print("Houston, we have a problem")
                # Add births, substract deaths (mortality = 1)
                nD = nD + np.random.poisson(fD*nD*dt) - np.random.poisson(nD*dt)            
                nW = nW + np.random.poisson(fW*nW*dt) - np.random.poisson(nW*dt)
                # Transform negative number of alleles into 0
                nD[np.where(nD<0)[0]]=0
                nW[np.where(nW<0)[0]]=0
                
                
                ### Migration  
                
                # Number of migrants in each site
                nD_mig = np.random.binomial(nD,m)
                nW_mig = np.random.binomial(nW,m)
                # Half migrate to the right, half to the left
                nD_mig_left = np.random.binomial(nD_mig,0.5); nD_mig_right = nD_mig - nD_mig_left
                nW_mig_left = np.random.binomial(nW_mig,0.5); nW_mig_right = nW_mig - nW_mig_left
                # Substract the migrants leaving
                nD -= nD_mig 
                nW -= nW_mig
                # ... except for those going outside the windows (they stay home)
                nD[0] += nD_mig_left[0]; nW[0] += nW_mig_left[0]
                nD[-1] += nD_mig_right[-1]; nW[-1] += nW_mig_right[-1]
                # Add the migrants in the neighboor sites
                nD[1:] += nD_mig_right[:-1]; nW[1:] += nW_mig_right[:-1] 
                nD[:-1] += nD_mig_left[1:]; nW[:-1] += nW_mig_left[1:]
                
                
chasing_matrix = np.ones((nb_l, nb_c))*(-1)
for i in range(nb_l):
    s = s_range[i]
    for j in range(nb_c): 
        m = m_range[j]        
        title = f"s_{s}_m_{m}"; print(title)
        chasing_matrix[i,j] = np.mean(np.loadtxt(f"chasing/chasing_{title}.txt"))
        
nb_drive_matrix = np.ones((nb_l, nb_c))*(-1)
for i in range(nb_l):
    s = s_range[i]
    for j in range(nb_c): 
        m = m_range[j]        
        title = f"s_{s}_m_{m}"; print(title)
        nb_drive_matrix[i,j] = np.mean(np.loadtxt(f"chasing/nb_drive_{title}.txt"))


fig, ax1 = plt.subplots(layout="constrained") 
ax2 = ax1.twinx() 
for i in range(nb_l):
    s = s_range[i]   
    ax1.plot(m_range, chasing_matrix[i,:], label=f"s = {s}", linestyle='--')
    ax2.plot(m_range, nb_drive_matrix[i,:], label=f"s = {s}")
ax2.legend(loc=0)
#ax1.set_yscale('log') 
ax1.set(xlabel = 'Migration rate', ylabel='Proportion of chasing') #, ylim = [1,100000])
ax2.set(ylabel='Nb drive', ylim = [0, 4000])
fig.savefig("chasing/chasing_nb_drive_s_m.png", format='png') 
plt.show()



