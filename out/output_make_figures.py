#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:54:30 2023
@author: Marine Laporte - Université Claude Bernard Lyon 1
@last-update: 24/10/2024
@contact: marinelaporte.ensg@gmail.com
"""

## Libraries we use 
## ------------------
import sys 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.pyplot import cm
import plotly.express as px
import numpy as np
import math
from math import exp 
from mpmath import mp
from datetime import date 
import scipy.signal as ss
import scipy.stats as st
## ------------------


# =============================================================================
# # MAIN FUNCTIONS
# =============================================================================
def read_outfiles_prop(infile):
    """Read the outfile prop.txt
    only one column"""
    try:
        fid = open(infile, "r")
    except:
        print( "Error : %s does not exist !" % infile)
        sys.exit()
    else:
        lines = fid.readlines()
        fid.close()
    prop = np.zeros(len(lines))
    flag=0
    for line in lines :
        prop[flag] = float(line.split()[0])
        flag+=1
    return prop


def read_dataset(infile):
    try:
        fid = open(infile, "r")
    except:
        print( "Error : %s does not exist !" % infile)
        sys.exit()
    else:
        lines = fid.readlines()
        fid.close() 
    m=[0]*len(lines)
    t=[0]*len(lines)
    flag=0
    for line in lines :
        t[flag] = float(line.split()[0])
        m[flag] = float(line.split()[1])
        flag+=1
    return t, m

def read_dataset_disc(infile,disc_calc, num):
    try:
        fid = open(infile, "r")
    except:
        print( "Error : %s does not exist !" % infile)
        sys.exit()
    else:
        lines = fid.readlines()
        fid.close() 
    m=[]
    t=[]
    flag=0
    for line in lines :
        time= float(line.split()[0])
        
        if num == 1 :
            disc= disc_calc[0]
            if time<=disc :
                t.append( time)
                m.append(float(line.split()[1]))
                flag+=1
        elif num == len(disc_calc)+1:
            disc= disc_calc[len(disc_calc)-1]
            if time>=disc :
                t.append(time)
                m.append( float(line.split()[1]))
                flag+=1
        else :
            disc= disc_calc[num-1]
            if time > disc_calc[num-2] and time <= disc :
                t.append(time)
                m.append( float(line.split()[1]))
                flag+=1
    print('There is %s events in the %s segment'%(str(flag),str(num)))
    return t, m

def adjR(x, y, degree):
    results = {}
    coeffs = np.polyfit(x, y, degree)
    p = np.poly1d(coeffs)
    yhat = p(x)
    ybar = np.sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((y - ybar)**2)
    results['r_squared'] = 1- (((1-(ssreg/sstot))*(len(y)-1))/(len(y)-degree-1))

    return results



def read_outfiles_like(infile):
    try:
        fid = open(infile, "r")
    except:
        print( "Error : %s does not exist !" % infile)
        sys.exit()
    else:
        lines = fid.readlines()
        fid.close()
    mat=np.zeros((len(lines), len(lines[0].split())))
    i=0
    for line in lines :
        line = np.array(line.split())
        line = [float( a) for a in line]
        for j  in range(len(line)):
            mat[i,j]=line[j]
            
        i+=1
    #mat = np.flip(mat, axis=0)
    return mat

def read_outfiles_acc(infile):
    try:
        fid = open(infile, "r")
    except:
        print( "Error : %s does not exist !" % infile)
        sys.exit()
    else:
        lines = fid.readlines()
        fid.close()
        
    N_acc = [int(a) for a in lines[0].split()]
    N_tot =  [int(b) for b in lines[1].split()]
    Pct = [float(c) for c in lines[2].split()]


    return N_acc, N_tot, Pct

def read_outfile_hypercube(infile):
    try:
        fid = open(infile, "r")
    except:
        print( "Error : %s does not exist !" % infile)
        sys.exit()
    else:
        lines = fid.readlines()
        fid.close()
        
    Mu = [float(a) for a in lines[0].split()]
    B =  [float(b) for b in lines[1].split()]
    S = [float(c) for c in lines[2].split()]
    I = [int(round(float(d)) )for d in lines[3].split()]
    return Mu, B, S, I


def figure_fit(Dataset, disc_calc, seg_i,  b_mle, m_mle,s_mle):
        times_di, mags_di = read_dataset_disc(Dataset, disc_calc, seg_i)
        beta_mle = b_mle * np.log(10)
        bins=np.linspace(0,6,31) # pourquoi jusqu'à 15 ?
        n,x=np.histogram(mags_di,bins=bins,density=True)
        print(len(x))
        print(len(n))
        titre_segment='Segment'+str(seg_i)
        n=n/sum(n)
        # COMPUTE BEST FITTING DISTRIBUTION 
        # detection function that maximizes the likelihood
        Q_max=np.zeros(30)
        E_max=np.zeros(30)
        
        for j in range(30):
            Q_max[j]= exp(-exp(-(x[j+1]-m_mle)/s_mle))* beta_mle * exp(-beta_mle*x[j+1])
            E_max[j]=(0.5+(math.erf((x[j+1]-m_mle)/(math.sqrt(2)*s_mle))/2)) *beta_mle * exp(-beta_mle*x[j+1])
        Q_max = Q_max/np.sum(Q_max)
        E_max = E_max/np.sum(E_max)
        
        Q_chi2=0
        E_chi2=0
        for j in range(30):

            if n[j]!=0:
                Q_chi2 = Q_chi2 + mp.mpf((Q_max[j]-n[j])**2/float(n[j]))
                E_chi2 = E_chi2 + mp.mpf((E_max[j]-n[j])**2/float(n[j]))
        print('CHI2_Gumbel :',Q_chi2)
        print('CHI2_Erf :', E_chi2)
        Q_chi2 = round(Q_chi2,2)
        E_chi2 = round(E_chi2,2)
        test = st.chisquare(n, f_exp=Q_max, ddof=25)

        print(test)
        # print(adjR(x[1:],n,2))
        # print(adjR(x[1:],n,3))
        # print(adjR(x[1:],n,4))
        # print(adjR(x[1:],n,5))
        # print( np.polyfit(x[1:],n,4))
        fig3, axs3 = plt.subplots(2,2, figsize=(10,10))
        plt.subplots_adjust(left=0.1,
                        bottom=0.2,
                        right=0.9,
                        top=0.9,
                        wspace=0.4,
                        hspace=0.3)
        # subplot1 : comparing mag distributions data vs model
        axs3[0,0].hist(mags_di,bins=bins, edgecolor='dimgrey', color='steelblue', alpha=0.5)
        
        axs3[0,0].set_xlim(0,6)
        axs3[0,0].set_xlabel('Magnitudes')
        axs3[0,0].set_ylabel('Nb of events (catalog)')
        axs30 = axs3[0,0].twinx()
        axs30.scatter(x[1:],n, c='grey',edgecolor='black',alpha=0.8,s=8)
        axs30.plot(x[1:],n, c='black',label='data')
        axs30.plot(x[1:], Q_max,c='red',label='Gumbel : Chi2:'+str(Q_chi2))
        axs30.plot(x[1:], E_max,linestyle='--', c='orange',label='Erf : Chi2:'+str(E_chi2))
        axs30.set_ylabel('Frequency ')
        axs30.legend()
        axs30.set_ylim(0,max(n)+0.008)
        axs3[0,0].grid(lw=0.4)
       
        # subplot2 : comparing Log(N) data vs model
        vec = np.arange(0,np.round(np.max(mags)),0.1)
        d = np.zeros(len(vec))
        for k in range(len(vec)):
         	ind = np.where(mags>vec[k])
         	d[k] = len(ind[0])/len(mags)
            
        axs3[0,1].scatter(x[1:],np.log10(n), c='grey', edgecolor='black',alpha=0.8,s=40,label='data')
        axs3[0,1].set_xlim(0,6)
        axs3[0,1].set_xlabel('Magnitudes')
        axs3[0,1].set_ylabel('Log(N)')
        #axs31 = axs3[0,1].twinx()
        axs3[0,1].plot(x[1:], np.log10(Q_max),c='red',label='Gumbel : Chi2:'+str(Q_chi2))
        axs3[0,1].plot(x[1:], np.log10(E_max),ls='--', c='orange',label='Erf : Chi2:'+str(E_chi2))
        axs3[0,1].set_ylabel('Log(F) ')
        axs3[0,1].legend()
        axs3[0,1].grid(lw=0.4)
       
        # subplot3 : comparing CNE data vs model
        cn=np.ones(30)
        cqmax=np.ones(30)
        cemax=np.ones(30)
        for i in range(1,30):
            cn[i]=sum(n[i:30])
            cqmax[i]=sum(Q_max[i:50])
            cemax[i]=sum(E_max[i:50])
        
        axs3[1,0].scatter(x[1:],cn, c='grey', edgecolor='black',alpha=0.8,s=40,label='data')
        axs3[1,0].set_xlim(0,6)
        axs3[1,0].set_xlabel('Magnitudes')
        axs3[1,0].set_ylabel('Cumulative number of magnitudes')
        #axs31 = axs3[1,0].twinx()
        axs3[1,0].plot(x[1:], cqmax,c='red',label='Gumbel : Chi2:'+str(Q_chi2))
        axs3[1,0].plot(x[1:], cemax,ls='--', c='orange',label='Erf : Chi2:'+str(E_chi2))

        axs3[1,0].legend()
        axs3[1,0].grid(lw=0.4) 
        
        axs3[1,1].scatter(x[1:],np.log10(cn), c='grey', edgecolor='black',alpha=0.8,s=40,label='data')
        axs3[1,1].set_xlim(0,6)
        axs3[1,1].set_xlabel('Magnitudes')
        axs3[1,1].set_ylabel('Log(Cumulative number)')
        axs3[1,1].plot(x[1:], np.log10(cqmax),c='red',label='Gumbel : Chi2:'+str(Q_chi2))
        axs3[1,1].plot(x[1:], np.log10(cemax),ls='--', c='orange',label='Erf : Chi2:'+str(E_chi2))
        axs3[1,1].legend()
        axs3[1,1].grid(lw=0.4) 
        fig3.suptitle('Fitting the Magnitude distribution for segment '+str(seg_i)+'\n'+'$Detection$ + $Richter$'+ '  :  ['+str(m_mle)+'-'+str(s_mle)+']'+' - ['+str(b_mle)+']' )
        fig3.savefig('Fit_'+titre_segment+'.pdf')
        plt.show()

def Bpos_function(mags, w):
    n_events =len(mags)
    difference = np.zeros(n_events-1)
    beta_positive= np.zeros(n_events-w)
    #difference[0] = 0
    
    for i in range (0, n_events-1):
        difference [i] = abs(mags[i+1] - mags[i])
    flag=0
    for m in range(len(difference)-(w-1)):
       flag+=1
       windowed_mean = sum(np.ones(w) * difference[m:m+w]) / w 
       Minimum_difference = min(difference[m:m+w])
       windowed_beta = 1/float(windowed_mean-Minimum_difference)
       beta_positive[m]=windowed_beta/np.log(10)
       #print(flag)
    return beta_positive






if __name__ == "__main__" :


# =============================================================================
#   Default outfiles names in OUT/
# =============================================================================
    Prop_out = 'prop.txt'
    Nc_out =   'NC_tot.txt'
    Dataset =  'input.txt'
    Mc_out =   'Likelihood_Mu.txt'
    b_out =    'Likelihood_b.txt'
    S_out =   'Likelihood_Sigma.txt'
    Acc_out =  'Acceptances.txt'
    Espace_out = 'hypercube.txt'
   
    synthetics='False'
# =============================================================================
#   If run on synthetics - give synthetic parameters to retrieve
# =============================================================================
    # T_vec = [0,20, 40, 55, 85,100]
    # N_vec = [2000,800,1500,500,2200]
    # b_vec =  np.array([0.7, 0.7, 0.7, 0.7, 0.7])
    
    
    # Mc_vec = np.array([2, 2.4, 1.3, 3, 4])
    # mu_vec = np.array([0.5, 1.3, 0.5, 1.3, 0.5])
    # sigma_vec = np.array([0.32,0.32,0.32,0.5,0.32])
    
    # T_vec = [0, 15, 20, 40, 60, 70, 88, 100]
    # N_vec = [1500,500,2000,200,2000,1500,500]




    # b_vec =  np.array([0.77, 0.75, 0.9, 0.55, 0.7,0.7,0.9])*np.log(10)
    # mu_vec = np.array([0.5, 0.8, 0.6, 0.75, 0.5,0.6,0.6])
    # sigma_vec = np.array([0.32,0.32,0.32,0.32,0.32,0.32,0.32])
    # # T_vec = [0,   20, 40,60, 65,100]
    # # N_vec = [1400,500,130,1000,1500]
    
    # # beta_vec =  np.array([0.77*np.log(10), 0.75*np.log(10), 0.9*np.log(10), 0.55*np.log(10), 0.7*np.log(10)])
    # # mu_vec = np.array([0.5, 0.8, 0.9, 0.75, 0.5])
    # # sigma_vec = np.array([0.32,0.32,0.32,0.32,0.32])
    
    # # 15.05
    # # T_vec = [0, 15, 20, 40, 60, 70, 88, 100]
    # # N_vec = [1500,500,2000,200,2000,1500,500]
    
    # # 31.08
    # T_vec = [0, 15, 20, 40, 60, 70, 88, 100]
    # N_vec = [750,250,1000,1000,500,900,600]  
    # b_vec =  np.array([0.77,1.04,1.04,0.8,0.8,0.75,0.7])*np.log(10)
    # mu_vec = np.array([0.2,0.26,0.26,0.15,0.15,0.1,0.12])
    # sigma_vec = np.array([0.32,0.32,0.24,0.32,0.32,0.32,0.32])
    
    
    # b_vec =  np.array([0.77, 0.75, 0.9, 0.55, 0.7,0.7,0.9])*np.log(10)
    # mu_vec = np.array([0.5, 0.8, 0.6, 0.75, 0.5,0.6,0.6])
    # sigma_vec = np.array([0.32,0.32,0.32,0.32,0.32,0.32,0.32])
    # name= 'TB_6DISC_0509_GUMB_5000_PY'
    # T_vec = [0, 15, 20, 40, 60, 70, 88, 100]
    # N_vec = [750,250,1000,1000,500,900,600]  
    # b_vec =  np.array([0.77, 1.04, 0.90, 0.8, 0.8, 0.75, 0.7])*np.log(10)
    # mu_vec = np.array([0.20, 0.26, 0.15, 0.1, 0.3, 0.10, 0.12])
    # sigma_vec = np.array([0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32])
    
    # name= 'TB_6DISC_2109_GUMB_5000_PY_BPOS'
    # T_vec = [0, 15, 20, 40, 60, 70, 88, 100]
    # N_vec = [750,250,1000,1000,500,900,600]  
    # b_vec =  np.array([0.77, 1.4, 0.95, 0.6, 1.2, 0.75, 0.7])*np.log(10)
    # mu_vec = np.array([0.20, 0.9, 0.15, 0.1, 0.8, 0.10, 0.6])
    # sigma_vec = np.array([0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32])
    
    
    # name= 'TB_6DISC_CLEAR_2909_EXP_5000_PY'
    # T_vec = [0, 15, 20, 40, 60, 70, 88, 100]
    # N_vec = [750,250,1000,1000,500,900,600]  
    # b_vec =  np.array([0.77, 1.4, 0.95, 0.6, 1.2, 0.75, 0.7])*np.log(10)
    # mu_vec = np.array([1.5, 1.0, -0.05, 0.8, -0.8, 1.3, -2])
    # sigma_vec = np.array([0.14, 0.14, 0.3, 0.14,0.14,0.14, 0.14])
    
    
    name= 'FIGURE_STAI_250_evts_mucst_v3'

    year0 = 0
    #name ='NorthChile_0111'
    #name= 'TB_6DISC_UNCLEAR_0611_Mmove_x4'
    # T_vec = [0, 15, 20, 40, 60, 70, 88, 100]
    # N_vec = [750,150,1000,1000,500,900,600]  
    # b_vec =  np.array([0.77, 0.8, 0.8, 0.9, 1.0, 0.7, 0.7])*np.log(10)
    # mu_vec = np.array([1.5, 1.0,1.1, 0.2, 0.4, 0.3, 0.6])
    # sigma_vec = np.array([0.14, 0.14, 0.3, 0.14,0.14,0.14, 0.14])
    # Ntot = 4384
    # T_vec = [0, 10, 20, 40, 60, 70, 88, 100]
    # N_vec = [348,597,109,1961,714,503,152]  
    # b_vec =  np.array([0.7,1.1, 1.0, 1.0, 0.8, 0.7, 0.65])*np.log(10)
    # mu_vec = np.array([0.5, 0.5,1.5, 0.2, 0.6, 0.9, 0.9])
    # sigma_vec = np.array([0.14, 0.14, 0.3, 0.14,0.3,0.4, 0.14])
    
    
    T_vec = [0, 10, 20, 40, 60, 70, 88, 100]
    N_vec = [514,647,261,2538,727,672,324]  
    Ntot = sum(N_vec)
    b_vec =  np.array([0.7,1.1, 1.0, 1.0, 0.8, 0.7, 0.65])*np.log(10)
    mu_vec = np.array([0.8, 0.8,1.5, 0.5, 0.75, 0.9, 0.9])
    sigma_vec = np.array([0.35, 0.3, 0.2, 0.15,0.3,0.4, 0.15])
    
    
# =============================================================================
#     Reading outfiles  
# =============================================================================
    Prop_vec = read_outfiles_prop(Prop_out)
    Nc_vecs =read_outfiles_like(Nc_out)
    times, mags = read_dataset(Dataset)
    times2 = [a+year0 for a in times]
    Like_mu = read_outfiles_like(Mc_out)
    Like_b  =read_outfiles_like(b_out)
    Like_Sigma = read_outfiles_like(S_out)
    N_acc, N_tot, Pct = read_outfiles_acc(Acc_out)
    Mu_hyp, B_hyp, Sigma_hyp , IT_hyp= read_outfile_hypercube(Espace_out)
    # On passe en b_value
    # B_hyp[0]=B_hyp[0]/np.log(10)
    # B_hyp[1]=B_hyp[1]/np.log(10)
    Nc_vec=Nc_vecs[:,0]
    N_iter = len(Nc_vec)
    N_tirs = int(IT_hyp[1])
    N_mags = len(mags)
   
    tmax = max(times)
    if N_mags ==1 :
        tmax= max(times)*2
        
        
        
        
# =============================================================================
#    B POSITIVE     
# =============================================================================

    window_1 = 50
    window_2 = 200
    window_3 = 500
    beta_positive_100 = Bpos_function(mags, window_1)
    #beta_positive_200 = Bpos_function(mags, window_2)
    #beta_positive_400 = Bpos_function(mags, window_3)

    
    
    
    #plt.plot(times[int(window/2):-int(window/2)], beta_positive)
# =============================================================================
#   Adapting the name of the figures in function of the run
# =============================================================================
    Num_run = 1 
    today = date.today()
    Date = today.strftime("%d%m")
    Title1 = 'Transdimensional RJ-MCMC - Synthetic dataset '+str(N_mags)+' evts'
    Subtitle =  today.strftime("%d/%m") +' : Run '+str(Num_run)+' - '+str(N_iter) +' models'  
    Outname = name+Date+'_'+str(N_iter)+'it_'+str(N_tirs)+'rd_'
    
    
    
    
    
# =============================================================================
#    Plotting variations of distribution of Mc 
# =============================================================================
    Mu_dist = np.linspace(Mu_hyp[0], Mu_hyp[1], int(Mu_hyp[2]))
   
    
    fig0 = plt.figure(figsize=(10,10))
    plt.subplots_adjust(left=0, bottom=0, right=1.5, top=1, wspace=0, hspace=0)
    ax0 = fig0.add_subplot(111, projection='3d')
    Mean_mu= []
    Std_mu=[]
    for i in range (len(Prop_vec)):
        ax0.plot(Mu_dist, Like_mu[:,i],i, zdir='x', c='black',alpha=0.5)
        ax0.scatter(Mu_dist, Like_mu[:,i],i, zdir='x', c=Like_mu[:,i],cmap='viridis',alpha=0.5)
        meanmu= np.sum([Like_mu[j,i]* Mu_dist[j] for j in range(int(Mu_hyp[2])) ])/sum(Like_mu[:,i])
        stdmu = math.sqrt(np.sum( [( Mu_dist[j]-meanmu)**2* Like_mu[j,i] for j in range(int(Mu_hyp[2])) ]) /sum(Like_mu[:,i]))
        #ax0.scatter(meanmu,max(Like_mu[:,i]),i,zdir='x',c='red')
        Mean_mu.append( meanmu)
        Std_mu.append(stdmu)
    ax0.set_ylim(Mu_hyp[0],Mu_hyp[1])
    ax0.set_xlabel('Times',fontsize=12)
    ax0.set_ylabel('Exploration of mu',fontsize=12)
    ax0.set_zlabel('Likelihood',fontsize=12)
    fig0.savefig('Exploration_mu.pdf')
    #plt.show()
    
    
    
    
    B_dist = np.linspace(B_hyp[0]/np.log(10), B_hyp[1]/np.log(10), int(B_hyp[2]))
    
    fig0b = plt.figure(figsize=(10,10))
    plt.subplots_adjust(left=0, bottom=0, right=1.5, top=1, wspace=0, hspace=0)
    ax0b = fig0b.add_subplot(111, projection='3d')
    Mean_b=[]
    Std_b=[]
    for i in range (len(Prop_vec)):
        ax0b.plot(B_dist, Like_b[:,i],i,zdir='x', c='navy',alpha=0.5)
        ax0b.add_collection3d(plt.fill_between(B_dist,Like_b[:,i],-0.1, color='orange', alpha=0.3,label="filled plot"),1, zdir='z')
        ax0b.scatter(B_dist, Like_b[:,i],i, zdir='x', c=Like_b[:,i],cmap='viridis',alpha=0.5)
        meanb=np.sum([Like_b[j,i]* B_dist[j] for j in range(int(B_hyp[2])) ])/sum(Like_b[:,i])
        stdb = math.sqrt(np.sum( [( B_dist[j]-meanb)**2* Like_b[j,i] for j in range(int(B_hyp[2])) ]) /sum(Like_b[:,i]))
        #ax0b.scatter(meanb,max(Like_b[:,i]),i,zdir='x',c='red')
        Mean_b.append( meanb)
        Std_b.append(stdb)
    ax0b.set_ylim(B_hyp[0]/np.log(10),B_hyp[1]/np.log(10))
    ax0b.set_xlabel('Times',fontsize=12)
    ax0b.set_ylabel('Exploration of beta',fontsize=12)
    ax0b.set_zlabel('Likelihood',fontsize=12)
    fig0b.savefig('Exploration_b.pdf')
    #plt.show()
    
    
    
    Sigma_dist = np.linspace(Sigma_hyp[0], Sigma_hyp[1], int(Sigma_hyp[2]))
    
    fig0c = plt.figure(figsize=(10,10))
    plt.subplots_adjust(left=0, bottom=0, right=1.5, top=1, wspace=0, hspace=0)
    ax0c = fig0c.add_subplot(111, projection='3d')
    Mean_S=[]
    Std_S=[]
    for i in range (len(Prop_vec)):
        ax0c.plot(Sigma_dist, Like_Sigma[:,i],i,zdir='x', c='navy',alpha=0.5)
        ax0c.scatter(Sigma_dist, Like_Sigma[:,i],i, zdir='x', c=Like_Sigma[:,i],cmap='viridis',alpha=0.5)
        meanSigma=np.sum([Like_Sigma[j,i]* Sigma_dist[j] for j in range(int(Sigma_hyp[2])) ])/sum(Like_Sigma[:,i])
        stdSigma = math.sqrt(np.sum( [( Sigma_dist[j]-meanSigma)**2* Like_Sigma[j,i] for j in range(int(Sigma_hyp[2])) ]) /sum(Like_Sigma[:,i]))
        #ax0b.scatter(meanb,max(Like_b[:,i]),i,zdir='x',c='red')
        Mean_S.append( meanSigma)
        Std_S.append(stdSigma)
    ax0c.set_ylim(Sigma_hyp[0]/np.log(10),Sigma_hyp[1]/np.log(10))
    ax0c.set_xlabel('Times',fontsize=12)
    ax0c.set_ylabel('Exploration of Sigma',fontsize=12)
    ax0c.set_zlabel('Likelihood',fontsize=12)
    fig0c.savefig('Exploration_Sigma.pdf')
    #plt.show()

# =============================================================================
# Bins we use for plotting 
# =============================================================================
    
    bins = np.linspace(0+year0,tmax+year0, len(Prop_vec))
    binsb = np.linspace(0+year0,tmax+year0, 20)
    bins2=np.linspace(0,len(Nc_vec), len(Nc_vec))
    bins4=np.arange(0,max(Nc_vec),1)
    
# =============================================================================
#   Retrieve the location of the discontinuities
# =============================================================================
    Nc_hist = plt.hist(Nc_vec, bins=np.arange(0,20,1))
    Nc_max = Nc_hist[1][np.where(Nc_hist[0]==max(Nc_hist[0]))][0]
    nbins_cat = float(max(times))/len(Prop_vec)
    disc_calc=[]
    liste_pics,dict_pics = ss.find_peaks(Prop_vec,height=1,prominence=1)
    prom_pics = dict_pics['prominences']
    for maxi in range(int(Nc_max)-1):
        index_1 = np.where(prom_pics == max(prom_pics))
        index_2 = liste_pics[index_1]
        prom_pics[index_1]=0
        disc_calc.append(index_2 * nbins_cat)
    



# =============================================================================
#  Look at the fits Gumbel+Richter at each discontinuities
# =============================================================================

    # figure_fit(Dataset, disc_calc, 1, 0.72,0.66,0.32)
    # figure_fit(Dataset, disc_calc, 2, 0.50,0.78,0.32)
    # #figure_fit(Dataset, disc_calc, 3, 0.6,0.72,0.32)
    # #figure_fit(Dataset, disc_calc, 4, 0.6,0.78,0.32)

    
    
# =============================================================================
#  Decompose the input catalogue at each discontinuities
# =============================================================================  
    # times_d1, mags_d1 = read_dataset_disc(Dataset, disc_calc, 1)
    # times_d2, mags_d2 = read_dataset_disc(Dataset, disc_calc, 2)
    # times_d3, mags_d3 = read_dataset_disc(Dataset, disc_calc, 3)

    
    
    
    
    import matplotlib.dates as mdates
    import datetime as datetime
    from cmcrameri import cm
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
    
    fig0, ax0= plt.subplots(1,1,figsize=(10,8))
    ax0b = ax0.twinx()
    
    ax0b.scatter(times,mags, c=times, s=[exp(magi)*2 for magi in mags],edgecolor='black', lw=0.25,cmap=cm.vik, alpha=0.8)
    ax0b.set_ylabel('Magnitudes',fontsize=17, c='dimgrey')

    ax0b.set_ylim(0,8)
    ax0.set_ylim(0,1)
    ax0.set_xlim(min(times), max(times))
    ax0.tick_params(axis='x', colors='black',labelrotation=20,labelsize=17,grid_markersize=17)
    found_disc = [y/sum(N_acc) for y in Prop_vec]
    Prop_norm = [y/max(Prop_vec) for y in Prop_vec]
 
    ax0.plot(bins, Prop_norm,lw=3, zorder=2,c='black', label='Inversion McMC')
    if synthetics=='True':
        for kk in range (len(T_vec)):
            if kk==0 :
                ax0.vlines(T_vec[kk]+year0, 0,1, linestyle=':', lw=2,color='green', label='Fixed discontinuities')
            ax0.vlines(T_vec[kk]+year0, 0, 1, linestyle=':', lw=2,color='green')
    ax0.set_xlabel('Time', fontsize=17)
    ax0.tick_params(axis='y', colors='black',labelsize=17,grid_markersize=17)
    ax0.tick_params(axis='x', colors='black',labelsize=17,grid_markersize=17)
    ax0.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax0b.tick_params(axis='y',colors='dimgrey',labelsize=17,grid_markersize=17)
    ax0.set_ylabel('Probability of discontinuities ',c='black', fontsize=17)
    #ax0.grid(lw=0.3)
    #ax0b.legend()
    ax0.legend(loc='upper right', fontsize=17)
    fig0.savefig(Outname+'Fig0.pdf')
    #plt.show()
# =============================================================================
#     FIGURE 1 : Informative plot on acceptance rates
# =============================================================================
    fig1, axs = plt.subplots(2,2,figsize=(10,8))
    ax1= axs[0,0]
    ax2= axs[1,0]
    ax3= axs[0,1]
    ax4= axs[1,1]
    plt.subplots_adjust(wspace=0.4, hspace=0.2)
    

    found_disc = [y/sum(N_acc) for y in Prop_vec]
    ax1.plot(bins, Prop_norm,zorder=2,c='black')
    if synthetics=='True':
        for kk in range (len(T_vec)):
            ax1b.vlines(T_vec[kk]+year0, 0, 1, lw=1.5, linestyle='--', color='green')
    ax1.set_xlabel('Dataset times')
    ax1.set_ylabel('Nb of acceptance for that time')
    ax1.grid(lw=0.3)
    
    #ax2.scatter(bins2, Nc_vec, c='firebrick',s=0.2)
    color = iter(cm.grayC(np.linspace(0, 1, len(Nc_vecs[0]))))
    for ni in range(len(Nc_vecs[0])):
        colori = next(color)
        ax2.plot(bins2,Nc_vecs[:,ni],c=colori, lw=0.6, alpha=0.5)
    mean_Nc = [np.mean(Nc_vecs[i,:]) for i in range(len(Nc_vecs))]
    # std_Nc_pos = [mean_Nc[i]+np.std(Nc_vecs[i,:]) for i in range(len(Nc_vecs))]
    # std_Nc_neg = [mean_Nc[i]-np.std(Nc_vecs[i,:]) for i in range(len(Nc_vecs))]
    # ax2.plot(bins2,mean_Nc,c='black', lw=0.8, alpha=0.9)
    #ax2.plot(bins2,std_Nc,c='black', lw=0.8, alpha=0.9)
    # ax2.fill_between(bins2, std_Nc_neg, std_Nc_pos,
    #              color='dimgrey', alpha=0.6,edgecolor=None)
    ax2.set_xlabel('Iteration number')
    ax2.set_ylabel('Nb of accepted discontinuities Nc')
    ax2.grid(lw=0.6)
    #ax2.set_ylim(0, max(Nc_vec)+2)
    
    
    lab = 'Nb total : %s models'%(str(sum(N_tot)))
    lab2 = 'Nb accepted : %s models'%(str(sum(N_acc)))
    ax3.bar(['Births','Deaths','Perts'], N_tot, color='darkslateblue',edgecolor='black',lw=0, alpha=0.5,label=lab)
    ax3.bar(['Births','Deaths','Perts'], N_acc,color='greenyellow',edgecolor='black',lw=2,alpha=0.6, label=lab2)
    ax3.legend(loc='upper center')
    for i in range(3):
        ax3.text(i-0.25,N_tot[i]/3,'Na='+str(N_acc[i]))
        ax3.text(i-0.25,N_tot[i]/5,'Ta='+str(round(Pct[i]*100,1))+u'%')
    ax3.set_ylabel('Nb of tot/acc models')
    
    
    ax4.hist(Nc_vecs.flatten(), bins=bins4, color='firebrick', alpha=0.6, edgecolor='darkslateblue',lw=2, orientation='horizontal',align='left',label='Nc')
    ax4.set_ylim(0,max(Nc_vec)+2)
    ax4.set_ylabel('Nb of accepted discontinuities Nc')
    ax4.set_xlabel('Number of models')
    ax4.legend(loc='upper right')
    fig1.suptitle(Title1 +'\n'+Subtitle,fontweight="bold")
    fig1.savefig(Outname+'Fig1.pdf')
    #plt.show()
    
# =============================================================================
#     FIGURE 2 - Likelihood Mc 
# =============================================================================
    fig2, ax = plt.subplots(1,figsize=(8,8))
    im =ax.imshow(Like_mu, aspect='auto', cmap='Reds', extent = (0+year0,tmax+year0,Mu_hyp[1],Mu_hyp[0]))


    
    ax.plot(bins, Mean_mu,c='red',linestyle='-.',lw=1,alpha=0.7)
    ax.plot(bins,[ Mean_mu[k]+Std_mu[k] for k in range(len(Mean_mu))],c='black',lw=1.5, linestyle='--',alpha=0.7)
    ax.plot(bins,[ Mean_mu[k]-Std_mu[k] for k in range(len(Mean_mu))],c='black',lw=1.5, linestyle='--',alpha=0.7)
    #ax.vlines(disc_calc2,Mu_hyp[0],Mu_hyp[1], color='gray',ls=':',lw=2)
    ax2 = ax.twinx()
    #ax2.scatter(times2,mags, s=[np.exp(magi) for magi in mags], alpha=0.05,color='white')
    ax2.set_ylim(0,7)
    ax2.set_ylabel('Magnitudes MLv',c='dimgrey',fontsize=18)
    #axbb.invert_yaxis()
    #axbb.set_ylim(30,0)
    #ax.invert_yaxis()
    ax.set_ylim(Mu_hyp[0],Mu_hyp[1])
    ax.grid(lw=0.3)
   #ax.vlines(disc_calc,0,4, color='white',linestyles=':',lw=2)
    if synthetics=='True':
        for di in range(len(T_vec)-1):
            ax.hlines(mu_vec[di],T_vec[di],T_vec[di+1], color='red',linestyles='--',lw=3)
            ax.vlines(T_vec[di]+year0, Mu_hyp[0],Mu_hyp[1], linestyle=':', lw=3,color='lightgreen', label='Fixed discontinuities')
    ax.spines['left'].set_color('red')
    ax.tick_params(axis='y', colors='firebrick',labelsize=18,grid_markersize=18)
    ax.tick_params(axis='x', colors='black',labelsize=18,grid_markersize=18)
    #ax2.tick_params(axis='y',colors='dimgrey',labelsize=14,grid_markersize=14)
    ax.set_ylabel('Values of $\mu$',c='firebrick',fontsize=18)
    ax.set_xlim(0,100)
    ax.set_xlabel('Time', fontsize=18)
    Title2 = 'Likelihood $\mu$ - Synthetic dataset '+str(N_mags)+' evts'
    fig2.suptitle(Title2 +'\n'+Subtitle,fontweight="bold")
    fig2.savefig(Outname+'_Fig2_likelihoodMu.pdf')
    #plt.show()
    
# # =============================================================================
# #     FIGURE 3 - Likelihood b
# # =============================================================================   
    b_true = Like_b/np.log(10)
    fig3, axb = plt.subplots(1,figsize=(8,11))
    imb = axb.imshow(b_true, aspect='auto', cmap='Reds', extent = (0+year0,tmax+year0,B_hyp[1]/np.log(10),B_hyp[0]/np.log(10)))

    axb.plot(bins, Mean_b,c='red',linestyle='-.', lw=1,alpha=0.9, label='b-value moyenne')
    
    axb.plot(bins,[ Mean_b[k]+Std_b[k] for k in range(len(Mean_b))],c='black',lw=1, linestyle='--',alpha=0.7, label="Ecart-type 1Sigma")
    axb.plot(bins,[ Mean_b[k]-Std_b[k] for k in range(len(Mean_b))],c='black',lw=1, linestyle='--',alpha=0.7)
    #axb.plot(times[int(window_1/2):-int(window_1/2)], beta_positive_100, c='lightblue', lw=2, label="Bpos - "+str(window_1))
    #axb.plot(times[int(window_2/2):-int(window_2/2)], beta_positive_200, c='lightblue', lw=2, label="Bpos - "+str(window_2))
    #axb.plot(times[int(window_3/2):-int(window_3/2)], beta_positive_400, c='darkblue', lw=2, label="Bpos - "+str(window_3))
    
    axb2 = axb.twinx()
    #axb2.scatter(times2,mags, s=[np.exp(magi) for magi in mags], alpha=0.05,color='white')
    axb2.set_ylim(0,7)
    axb2.set_ylabel('Magnitudes MLv',c='dimgrey',fontsize=18)
    #axb2.tick_params(axis='y',colors='dimgrey',labelsize=14,grid_markersize=14)
    #/axb.scatter(t,m, s=3, alpha=0.2,color='white')
    #axb.invert_yaxis()
    axb.set_ylim(B_hyp[0]/np.log(10),B_hyp[1]/np.log(10))
    axb.grid(lw=0.3)
    
    axb.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    axb.spines['left'].set_color('firebrick')
    axb.tick_params(axis='y', colors='firebrick',labelsize=18,grid_markersize=18)
    axb.tick_params(axis='x', colors='black',labelsize=18,grid_markersize=18)
    
    axb.set_ylabel('b values',c='firebrick', fontsize=14)
    axb.set_xlabel('Time',fontsize=18)
    #axb.vlines(disc_calc2,B_hyp[0]/np.log(10),B_hyp[1]/np.log(10), color='grey',ls=':',lw=2)
    axb.set_xlim(0,100)
    if synthetics=='True':
        for di in range(len(T_vec)-1):
            if di==0:
                axb.hlines(b_vec[di]/np.log(10),T_vec[di],T_vec[di+1], color='red',linestyles='--',lw=3, label='Valeurs imposees')
            axb.hlines(b_vec[di]/np.log(10),T_vec[di],T_vec[di+1], color='red',linestyles='--',lw=3)
            axb.vlines(T_vec[di]+year0, B_hyp[0]/np.log(10),B_hyp[1]/np.log(10), linestyle=':', lw=3,color='lightgreen', label='Fixed discontinuities')
    Title3 = 'Likelihood b - Synthetic dataset '+str(N_mags)+' evts'
    fig3.colorbar(imb, ax=axb,orientation = 'horizontal')
    fig3.suptitle(Title3 +'\n'+Subtitle,fontweight="bold")
    fig3.savefig(Outname+'_Fig3_likelihoodb.pdf')
    #plt.show()
    
    
# # =============================================================================
# #     FIGURE 4 - Likelihood S
# # =============================================================================   

    fig4, axs = plt.subplots(1,figsize=(8,8))
    ims = axs.imshow(Like_Sigma, aspect='auto', cmap='Reds', extent = (0+year0,tmax+year0, Sigma_hyp[1],Sigma_hyp[0]))
    
    axs.plot(bins, Mean_S,c='red',linestyle='-.',lw=1)
    axs.plot(bins,[ Mean_S[k]+Std_S[k] for k in range(len(Mean_S))],c='black',lw=1, linestyle='--',alpha=0.7)
    axs.plot(bins, [ Mean_S[k]-Std_S[k] for k in range(len(Mean_S))],c='black',lw=1, linestyle='--',alpha=0.7)
    axs.set_xlim(0,100)
    axs2 = axs.twinx()
    #axs2.scatter(times2,mags, s=[np.exp(magi) for magi in mags], alpha=0.05,color='white')
    axs2.set_ylim(0,7)
    axs2.set_ylabel('Magnitudes MLv',c='dimgrey', fontsize=18)
    
    #axs.vlines(disc_calc2,Sigma_hyp[0],Sigma_hyp[1], color='gray',ls=':',lw=2)
    #/axb.scatter(t,m, s=3, alpha=0.2,color='white')
    #axb.invert_yaxis()
    axs.set_ylim(Sigma_hyp[0],Sigma_hyp[1])
    axs.grid(lw=0.3)
    
    
    axs.spines['left'].set_color('firebrick')
    axs.tick_params(axis='y', colors='firebrick',labelsize=18,grid_markersize=18)
    axs.tick_params(axis='x', colors='black',labelsize=18,grid_markersize=18)
    #axs2.tick_params(axis='y',colors='dimgrey',labelsize=14,grid_markersize=14)
    axs.set_ylabel('Values of Sigma',c='firebrick',fontsize=18)
    axs.set_xlabel('Time',fontsize=18)
    #axb.vlines(disc_calc,0,8, color='white',ls=':',lw=2)
    if synthetics=='True':
        for di in range(len(T_vec)-1):
            axs.hlines(sigma_vec[di],T_vec[di],T_vec[di+1], color='red',linestyles='--',lw=3)
            axs.vlines(T_vec[di]+year0, Sigma_hyp[0],Sigma_hyp[1], linestyle=':', lw=3,color='lightgreen', label='Fixed discontinuities')
    Title4 = 'Likelihood Sigma - Synthetic dataset '+str(N_mags)+' evts'
    fig4.suptitle(Title4 +'\n'+Subtitle,fontweight="bold")
    fig4.savefig(Outname+'_Fig4_likelihood_Sigma.pdf')
    #plt.show()