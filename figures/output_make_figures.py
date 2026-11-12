#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:54:30 2023
@author: Marine Laporte - Université Claude Bernard Lyon 1
@last-update: 26/10/2026
@contact: marinelaporte.ensg@gmail.com
"""

## Libraries we use 
## ------------------
import sys 
import os
import time

import plotly.express as px
import numpy as np
import pandas as pd
from math import exp,sqrt, erf
from mpmath import mp
from datetime import date 
import scipy.signal as ss
import scipy.stats as st
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

from matplotlib import cm
from matplotlib.pyplot import cm
from cmcrameri import cm
## ------------------
from imports_outfiles_ import * 



def datetime_to_decimal_year_fast(series):
    year = series.dt.year
    start_of_year = pd.to_datetime({'year': year, 'month': 1, 'day': 1})
    start_of_next = pd.to_datetime({'year': year + 1, 'month': 1, 'day': 1})
    fraction = (series - start_of_year) / (start_of_next - start_of_year)
    return year + fraction


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
            E_max[j]=(0.5+(erf((x[j+1]-m_mle)/(sqrt(2)*s_mle))/2)) *beta_mle * exp(-beta_mle*x[j+1])
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



def B_Positive_function(times, mags, w):
    """Fonction qui calcule la valeur de la b-value
    a partir de la formule B positive du papier de Van Der Elst (2021)
    utilisant les différences positives de magnitudes suivant la distribution de Laplace
    Inputs :
        - mags : liste ou array de magnitudes (floats)
        - times : liste ou array des temps (floats) associés aux magnitudes 
        - w : taille de la fenetre (nb d'evenements/fenetre)
    Outputs :
        - b_positive : liste ou array de valeurs de b-positives
        - times_b_positive : temps associés aux valeurs de b_positive """
    # INITIALISATION
    n_events =len(mags)
    differences_magnitudes = np.zeros(n_events-1)
    differences_absolues = np.zeros(n_events-1)
    # init outputs
    b_positive= np.zeros(n_events-w)
    times_b_positive= np.zeros(n_events-w)
    
    # Calcul des differences de magnitudes sur tout le catalogue
    for i in range (0, n_events-1):
        differences_magnitudes[i] = mags[i+1] - mags[i]       
        differences_absolues[i] = abs(mags[i+1] - mags[i])
        
    for m in range(len(differences_magnitudes)-(w-1)):
       differences_positives=[]
       flag_pos =0
       flag =0
       Minimum_difference = min(differences_absolues[m:m+w]+0.4)
       for valeur in differences_magnitudes[m:m+w]:
           flag+=1
           if valeur > Minimum_difference:
               flag_pos+=1 # Pour savoir sur combien d evenements on calcule effectivement la b-value
               differences_positives.append(valeur)
       windowed_mean = sum(np.ones(len(differences_positives)) * differences_positives) / len(differences_positives)
       
       
       # d'apres la formule de Laplace (Van der Elst 2021)
       windowed_beta = 1/float(windowed_mean-Minimum_difference) 
       b_positive[m]=windowed_beta/np.log(10)  # On retrouve b_value = beta/2.3

    print("Valeur de b-value obtenue sur %s events sur fenetre de %s events" %(str(flag_pos),str(w)))
    times_b_positive = times[int(w/2):-int(w/2)]
    return times_b_positive, b_positive



if __name__ == "__main__" :
    
    
    today = time.strftime('%d_%m_%Y_%Hh')
    name = 'Test'
    year0 = 0      # input.txt contient deja les temps absolus : ne rien ajouter
    synthetics='True'
    out_path = 'Run_'+ today
    os.makedirs(out_path, exist_ok=True)
# =============================================================================
#   Default outfiles names in OUT/
# =============================================================================
    Prop_out ='../outputs/prop.txt'
    Nc_out ='../outputs/NC_tot.txt'
    Dataset ='../outputs/input.txt'
    Mc_out ='../outputs/Likelihood_Mu.txt'
    b_out ='../outputs/Likelihood_b.txt'
    S_out ='../outputs/Likelihood_Sigma.txt'
    Acc_out ='../outputs/Acceptances.txt'
    Espace_out ='../outputs/hypercube.txt'
    Beta_tot ='../outputs/Beta_tot.txt'
# =============================================================================
#   If run on synthetics - give synthetic parameters to retrieve
# =============================================================================

    T_vec = [0, 10, 20, 40, 60, 70, 88, 100]
    N_vec = [394,467,260,2723,657,1070,431]  
    
    b_vec =  np.array([0.7,1.1, 1.0, 1.0, 0.9, 0.7, 0.65])*np.log(10)
    mu_vec = np.array([0.8, 0.8,1.5, 0.5, 0.7, 0.9, 0.9])
    sigma_vec = np.array([0.35, 0.3, 0.2, 0.15,0.3,0.4, 0.15])
    
    # T_vec = [0,  100]
    # N_vec = [2987]  
    
    # b_vec =  np.array([0.7])*np.log(10)
    # mu_vec = np.array([0.8])
    # sigma_vec = np.array([0.35])
    
    


    
    Ntot = sum(N_vec)
# =============================================================================
#     Reading outfiles  
# =============================================================================

    Prop_vec = read_outfiles_prop(Prop_out)
    Nc_vecs =read_outfiles_like(Nc_out)
    times, mags = read_dataset(Dataset)
    times2 = [a+year0 for a in times]
    mags = np.asarray(mags)
    Like_mu = read_outfiles_like(Mc_out)
    Like_b  =read_outfiles_like(b_out)
    Like_Sigma = read_outfiles_like(S_out)
    N_acc, N_tot, N_pct = read_outfiles_acc(Acc_out)
    with open(Espace_out) as f:
        hyp_lines = [np.array(l.split(), dtype=float) for l in f if l.strip()]
    Mu_hyp, B_hyp, Sigma_hyp, IT_hyp = hyp_lines[0], hyp_lines[1], hyp_lines[2], hyp_lines[3]
    t_offset, time_max = hyp_lines[4]
    beta_store, mu_store, sigma_store = hyp_lines[5]
    #Beta_vecs = read_outfiles_like(Beta_tot)
    
    Nc_vec=Nc_vecs[:,0]
    N_iter = len(Nc_vec)
    N_tirs = int(IT_hyp[1])
    N_mags = len(mags)
    T_bins = len(Prop_vec)
    dt_bin = time_max / T_bins
    t_edges = (t_offset + year0, t_offset + time_max + year0)   # bords de l'axe des temps des inversions

    
# =============================================================================
#   Adapting the name of the figures in function of the run
# =============================================================================
    Num_run = 1 
    today = date.today()
    Date = today.strftime("%d%m")
    Title1 = 'Transdimensional RJ-MCMC - Synthetic dataset '+str(N_mags)+' evts'
    Subtitle =  today.strftime("%d/%m") +' : Run '+str(Num_run)+' - '+str(N_iter) +' models'  
    out_name = name+'_'+str(N_iter)+'it_'
    
# =============================================================================
#    Estimating means and stds 
# =============================================================================
    def weighted_mean_std(values, weights):
     """Moyenne et ecart-type ponderes par colonne (un bin de temps par colonne)."""
     totals = weights.sum(axis=0)
     totals = np.where(totals > 0, totals, np.nan)          # bin vide -> NaN au lieu d'une erreur
     mean = (weights * values[:, None]).sum(axis=0) / totals
     std = np.sqrt(((values[:, None] - mean) ** 2 * weights).sum(axis=0) / totals)
     return mean, std

    Mu_dist    = Mu_hyp[0] + (np.arange(int(Mu_hyp[2])) + 0.5) * mu_store
    B_dist     = (B_hyp[0] + (np.arange(int(B_hyp[2])) + 0.5) * beta_store) / np.log(10)
    Sigma_dist = Sigma_hyp[0] + (np.arange(int(Sigma_hyp[2])) + 0.5) * sigma_store

    Mean_mu, Std_mu = weighted_mean_std(Mu_dist, Like_mu)
    Mean_b,  Std_b  = weighted_mean_std(B_dist, Like_b)
    Mean_S,  Std_S  = weighted_mean_std(Sigma_dist, Like_Sigma)

# =============================================================================
# Bins we use for plotting 
# =============================================================================
    
    bins = t_offset + year0 + (np.arange(T_bins) + 0.5) * dt_bin    # centres des bins temporels
    bins2=np.linspace(0,len(Nc_vec), len(Nc_vec))
    bins4=np.arange(-0.5,max(Nc_vec)+1.5,1)
    #bins5=np.linspace(0,len(Beta_vecs[:,0]), len(Beta_vecs[:,0]))
    
# =============================================================================
#   Retrieve the location of the discontinuities
# =============================================================================
    burnin_frac = 0.5                                    # = it_burnin / it_max
    Nc_post = Nc_vec[int(len(Nc_vec)*burnin_frac):].astype(int)
    Nc_max = np.bincount(Nc_post).argmax()               # Nc le plus frequent apres burn-in
    liste_pics, dict_pics = ss.find_peaks(Prop_vec, height=1, prominence=1)
    best = np.argsort(dict_pics['prominences'])[::-1][:Nc_max]
    disc_calc = np.sort(t_offset + year0 + (liste_pics[best] + 0.5) * dt_bin)
    print('Nc le plus probable :', Nc_max)
    print('Discontinuites (annees) :', np.round(disc_calc, 2))
    

    
    
# =============================================================================
#  Decompose the input catalogue at each discontinuities
# =============================================================================  
    # times_d1, mags_d1 = read_dataset_disc(Dataset, disc_calc, 1)
    # times_d2, mags_d2 = read_dataset_disc(Dataset, disc_calc, 2)
    # times_d3, mags_d3 = read_dataset_disc(Dataset, disc_calc, 3)

    
    fig0, ax0= plt.subplots(1,1,figsize=(10,8))
    ax0b = ax0.twinx()
    
    ax0b.scatter(times2,mags, c=times, s=[exp(magi)*2 for magi in mags],edgecolor='black', lw=0.25,cmap=cm.vik, alpha=0.8)
    ax0b.set_ylabel('Magnitudes',fontsize=17, c='dimgrey')

    ax0b.set_ylim(mags.min()-0.5, mags.max()+0.5)
    
    ax0.set_xlim(min(times2), max(times2))
    ax0.tick_params(axis='x', colors='black',labelrotation=20,labelsize=17,grid_markersize=17)
    found_disc = [y/sum(N_acc) for y in Prop_vec]
    n_models_post = (IT_hyp[0] * (1 - burnin_frac)) / 100     # (it_max - it_burnin) / it_thin
    Prop_norm = Prop_vec / n_models_post                     # probabilite d'une discontinuite par bin
 
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
    ax0.legend(loc='upper right', fontsize=17)
    ax0.set_ylim(0, max(0.05, Prop_norm.max()*1.1))
    fig0.savefig(os.path.join(out_path , out_name +'Fig0.pdf'))
    plt.show()
# =============================================================================
#     FIGURE 1 : Informative plot on acceptance rates
# =============================================================================
    fig1, axs = plt.subplots(2,2,figsize=(15,12))
    ax1= axs[0,0]
    ax2= axs[1,0]
    ax3= axs[0,1]
    ax4= axs[1,1]
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    

    found_disc = [y/sum(N_acc) for y in Prop_vec]
    ax1.plot(bins, Prop_vec,zorder=2,c='black')
    if synthetics=='True':
        ax1b = ax1.twinx()
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
    
    
    lab = 'Nb total : %s models'%(str(int(N_tot[0]+N_tot[4])))
    lab2 = 'Nb accepted : %s models'%(str(int(N_acc[0]+N_acc[4])))
    ax3.bar(['GEO', 'Births','Deaths','Perts','GR','Beta','Mu','Sigma'], N_tot, color='darkslateblue',edgecolor='black',lw=0, alpha=0.5,label=lab)
    ax3.bar(['GEO', 'Births','Deaths','Perts','GR','Beta','Mu','Sigma'], N_acc,color='greenyellow',edgecolor='black',lw=2,alpha=0.6, label=lab2)
    ax3.legend(loc='upper center')
    for i in range(8):
        ax3.text(i-0.5,N_tot[i]/3,str(int(N_acc[i])))
        ax3.text(i-0.5,N_tot[i]/5,str(round(N_pct[i],1))+u'%')
    ax3.set_ylabel('Nb of tot/acc models')
    ax3.tick_params(axis='x', rotation=45)
    
    ax4.hist(Nc_vecs, bins=bins4, color='firebrick', alpha=0.6, edgecolor='darkslateblue',lw=2, orientation='horizontal',align='mid',label='Nc')
    ax4.set_ylim(0,max(Nc_vec)+2)
    ax4.set_ylabel('Nb of accepted discontinuities Nc')
    ax4.set_xlabel('Number of models')
    ax4.legend(loc='upper right')
    fig1.suptitle(Title1 +'\n'+Subtitle,fontweight="bold")
    fig1.savefig(os.path.join(out_path , out_name +'Fig1.pdf'))
    plt.show()
    
# =============================================================================
#     FIGURE 2 - Likelihood Mu
# =============================================================================
    fig2, ax = plt.subplots(1,figsize=(8,8))
    im =ax.imshow(Like_mu, aspect='auto', cmap='Reds', extent = (t_edges[0], t_edges[1], Mu_hyp[1], Mu_hyp[0]))


    
    ax.plot(bins, Mean_mu,c='red',linestyle='-.',lw=1,alpha=0.7)
    ax.plot(bins,[ Mean_mu[k]+Std_mu[k] for k in range(len(Mean_mu))],c='black',lw=1.5, linestyle='--',alpha=0.7)
    ax.plot(bins,[ Mean_mu[k]-Std_mu[k] for k in range(len(Mean_mu))],c='black',lw=1.5, linestyle='--',alpha=0.7)
    ax.vlines(disc_calc, Mu_hyp[0], Mu_hyp[1], color='gray', ls=':', lw=2)
    ax2 = ax.twinx()
    ax2.scatter(times2,mags, s=[5*np.exp(magi) for magi in mags], alpha=0.4,color='grey', edgecolor='black')
    ax2.set_ylim(mags.min()-0.5, mags.max()+0.5)
    ax2.set_ylabel('Magnitudes ML',c='dimgrey',fontsize=18)
    ax.set_ylim(Mu_hyp[0],Mu_hyp[1])
    ax.grid(lw=0.3)
    if synthetics=='True':
        for di in range(len(T_vec)-1):
            ax.hlines(mu_vec[di],T_vec[di]+year0,T_vec[di+1]+year0, color='red',linestyles='--',lw=3)
            ax.vlines(T_vec[di]+year0, Mu_hyp[0],Mu_hyp[1], linestyle=':', lw=3,color='lightgreen', label='Fixed discontinuities')
    ax.spines['left'].set_color('red')
    ax.tick_params(axis='y', colors='firebrick',labelsize=18,grid_markersize=18)
    ax.tick_params(axis='x', colors='black',labelsize=18,grid_markersize=18)
    ax.set_ylabel(r'Values of $\mu$',c='firebrick',fontsize=18)
    ax.set_xlabel('Time', fontsize=18)
    Title2 = r'Likelihood $\mu$ - Synthetic dataset '+str(N_mags)+' evts'
    fig2.suptitle(Title2 +'\n'+Subtitle,fontweight="bold")
    fig2.savefig(os.path.join(out_path , out_name +'_Fig2_likelihoodMu.pdf'))
    plt.show()
    
# # =============================================================================
# #     FIGURE 3 - Likelihood b
# # =============================================================================   
    b_true = Like_b/np.log(10)
    
    # --- Préparation des données d'injection ---

    injection_data = '/home/marine/Documents/2-Catalogues/12-Lacq/injection.txt'
    if os.path.exists(injection_data):
        df_inject = pd.read_csv(
            injection_data,
            sep=r'\s+', header=None, names=['date', 'value'])
        df_inject['date'] = pd.to_datetime(df_inject['date'], format='%d-%b-%Y')
        df_inject  = df_inject.set_index('date').sort_index()
        df_inject['decimal_year_abs'] = datetime_to_decimal_year_fast(df_inject.index.to_series())
        df_inject['decimal_year'] = (
        df_inject['decimal_year_abs'] - df_inject['decimal_year_abs'].iloc[0])
        
        window    = 30  # fenêtre de la moyenne glissante en jours
        threshold = df_inject['value'].quantile(0.99)  # seuil valeurs extrêmes (95ème percentile)
        df_inject_smooth  = df_inject['value'].rolling(window=window, center=True).mean()
        df_inject_extreme = df_inject['value'].where(df_inject['value'] >= threshold)  # NaN en dessous du seuil
        cumsum_extreme    = df_inject_extreme.dropna().cumsum()
    
    
    
    fig3, axb = plt.subplots(1,figsize=(8,10))
    imb = axb.imshow(b_true, aspect='auto', cmap='Reds', extent = (t_edges[0], t_edges[1], B_hyp[1]/np.log(10), B_hyp[0]/np.log(10)))

    axb.plot(bins, Mean_b,c='red',linestyle='-.', lw=1,alpha=0.9, label='b-value moyenne')
    
    axb.plot(bins,[ Mean_b[k]+Std_b[k] for k in range(len(Mean_b))],c='black',lw=1, linestyle='--',alpha=0.7, label="Ecart-type 1Sigma")
    axb.plot(bins,[ Mean_b[k]-Std_b[k] for k in range(len(Mean_b))],c='black',lw=1, linestyle='--',alpha=0.7)
    #axb.plot(times[int(window_1/2):-int(window_1/2)], beta_positive_100, c='lightblue', lw=2, label="Bpos - "+str(window_1))
    #axb.plot(times[int(window_2/2):-int(window_2/2)], beta_positive_200, c='lightblue', lw=2, label="Bpos - "+str(window_2))
    #axb.plot(times[int(window_3/2):-int(window_3/2)], beta_positive_400, c='darkblue', lw=2, label="Bpos - "+str(window_3))
    
    # axb2 = axb.twinx()
    # axb2.scatter(times2,mags, s=[4*np.exp(magi) for magi in mags], alpha=0.05,color='black')
    # axb2.set_ylim(0,4.5)
    # axb2.set_ylabel('Magnitudes MLv',c='dimgrey',fontsize=18)
    # axb2.tick_params(axis='y',colors='dimgrey',labelsize=14,grid_markersize=14)
    
    
    axb.set_ylim(B_hyp[0]/np.log(10),B_hyp[1]/np.log(10))
    axb.grid(lw=0.3)
    
    #axb.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    axb.spines['left'].set_color('firebrick')
    axb.tick_params(axis='y', colors='firebrick',labelsize=18,grid_markersize=18)
    axb.tick_params(axis='x', colors='black',labelsize=18,grid_markersize=18)
    
    axb.set_ylabel('b values',c='firebrick', fontsize=14)
    axb.set_xlabel('Time',fontsize=18)
    axb.vlines(disc_calc, B_hyp[0]/np.log(10), B_hyp[1]/np.log(10), color='grey', ls=':', lw=2)
    if synthetics=='True':
        for di in range(len(T_vec)-1):
            if di==0:
                axb.hlines(b_vec[di]/np.log(10),T_vec[di]+year0,T_vec[di+1]+year0, color='red',linestyles='--',lw=3, label='Valeurs imposees')
            axb.hlines(b_vec[di]/np.log(10),T_vec[di]+year0,T_vec[di+1]+year0, color='red',linestyles='--',lw=3)
            axb.vlines(T_vec[di]+year0, B_hyp[0]/np.log(10),B_hyp[1]/np.log(10), linestyle=':', lw=3,color='lightgreen', label='Fixed discontinuities')
    Title3 = 'Likelihood b - Synthetic dataset '+str(N_mags)+' evts'
    fig3.colorbar(imb, ax=axb,orientation = 'horizontal')
    fig3.suptitle(Title3 +'\n'+Subtitle,fontweight="bold")
    fig3.savefig(os.path.join(out_path , out_name +'_Fig3_likelihoodb.pdf'))
    plt.show()
    
    
# # =============================================================================
# #     FIGURE 4 - Likelihood S
# # =============================================================================   

    fig4, axs = plt.subplots(1,figsize=(8,8))
    ims = axs.imshow(Like_Sigma, aspect='auto', cmap='Reds', extent = (t_edges[0], t_edges[1], Sigma_hyp[1], Sigma_hyp[0]))
    
    axs.plot(bins, Mean_S,c='red',linestyle='-.',lw=1)
    axs.plot(bins,[ Mean_S[k]+Std_S[k] for k in range(len(Mean_S))],c='black',lw=1, linestyle='--',alpha=0.7)
    axs.plot(bins, [ Mean_S[k]-Std_S[k] for k in range(len(Mean_S))],c='black',lw=1, linestyle='--',alpha=0.7)
    axs2 = axs.twinx()
    axs2.set_ylim(0,7)
    axs2.set_ylabel('Magnitudes MLv',c='dimgrey', fontsize=18)
    
    axs.vlines(disc_calc, Sigma_hyp[0], Sigma_hyp[1], color='gray', ls=':', lw=2)
    axs.set_ylim(Sigma_hyp[0],Sigma_hyp[1])
    axs.grid(lw=0.3)
    
    
    axs.spines['left'].set_color('firebrick')
    axs.tick_params(axis='y', colors='firebrick',labelsize=18,grid_markersize=18)
    axs.tick_params(axis='x', colors='black',labelsize=18,grid_markersize=18)
    #axs2.tick_params(axis='y',colors='dimgrey',labelsize=14,grid_markersize=14)
    axs.set_ylabel('Values of Sigma',c='firebrick',fontsize=18)
    axs.set_xlabel('Time',fontsize=18)
    if synthetics=='True':
        for di in range(len(T_vec)-1):
            axs.hlines(sigma_vec[di],T_vec[di]+year0,T_vec[di+1]+year0, color='red',linestyles='--',lw=3)
            axs.vlines(T_vec[di]+year0, Sigma_hyp[0],Sigma_hyp[1], linestyle=':', lw=3,color='lightgreen', label='Fixed discontinuities')
    Title4 = 'Likelihood Sigma - Synthetic dataset '+str(N_mags)+' evts'
    fig4.suptitle(Title4 +'\n'+Subtitle,fontweight="bold")
    fig4.savefig(os.path.join(out_path , out_name +'_Fig4_likelihood_Sigma.pdf'))
    plt.show()
