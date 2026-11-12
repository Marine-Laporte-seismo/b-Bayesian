#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 11:18:29 2024

@author: marine
"""
import numpy as np
import sys
# =============================================================================
# # MAIN FUNCTIONS FOR IMPORT
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
        

        
    N_geometry_acc = float(lines[2].split()[3])
    N_Births_acc = float(lines[3].split()[2])
    N_Deaths_acc = float(lines[4].split()[2])
    N_Moves_acc = float(lines[5].split()[2])
    
    N_GR_acc = float(lines[6].split()[4])
    N_beta_acc = float(lines[7].split()[2])
    N_mu_acc = float(lines[8].split()[2])
    N_sigma_acc = float(lines[9].split()[2])
    
    
    
    
    N_geometry = float(lines[2].split()[4])
    N_Births = float(lines[3].split()[3])
    N_Deaths = float(lines[4].split()[3])
    N_Moves = float(lines[5].split()[3])
    
    N_GR = float(lines[6].split()[5])
    N_beta = float(lines[7].split()[3])
    N_mu = float(lines[8].split()[3])
    N_sigma = float(lines[9].split()[3])
    
    
    
    
    
    N_geometry_pct = float(lines[2].split()[5])
    N_Births_pct = float(lines[3].split()[4])
    N_Deaths_pct = float(lines[4].split()[4])
    N_Moves_pct = float(lines[5].split()[4])
    
    N_GR_pct = float(lines[6].split()[6])
    N_beta_pct = float(lines[7].split()[4])
    N_mu_pct = float(lines[8].split()[4])
    N_sigma_pct = float(lines[9].split()[4])

    N_acc = [N_geometry_acc, N_Births_acc, N_Deaths_acc, N_Moves_acc, N_GR_acc, N_beta_acc, N_mu_acc, N_sigma_acc]
    N_tot = [N_geometry, N_Births, N_Deaths, N_Moves, N_GR, N_beta, N_mu, N_sigma]
    N_pct = [N_geometry_pct, N_Births_pct, N_Deaths_pct, N_Moves_pct, N_GR_pct, N_beta_pct, N_mu_pct, N_sigma_pct]
    
    AR_GEO = [N_geometry_acc, N_geometry, N_geometry_pct]
    AR_Births = [N_Births_acc, N_Births, N_Births_pct]
    AR_Deaths=[N_Deaths_acc, N_Deaths, N_Deaths_pct]
    AR_Moves = [N_Moves_acc, N_Moves, N_Moves_pct]
    AR_GR = [N_GR_acc, N_GR, N_GR_pct]
    AR_beta = [N_beta_acc, N_beta, N_beta_pct]
    AR_mu = [N_mu_acc, N_mu, N_mu_pct]
    AR_sigma = [N_sigma_acc, N_sigma, N_sigma_pct]
    return N_acc, N_tot, N_pct

def read_outfile_hypercube_b(infile):
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

def read_outfile_hypercube_b_a(infile):
    import sys

    try:
        with open(infile, "r") as fid:
            lines = fid.readlines()
    except FileNotFoundError:
        print(f"Erreur : {infile} n'existe pas !")
        sys.exit(1)

    # Vérification : on attend 5 lignes
    if len(lines) < 5:
        print(f"Erreur : {infile} contient seulement {len(lines)} lignes (5 attendues).")
        sys.exit(1)

    # Lecture ligne par ligne
    Mu     = [float(x) for x in lines[0].split()]  # mu_min, mu_max, mu_bins, mu_init, mu_stdi
    B      = [float(x) for x in lines[1].split()]  # beta_min, beta_max, ...
    A      = [float(x) for x in lines[2].split()]  # a_min, a_max, ...
    Sigma  = [float(x) for x in lines[3].split()]  # sigma_min, sigma_max, ...
    
    # it_max et T_bins sont des entiers (tu peux adapter selon ton besoin)
    I_vals = [int(round(float(x))) for x in lines[4].split()]  # [it_max, T_bins, comp_time]

    return Mu, B, A, Sigma, I_vals
