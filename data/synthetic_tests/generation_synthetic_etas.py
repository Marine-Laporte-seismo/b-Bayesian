import numpy as np
import matplotlib.pyplot as plt
import time
from math import exp, sqrt, erf
from cmcrameri import cm
### implement the probability of one change in time



def synthetic_dataset(N_samples, mu, sigma, b_val, function='gumbel'):
    """ Creates a synthetic catalog of N_sampels values of magnitudes
    following the conjoint detection distribution (mu, sigma) 
    and the GR distribution associated to the input b_val
    The detection function can be chosen from :
        - full
        - error
        - gumbel
        - exp
        - expdouble"""
    
    i = 0
    M_distribution = np.zeros(N_samples)
    b_vec = np.linspace(b_val,b_val+0.4,N_samples)
    while i < N_samples :
        # probability of having a magnitude (GR law)
        x = np.random.exponential(1./b_vec[i])
        # modulation by the detection function
        if function == 'gumbel' :
            q = exp(-exp(-((x-mu)/sigma))) 
        elif function == 'full':
            q = 1
        elif function =='error':
            q=0.5+0.5*erf((x-mu)/(sigma*sqrt(2)))
        elif function == 'exp':
            if x > mu :
                q = 1 -( exp(-(x-mu)/sigma)*(1+(x-mu)/sigma +((x-mu)**2)/(2*sigma**2)))
            else : 
                q = 0
        else :
            print('We do not recognize the detection function ! Choose between : gumbel, exp, full, error')


        if np.random.rand() < q :
            i+=1
            M_distribution[i-1] = x
    return M_distribution


def etas(mu, alpha, barN, p, c, beta, tmax):
    """Draw all background earthquakes, and then the direct aftershocks of each consecutive earthquake.
   
    Il y a +sieurs taux de fond mu(1), mu(2), ..., identifiés comme parent=-1, -2, ...
    #Paramètres : (mu) nombre moyen de séismes / unité de temps en l'absence d'interactions
    # (alpha) en présence d'interactions, le nombre moyen de répliques pour un choc principal
    # de magnitude m est exp(alpha*m) [typiquement alpha~2], (barN) paramètre d'interaction
    # entre 0 (pas d'interations) et 1 (interaction maximale) [typiquement entre 0.5 et 0.9
    # pour la sismicité tectonique], (p) exposant de la loi d'Omori : nombre de répliques
    # décroit en temps (t) en 1/(t+c)^{-p} [il faut ici p>1, typiquement p~1.1 ou 1.2],
    # (c) temps de coupure de la loi d'Omori, (beta) valeur de b de Gutenberg-Richter * log(10),
    # (tmax) durée de la série en unité de temps"""
    
    #A = barN/(beta* c**(1-p))/((p-1)*(beta-alpha))

    A = barN / (beta * c**(1-p)/((p-1)*(beta-alpha)))
    tc = 0
    n = 0

    t, m, parent = [], [], []

    # background
    if mu:
        while tc < tmax:
            dt = -1 / np.sum(mu) * np.log(np.random.rand())
            tc = tc + dt
            if tc < tmax:
                n = n + 1
                t.append(tc)
                m.append(-1 / beta * np.log(np.random.rand()))
                I = np.where(np.cumsum(mu) / np.sum(mu) >= np.random.rand())[0]
                parent.append(-I[0])

    # cascade
    if A > 0:
        nc = 1
        #mmax = max(m) if m else 0
        while t :
            tc = 0
            while True:
                tmp = (tc + c)**(1 - p) + (p - 1) / (A * np.exp(alpha * m[nc - 1])) * np.log(np.random.rand())
                if tmp > 0:
                    dt = tmp**(1 / (1 - p)) - tc - c
                    tc = tc + dt
                    if tc + t[nc - 1] < tmax:
                        n = n + 1
                        t.append(tc + t[nc - 1])
                        new_m = -1 / beta * np.log(np.random.rand())
                        m.append(new_m)
                        # if new_m > mmax:
                        #     mmax = new_m
                        parent.append(nc)
                    else:
                        break
                else:
                    break

            perm = np.argsort(t)
            t = np.array(t)[perm].tolist()
            m = np.array(m)[perm].tolist()
            parent = np.array(parent)[perm].tolist()

            nc = nc + 1
            if nc > len(t):
                break

    if not t:
        t = m = parent = -1
    return t, m, parent



def ETAS_synthetic_dataset_vntot(N_vec, T_vec, b_vec, mu_vec, sigma_vec, function='gumbel'):
    """ Creates a synthetic catalog of N_sampels values of magnitudes
    following the conjoint detection distribution (mu, sigma) 
    and the GR distribution associated to the input b_val
    The detection function can be chosen from :
        - full
        - error
        - gumbel
        - exp
        - expdouble"""
    
    ntot = 0
    tmax = 0 
    M_distribution = []
    T_distribution = []
    mu_rate = [1]  # Example values for mu
    alpha = 2 # attention alpha < beta obligatoire sinon explosion
    barN = 0.95
    p = 1.1
    c = 0.0001
    
    for i in range (len(b_vec)):
        nsub = 0 
        while nsub < int(round(N_vec[i])) :
            beta = b_vec[i]
            mu = mu_vec[i]
            sigma = sigma_vec[i]
            nt = int(round(N_vec[i]))-nsub
            #print(nt)
            t_full, m_full, parent = etas(mu_rate, alpha, barN, p, c, beta, nt)
            print(t_full)
            try : 
                ncur = len(m_full)
            except TypeError :
                ncur = 0 
                x = m_full
                t = t_full
            else :
                for j in range(len(m_full)): 
                    x = m_full[j]
                    t = t_full[j]
                    # modulation by the detection function
                    if function == 'gumbel' :
                        q = exp(-exp(-((x-mu)/sigma))) 
                    elif function == 'full':
                        q = 1
                    elif function =='error':
                        q=0.5+0.5*erf((x-mu)/(sigma*sqrt(2)))
                    elif function == 'exp':
                        if x > mu :
                            q = 1 -( exp(-(x-mu)/sigma)*(1+(x-mu)/sigma +((x-mu)**2)/(2*sigma**2)))
                        else : 
                            q = 0
                    else :
                        print('We do not recognize the detection function ! Choose between : gumbel, exp, full, error')
            
            
                    if np.random.rand() < q :
                        ntot+=1
                        nsub+=1
                        t_q = t + tmax
                        
                        M_distribution.append( x)
                        T_distribution.append(t_q)
                        tmax= max(T_distribution)
    return  T_distribution, M_distribution

def ETAS_synthetic_dataset(N_vec, T_vec, b_vec, mu_vec, sigma_vec, function='error'):
    """ Creates a synthetic catalog of N_sampels values of magnitudes
    following the conjoint detection distribution (mu, sigma) 
    and the GR distribution associated to the input b_val
    The detection function can be chosen from :
        - full
        - error
        - gumbel
        - exp
        - expdouble"""
    
    ntot = 0
    tmax = 0 
    M_distribution = []
    T_distribution = []
    Liste_M = [[]] * len(b_vec)
    Liste_T = [[]] * len(b_vec)
    mu_rate = [1]  # Example values for mu
    alpha = 1.5 # attention alpha < beta obligatoire sinon explosion
    barN = 0.95
    p = 1.1
    c = 0.0001
    
    for i in range (len(b_vec)):
        print(i)
        # -- init --
        M_subset = []
        T_subset = []
        
        # Values b, mu sigma du subset
        t0 = T_vec[i]
        tvec = T_vec[i+1]-T_vec[i]
        beta = b_vec[i]
        mu = mu_vec[i]
        sigma = sigma_vec[i]
        nt = int(round(N_vec[i]))
        
        # ETAS subset 
        t_full, m_full, parent = etas(mu_rate, alpha, barN, p, c, beta, nt)
        tc_full = [t0+ ti * tvec/nt for ti in t_full]
        # Modulation by detection
        for j in range(len(m_full)): 
            x = m_full[j]
            t = tc_full[j]
            if function == 'gumbel' :
                q = exp(-exp(-((x-mu)/sigma))) 
            elif function == 'full':
                q = 1
            elif function =='error':
                q=0.5+0.5*erf((x-mu)/(sigma*sqrt(2)))
            elif function == 'exp':
                if x > mu :
                    q = 1 -( exp(-(x-mu)/sigma)*(1+(x-mu)/sigma +((x-mu)**2)/(2*sigma**2)))
                else : 
                    q = 0
            else :
                print('We do not recognize the detection function ! Choose between : gumbel, exp, full, error')
    
    
            if np.random.rand() < q :
                ntot+=1
                M_distribution.append(x)
                T_distribution.append(t)
                M_subset.append( x)
                T_subset.append(t)
                tmax= max(T_distribution)
                
            Liste_M[i] = M_subset
            Liste_T[i] = T_subset
            
            
    return T_distribution, M_distribution, Liste_M, Liste_T
from cmcrameri import cm


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
       for valeur in differences_magnitudes[m:m+w]:
           flag+=1
           if valeur > 0:
               flag_pos+=1 # Pour savoir sur combien d evenements on calcule effectivement la b-value
               differences_positives.append(valeur)
       windowed_mean = sum(np.ones(len(differences_positives)) * differences_positives) / len(differences_positives)
       
       Minimum_difference = min(differences_absolues[m:m+w])
       # d'apres la formule de Laplace (Van der Elst 2021)
       windowed_beta = 1/float(windowed_mean-Minimum_difference) 
       b_positive[m]=windowed_beta/np.log(10)  # On retrouve b_value = beta/2.3

    print("Valeur de b-value obtenue sur %s events sur fenetre de %s events" %(str(flag_pos),str(w)))
    times_b_positive = times[int(w/2):-int(w/2)]
    return times_b_positive, b_positive

def Beta_Laplace_function(mags, w):
    n_events =len(mags)
    difference = np.zeros(n_events-1)
    beta_positive= np.zeros(n_events-w)
    for i in range (0, n_events-1):
        difference [i] = abs(mags[i+1] - mags[i])
    flag=0
    for m in range(len(difference)-(w-1)):
       flag+=1
       windowed_mean = sum(np.ones(w) * difference[m:m+w]) / w 
       Minimum_difference = min(difference[m:m+w])
       windowed_beta = 1/float(windowed_mean-Minimum_difference)
       beta_positive[m]=windowed_beta/np.log(10)
    return beta_positive

def Beta_GR_function(mags,w, Mc):
    
    
    #difference[0] = 0
    flag=0
    mags_sup = []
    for mi in range(len(mags)):
       flag+=1
       if mags[mi] > Mc:
           mags_sup.append(mags[mi])


    mean = np.mean(mags_sup)
       
    beta_GR = 1/float(mean- Mc)
    sigma_GR =  beta_GR/sqrt(len(mags_sup))
    return beta_GR, sigma_GR

def cumulative_distribution(magnitudes1, min_mag = None, max_mag = None, interval = 0.05):

    """
    This function reads an earthquake catalogue file and creates the cumulative distribution
    Funtion arguments:
        infile: file containing earthquake catalogue
                Expected input file format: csv
                One header line
                Year of earthquake in 3rd column, magnitude in 6th column.
        min_mag: minimum magnitude for which data will be used - i.e. catalogue
                completeness
        max_mag: maximum magnitude used in bounded G-R curve. If not specified,
                defined as the maximum magnitude in the catlogue + 0.1 magnitude
                units.
        interval: Width of magnitude bins for generating cumulative histogram
                of earthquake recurrence. Default avlue 0.1 magnitude units.  
    
    """
        # If minimum magnitude is not specified, read all magnitudes
    if min_mag is not None:
        pass
    else:
        min_mag = -1.0
        

    # If minimum magnitude is not specified default value to minimum in catalogue
    if min_mag == -1.0:
        min_mag = min(magnitudes1)
            

    # If maximum magnitude is not specified default value to maximum in catalogue
    if max_mag is not None:
        pass
    else:
        max_mag = max(magnitudes1) + 0.1

    num_eq = len(magnitudes1)
    max_mag_bin = max(magnitudes1) + 0.15
    
    # Magnitude bins
    bins = np.arange(min_mag, max_mag_bin, interval)

    # Magnitude bins for plotting - we will re-arrange bins later
    plot_bins = np.arange(min_mag, max_mag, interval)

    ###########################################################################
    # Generate distribution
    ###########################################################################
    # Generate histogram
    hist1 = np.histogram(magnitudes1,bins=bins)
    # Reverse array order
    hist1 = hist1[0][::-1]
    bins = bins[::-1]

    # Calculate cumulative sum
    cum_hist1 = hist1.cumsum()
    # Ensure bins have the same length has the cumulative histogram.
    # Remove the upper bound for the highest interval.
    bins = bins[1:]
    return cum_hist1[::-1],bins[::-1]




if __name__ == "__main__" :
### First I build a synthetic distribution
    # =============================================================================
    # # Construction du catalogue ETAS
    # =============================================================================
    
    # # Example usage:
    # Transdimensional synth
    # T_vec = [0, 15, 20, 40, 60, 70, 88, 100]
    # N_vec = [1500,500,2000,2000,1000,1800,1200]  
    # b_vec =  np.array([0.77,1,1.2,0.55,0.7,0.7,0.75])*np.log(10)
    # mu_vec = np.array([-0.2,0.01,-0.1,-0.15,-0.15,0.1,0.12])
    # sigma_vec = np.array([0.32,0.32,0.32,0.32,0.32,0.32,0.32])
    
    
    
    # Simple - no transdimensional
    # T_vec = [0, 100]
    # N_vec = [1000]  
    # b_vec =  np.array([1.05])*np.log(10)
    # mu_vec = np.array([-0.14])
    # sigma_vec = np.array([0.19])
    

    # T_vec = [0, 100]
    # N_vec = [10000]  
    # b_vec =  np.array([0.97])*np.log(10)
    # mu_vec = np.array([-0.02])
    # sigma_vec = np.array([0.29])   
    # name= 'TB_6DISC_2109_GUMB_5000_PY_v2'
    # T_vec = [0, 15, 20, 40, 60, 70, 88, 100]
    # N_vec = [750,250,1000,1000,500,900,600]  
    # b_vec =  np.array([0.77, 1.4, 0.95, 0.6, 1.2, 0.75, 0.7])*np.log(10)
    # mu_vec = np.array([0.20, 0.9, 0.15, 0.1, 0.8, 0.10, 0.6])
    # sigma_vec = np.array([0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32])
    
    
    name= 'TB_6DISC_UNCLEAR_0403_ERROR_ETAS_highmu_d'
    Ntot = 8000
    T_vec = [0, 10, 20, 40, 60, 70, 88, 100]
    N_vec = [T_vec[1]*Ntot/100,(T_vec[2]-T_vec[1])*Ntot/100,(T_vec[3]-T_vec[2])*Ntot/100,(T_vec[4]-T_vec[3])*Ntot/100,(T_vec[5]-T_vec[4])*Ntot/100,(T_vec[6]-T_vec[5])*Ntot/100,(T_vec[7]-T_vec[6])*Ntot/100]  
    b_vec =  np.array([0.7,1.1, 1.0, 1.0, 0.8, 0.7, 0.65])*np.log(10)
    mu_vec = np.array([0.8, 0.8,1.5, 0.5, 0.75, 0.9, 0.9])
    sigma_vec = np.array([0.35, 0.3, 0.2, 0.15,0.3,0.4, 0.15])
    Mc = 1.8

    
    
    times,mags , lm, lt = ETAS_synthetic_dataset(N_vec, T_vec, b_vec, mu_vec, sigma_vec, function='error')

    print("Longueur du catalogue %f :"%(len(mags)))
        
    N_samples=len(mags)
    
    bins1 = np.arange(0,6,0.1)
    bins2 = np.arange(0,6,0.05)
    x_count = np.arange(0,N_samples, 1)
    

# =============================================================================
#     FIGURE CONSTRUCTION
# =============================================================================  
    color = cm.roma(np.linspace(0,1,len(b_vec)))
    fig, axs = plt.subplots(1,3,figsize=(18,5))
    plt.subplots_adjust(wspace=0.3, hspace=0.2)
    
    
    axs[0].set_xlim(0,100)
    axs[0].set_ylim(0,max(mags)+0.5)
    for kk in range(len(T_vec)-1):
        axs[0].vlines(T_vec[kk+1], 0,max(mags), linestyle='--', lw=2,color='grey')
        axs[0].scatter(lt[kk],lm[kk],color='grey', edgecolor=color[kk], s=[exp(magi)*1.5 for magi in lm[kk]], alpha=0.6, label = 'N=%s events'%(str(len(lm[kk]))))
    axs[0].set_ylabel('Magnitudes',fontsize=16)
    axs[0].set_xlabel('Times',fontsize=16)
    #axs[0].legend(loc = 'upper right', fontsize=12)



    
    for k in range(len(T_vec)-1):
        mags_subset = lm[k]
        n_subset = len(mags_subset)
        label = 'S'+str(k)+ ' - b = %s '%(str(round(b_vec[k]/np.log(10),2)))
        y, x = cumulative_distribution(mags_subset, min_mag = None, max_mag = None, interval = 0.05)

        lmag =[round(magi,1) for magi in mags_subset]
        axs[1].scatter(x, y, c=color[k], edgecolor='black',linewidth=0.7, label=label)
        # Création d'un dictionnaire de magnitude:frequence pour chaque décimale (peut aussi être réalisé avec pd.Dateframe())
        X = {x:lmag.count(x) for x in lmag}  ## Dictionnaire of magnitude:frequency for MLv - Seiscomp
        lmag_round, fmag_round = X.keys(), X.values()
    
    axs[1].set_yscale('log')
    axs[1].set_xlabel('Magnitudes',fontsize=16)
    axs[1].set_ylabel('Log(N)',fontsize=16)
    axs[1].grid(lw=0.5)
    axs[1].set_xlim(0,6)
    axs[1].legend(loc='upper right', fontsize=12)



    w=200
    Mc= 1.8
    t_pos, b_pos  = B_Positive_function( times,mags, w)
    b_lap  = Beta_Laplace_function( mags, w)
    b_aki , s_aki = Beta_GR_function(mags,w, Mc) 
    axs[2].plot(times[int(w/2):-int(w/2)],b_pos, color='navy', linestyle='--',label = 'b Positive| w='+str(w))
    axs[2].hlines(b_aki/np.log(10),0,100, color='black',linestyles='--',lw=2,label = 'b Aki (full) = %s+/-%s'%(str(round(b_aki/np.log(10),2)), round(s_aki/np.log(10),2)))
    axs[2].hlines((b_aki+s_aki)/np.log(10),0,100, color='grey',linestyles='--',lw=1,label = 'B Aki (subsets)')
    axs[2].hlines((b_aki-s_aki)/np.log(10),0,100, color='grey',linestyles='--',lw=1)
        
    for k in range(len(T_vec)-1):
        mags_subset = lm[k]
        axs[2].hlines(b_vec[k]/np.log(10),T_vec[k],T_vec[k+1], color='red',linestyles='--',lw=3)
        axs[2].vlines(T_vec[k+1], 0.4,1.5, linestyle='--', lw=2,color='grey')
        Mc = mu_vec[k]+sigma_vec[k]
        b_GR , s_GR = Beta_GR_function(mags_subset, w, Mc)
        axs[2].hlines(b_GR/np.log(10),T_vec[k],T_vec[k+1], color='goldenrod',linestyles='--',lw=2)
        axs[2].hlines((b_GR+s_GR)/np.log(10),T_vec[k],T_vec[k+1], color='wheat',linestyles='--',lw=1)
        axs[2].hlines((b_GR-s_GR)/np.log(10),T_vec[k],T_vec[k+1], color='wheat',linestyles='--',lw=1)
    axs[2].legend(loc='upper right', fontsize= 12)
        
        
    axs[2].set_xlim(0,100)
    axs[2].set_ylim(0.4,1.5)
    axs[2].set_xlabel('b-value',fontsize=16)
    axs[2].set_ylabel('Log(N)',fontsize=16)
    
    # axs[3].set_xlim(0,100)
    # axs[3].set_ylim(0.1,1.8)
    # for kk in range(len(T_vec)):
    #     axs[3].vlines(T_vec[kk], 0.1,1.8, linestyle='--', lw=2,color='grey', label='Fixed discontinuities')
    # for di in range(len(T_vec)-1):
    #     axs[3].hlines(mu_vec[di],T_vec[di],T_vec[di+1], color='blue',linestyles='--',lw=3)
    # axs3b = axs[3].twinx()
    # axs3b.plot(times[int(w/2):-int(w/2)],b_pos, color='green', linestyle='--',label = 'B Positive| w='+str(w))

    fig.suptitle('Synthetic seismic dataset',fontsize=17)
    fig.savefig(name+'_Fig2_Synthetic_dataset.pdf', format = 'pdf')
    
    plt.show()
    
    
    
    
# =============================================================================
#     Writting catalog
# =============================================================================     
    import csv
    
    with open(name+'.in', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(times,mags))
        
    with open(name+'_inputs.txt', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(T_vec,N_vec, b_vec, mu_vec, sigma_vec))
