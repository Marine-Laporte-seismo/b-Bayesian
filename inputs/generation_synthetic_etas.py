import numpy as np
import matplotlib.pyplot as plt
import time
from math import exp, sqrt, erf
from datetime import date
from scipy.special import ndtr, log_ndtr
from scipy.optimize import minimize
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
    # (l'ordre de traitement des evenements n'a pas d'importance dans un processus de
    #  branchement : on ne trie qu'une seule fois a la fin, ce qui evite un tri par
    #  evenement, tres lent pour les gros catalogues)
    if A > 0:
        nc = 1
        while nc <= len(t):
            tc = 0
            while True:
                tmp = (tc + c)**(1 - p) + (p - 1) / (A * np.exp(alpha * m[nc - 1])) * np.log(np.random.rand())
                if tmp > 0:
                    dt = tmp**(1 / (1 - p)) - tc - c
                    tc = tc + dt
                    if tc + t[nc - 1] < tmax:
                        n = n + 1
                        t.append(tc + t[nc - 1])
                        m.append(-1 / beta * np.log(np.random.rand()))
                        parent.append(nc)
                    else:
                        break
                else:
                    break
            nc = nc + 1

    perm = np.argsort(t)
    t = np.array(t)[perm].tolist()
    m = np.array(m)[perm].tolist()
    parent = np.array(parent)[perm].tolist()
    return t, m, parent



def ETAS_synthetic_dataset(N_vec, T_vec, b_vec, mu_vec, sigma_vec, function='error',
                           mu_rate=[1], alpha=1.5, barN=0.95, p=1.1, c=0.0001):
    """ Creates a synthetic catalog of N_sampels values of magnitudes
    following the conjoint detection distribution (mu, sigma)
    and the GR distribution associated to the input b_val
    The detection function can be chosen from :
        - full
        - error
        - gumbel
        - exp
        - expdouble
    N_vec : duree de chaque segment en nombre moyen de seismes de fond
            (le nombre final d'evenements detectes est different)
    b_vec : valeurs de beta = b*ln(10)"""

    ntot = 0
    tmax = 0
    M_distribution = []
    T_distribution = []
    Liste_M = [[] for _ in range(len(b_vec))]
    Liste_T = [[] for _ in range(len(b_vec))]
    Alpha_used = []
    N_raw = []

    for i in range (len(b_vec)):
        print('Segment', i+1)
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

        # ETAS exige alpha < beta, sinon A <= 0 et AUCUNE replique n'est generee
        alpha_i = alpha
        if alpha >= 0.9 * beta:
            alpha_i = 0.9 * beta
            print('  ATTENTION : alpha=%.2f trop proche de beta=%.2f -> alpha=%.2f pour ce segment'
                  % (alpha, beta, alpha_i))
        Alpha_used.append(alpha_i)

        # ETAS subset
        t_full, m_full, parent = etas(mu_rate, alpha_i, barN, p, c, beta, nt)
        N_raw.append(len(m_full))
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
        print('  %d evenements generes, %d detectes' % (len(m_full), len(M_subset)))


    return T_distribution, M_distribution, Liste_M, Liste_T, Alpha_used, N_raw


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


def neg_loglike(params, m, Mmin):
    """- log-vraisemblance GR x detection 'error' (meme formule que le Fortran)."""
    b, mu, sigma = params
    beta = b * np.log(10)
    K = ndtr((Mmin - mu) / sigma) + np.exp(0.5 * (beta * sigma) ** 2 - beta * (mu - Mmin)) \
        * ndtr(-(Mmin + beta * sigma ** 2 - mu) / sigma)
    ll = log_ndtr((m - mu) / sigma).sum() + len(m) * (np.log(beta) - np.log(K)) \
        - beta * (m - Mmin).sum()
    return -ll


def fit_segment(m, Mmin):
    """Maximum de vraisemblance (b, mu, sigma) sur un segment, dans les bornes des priors Fortran."""
    x0 = [1.0, np.median(m) - 0.3, 0.2]
    bounds = [(0.5 * 2.3 / np.log(10), 2.0 * 2.3 / np.log(10)), (0.0, 4.0), (0.0001, 0.5)]
    res = minimize(neg_loglike, x0, args=(m, Mmin), method='L-BFGS-B', bounds=bounds)
    return res.x




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


    seed = 11                     # graine aleatoire (meme graine = meme catalogue)
    np.random.seed(seed)
    function = 'error'             # 'error' = fonction de detection du Fortran

    # Parametres ETAS
    mu_rate = [1]                  # taux de fond (evenements / unite de temps ETAS)
    alpha = 1.5                   # attention alpha < beta obligatoire sinon explosion
    barN = 0.8                    # taux de branchement (0.8 = cascades plus raisonnables)
    p = 1.1
    c = 0.0001

    name= 'TB_6DISC_UNCLEAR_2710_ERROR_ETAS_erf'
    Ntot = 6000
    T_vec = [0, 10, 20, 40, 60, 70, 88, 100]
    N_vec = [(T_vec[k+1]-T_vec[k])*Ntot/100 for k in range(len(T_vec)-1)]
    b_vec =  np.array([0.7,1.1, 1.0, 1.0, 0.9, 0.7, 0.65])*np.log(10)
    mu_vec = np.array([0.8, 0.8,1.5, 0.5, 0.7, 0.9, 0.9])
    sigma_vec = np.array([0.35, 0.3, 0.2, 0.15,0.3,0.4, 0.15])
    Mc = 1.8



    times,mags , lm, lt, alpha_used, N_raw = ETAS_synthetic_dataset(
        N_vec, T_vec, b_vec, mu_vec, sigma_vec, function=function,
        mu_rate=mu_rate, alpha=alpha, barN=barN, p=p, c=c)

    print("Longueur du catalogue : %d"%(len(mags)))

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


    axs[0].set_xlim(T_vec[0],T_vec[-1])
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
    axs[2].hlines(b_aki/np.log(10),T_vec[0],T_vec[-1], color='black',linestyles='--',lw=2,label = 'b Aki (full) = %s+/-%s'%(str(round(b_aki/np.log(10),2)), round(s_aki/np.log(10),2)))
    axs[2].hlines((b_aki+s_aki)/np.log(10),T_vec[0],T_vec[-1], color='grey',linestyles='--',lw=1,label = 'B Aki (subsets)')
    axs[2].hlines((b_aki-s_aki)/np.log(10),T_vec[0],T_vec[-1], color='grey',linestyles='--',lw=1)

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


    axs[2].set_xlim(T_vec[0],T_vec[-1])
    axs[2].set_ylim(0.4,1.5)
    axs[2].set_xlabel('Times',fontsize=16)
    axs[2].set_ylabel('b-value',fontsize=16)

    fig.suptitle('Synthetic seismic dataset',fontsize=17)
    fig.savefig(name+'_Fig2_Synthetic_dataset.pdf', format = 'pdf')

    plt.show()




# =============================================================================
#     Writting catalog (2 colonnes : temps, magnitude, trie en temps -> Fortran)
# =============================================================================
    times = np.asarray(times)
    mags = np.asarray(mags)
    order = np.argsort(times)
    times, mags = times[order], mags[order]
    np.savetxt(name+'.in', np.column_stack([times, mags]), fmt='%.6f  %.4f')

# =============================================================================
#     Writting parameters file
# =============================================================================
    Mmin = mags.min()
    fits = [fit_segment(np.asarray(lm[k]), Mmin) if len(lm[k]) > 10 else [np.nan]*3
            for k in range(len(b_vec))]

    with open(name+'_inputs.txt', 'w') as f:
        f.write("# Test synthetique ETAS + detection - RJ-MCMC\n")
        f.write("# Genere le %s\n" % date.today().isoformat())
        f.write("catalogue          = %s.in\n" % name)
        f.write("seed               = %d\n" % seed)
        f.write("fonction_detection = %s\n" % function)
        f.write("ETAS : mu_rate=%s  alpha=%s  barN=%s  p=%s  c=%s\n" % (mu_rate, alpha, barN, p, c))
        f.write("nb_segments        = %d\n" % (len(b_vec)))
        f.write("N_total (detectes) = %d\n" % len(mags))
        f.write("mag_min / mag_max  = %.4f / %.4f\n" % (mags.min(), mags.max()))
        f.write("temps min / max    = %.4f / %.4f\n\n" % (times.min(), times.max()))

        f.write("# Parametres vrais par segment\n")
        f.write("# N_fond = duree ETAS (nb moyen de seismes de fond), N_gen = evenements generes, N_det = detectes\n")
        f.write("# %3s %8s %8s %8s %7s %7s %7s %7s %7s %7s %8s %8s\n" %
                ('seg', 't_debut', 't_fin', 'N_fond', 'N_gen', 'N_det', 'b', 'beta', 'mu', 'sigma', 'alpha', 'mag_max'))
        for k in range(len(b_vec)):
            mmax = max(lm[k]) if len(lm[k]) > 0 else np.nan
            f.write("  %3d %8.3f %8.3f %8.0f %7d %7d %7.3f %7.4f %7.3f %7.3f %8.3f %8.3f\n" %
                    (k+1, T_vec[k], T_vec[k+1], N_vec[k], N_raw[k], len(lm[k]),
                     b_vec[k]/np.log(10), b_vec[k], mu_vec[k], sigma_vec[k], alpha_used[k], mmax))

        f.write("\n# Verification : maximum de vraisemblance par segment (bornes vraies connues)\n")
        f.write("# (valable pour function='error' ; Mmin = minimum du catalogue entier, comme le Fortran)\n")
        f.write("# %3s %7s %7s %8s %7s %8s %7s\n" % ('seg', 'b_vrai', 'b_MLE', 'mu_vrai', 'mu_MLE', 'sig_vrai', 'sig_MLE'))
        for k in range(len(b_vec)):
            f.write("  %3d %7.3f %7.3f %8.3f %7.3f %8.3f %7.3f\n" %
                    (k+1, b_vec[k]/np.log(10), fits[k][0], mu_vec[k], fits[k][1], sigma_vec[k], fits[k][2]))

        f.write("\n# --- A copier dans plot_results_RJMCMC.py (avec synthetics='True') ---\n")
        f.write("T_vec = %s\n" % list(T_vec))
        f.write("N_vec = %s\n" % [len(lm[k]) for k in range(len(b_vec))])
        f.write("b_vec = np.array(%s)*np.log(10)\n" % [float(round(float(b)/np.log(10), 4)) for b in b_vec])
        f.write("mu_vec = np.array(%s)\n" % [float(x) for x in mu_vec])
        f.write("sigma_vec = np.array(%s)\n" % [float(x) for x in sigma_vec])
        f.write("\n# --- Dans le Fortran ---\n")
        f.write("# Magnitude_file = '../inputs/%s.in'\n" % name)

    print("Catalogue  : %s.in" % name)
    print("Parametres : %s_inputs.txt" % name)