import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def QuantBond(df, char='num_pop', epoch=20):
    df_sum = np.zeros(shape=np.array(df[df.Epoch==1])[:, :, np.newaxis].shape)
    for i in np.arange(epoch)+1:
        df_array = np.array(df[df.Epoch==i])[:, :, np.newaxis]
        df_sum = np.concatenate([df_sum, df_array], 2)
    df_sum = df_sum[:, :, 1:]
    df_temp = df_sum[:, 4:6, 0]
    df_temp_df = pd.DataFrame(df_temp, columns=['Finput', "CVg"])
    df_min = np.concatenate([np.min(df_sum[:, 0:4, :], axis=2), df_temp], 1)
    df_mean = np.concatenate([np.mean(df_sum[:, 0:4, :], axis=2), df_temp], 1)
    df_max = np.concatenate([np.max(df_sum[:, 0:4, :], axis=2), df_temp], 1)
    
    df_min = pd.DataFrame(df_min, columns=df.columns[:-1])
    df_mean= pd.DataFrame(df_mean, columns=df.columns[:-1])
    df_max = pd.DataFrame(df_max, columns=df.columns[:-1])
    return(df_min, df_mean, df_max)
  
  
def gradient_start2end (df_list, Vtrait='SumBiomass', scenario='evo'):
  Age_=np.array([10, 4, 8, 6, 7, 9])
  Finput_list = df_list[0]['Finput'].unique()
  CVg_list = df_list[0]['CVg'].unique()
  Gradient_List = []
  for spe_idx, spe in enumerate(df_list):
      spe_list = []
      for i, j in enumerate(Finput_list):
          low_ = 50-Age_[spe_idx]-10
          upper_ = 100-Age_[spe_idx]
          start_ref = spe[(spe.Finput==j)&(spe.CVg==CVg_list[0])].iloc[low_:upper_,:][Vtrait].mean()
          if scenario == 'evo':
              end_ref = spe[(spe.Finput==j)&(spe.CVg==CVg_list[1])][Vtrait][240:].mean()
          else:
              end_ref = spe[(spe.Finput==j)&(spe.CVg==CVg_list[0])][Vtrait][240:].mean()
          recover_rate = end_ref/start_ref
          spe_list.append(recover_rate)
      Gradient_List.append(spe_list)
  return(Gradient_List)


def gradient_half2end (df_end_list, df_half_list, Vtrait = 'SumBiomass', scenario = 'evo'):
    Finput_list = df_end_list[0]['Finput'].unique()
    CVg_list = df_end_list[0]['CVg'].unique()
    Gradient_end_List = []
    Gradient_half_List = []
    for spe_idx, spe in enumerate(zip(df_end_list, df_half_list)):
        spe_end_list = []
        spe_half_list = []
        for i, j in enumerate(Finput_list):
            if scenario == 'evo':
                end_ref = spe[0][(spe[0].Finput==j)&(spe[0].CVg==CVg_list[1])][Vtrait].mean()
                half_ref = spe[1][(spe[1].Finput==j)&(spe[1].CVg==CVg_list[1])][Vtrait].mean()
            else:
                end_ref = spe[0][(spe[0].Finput==j)&(spe[0].CVg==CVg_list[0])][Vtrait].mean()
                half_ref = spe[1][(spe[1].Finput==j)&(spe[1].CVg==CVg_list[0])][Vtrait].mean()
            spe_end_list.append(end_ref)
            spe_half_list.append(half_ref)
        Gradient_end_List.append(spe_end_list)
        Gradient_half_List.append(spe_half_list)
    return(Gradient_end_List, Gradient_half_List)


def do_stuff(cell): #show up the plots
    ax = plt.subplot(cell)
    ax.plot()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    
def Matist(Linf, K, t0, Lmin, Lmax):
    L1 = (Lmin+Lmax)/2
    firstMat = t0+1/K*np.log(Linf/(Linf-L1))
    return firstMat
  
  
def logis(a, b, linf):
    l = np.arange(linf)
    S_l = np.exp(a+b*l)/(1+np.exp(a+b*l))
    return (l, S_l)

def logis_es(lmax, alpha, linf):
    rand_alpha = alpha*np.random.normal(1, 0.1, 1)
    rand_lmax = lmax*np.random.normal(1, 0.1, 1)
    l1 = np.arange(linf)
    S_l1 = 1/(1+np.exp(rand_alpha*(rand_lmax-l1)))
    return (l1, S_l1)

def Logit(linfty, l50, SR):
    l = np.arange(linfty)
    r = (np.exp(np.log(9)/SR*(l-l50)))/(1.0+np.exp(np.log(9)/SR*(l-l50)))
    return (l*10, r)
  
  
def Sens(Fish_list, ax_i, para, stat, ax):
    Age_ = [10, 4, 8, 6, 7, 9]
    Fish_name = ['JSM', 'JA', 'PCM', 'SYC', 'LH', 'BS']
    fontsize=16
    for spe_idx, spe in enumerate(Fish_list):  
        Para_list = spe[para].unique()
        spe_list = []
        for para_idx, para_value in enumerate(Para_list):
            if para == "Finput":  ## Finput a.k.a MIHR 
                spe_biomass = spe[spe[para]==para_value][110:160][stat]
            else:
                spe_biomass = spe[spe[para]==para_value][50:100-Age_[spe_idx]][stat]
            spe_list.append(spe_biomass)
        ax[ax_i, spe_idx].boxplot(spe_list, labels=Para_list)
        ax[ax_i, spe_idx].set_xticklabels(Para_list, fontsize=16, rotation=90)
        ax[ax_i, spe_idx].ticklabel_format(style='sci', scilimits=(-2, 2), axis='y')
        if ax_i==0:
            ax[ax_i, spe_idx].set_title(Fish_name[spe_idx], fontsize=fontsize+13, fontweight='bold')
            if spe_idx==0:
                ax[ax_i, spe_idx].set_ylabel('Biomass', fontsize=fontsize+13)
