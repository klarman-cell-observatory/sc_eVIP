
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import sys 
import seaborn as sns
import os

def get_conf_interval(data,ci_size=0.95):
    from scipy.stats import t, sem
    ci1,ci2=t.interval(alpha=ci_size, df=len(data)-1, 
               loc=np.mean(data), scale=sem(data)) 
    return(ci1,ci2)

def empirical_pvalues(values_df,empirical_dist):

    values=list(values_df)
    from scipy import stats
    ps=[]
    for i in range(len(values)):
        p=(100.0-stats.percentileofscore(empirical_dist,values[i]))/100.0
        ps.append(p)
    return(np.array(ps))

def qvalue(ps,mini=1e-5,thresh=0.5):

    total=ps.shape[0]
    above=(ps>thresh).sum()
    exp_fp=int(above*(1.0/thresh))
    pi_0=1.0*exp_fp/total

    qs=[]
    for p_idx in range(ps.shape[0]):
        p=ps[p_idx]
        R=(ps<=p).sum()
        q=max(1.0*p*total*pi_0/R,mini)
        qs.append(q)
    return(qs)


def compare_2_groups_df(cellsxvalues,
                     group1_cells,group2_cells,
                     method):
    
    g1_data=cellsxvalues.loc[group1_cells,:]
    g2_data=cellsxvalues.loc[group2_cells,:]

    #T2 Hotelling
    #============
    if method=='HotellingT2':
        import spm1d
        
        T2=spm1d.stats.hotellings2(g1_data,g2_data)
        value=T2.z #the T2 statistic

    #bulk analysis
    #==============
    if 'avg' in method:
        #get the averaged data
        g1_bulk=np.mean(np.array(g1_data),axis=0)
        g2_bulk=np.mean(np.array(g2_data),axis=0)

    if method=='avg.pearson':
        from scipy.stats import pearsonr
        value=1-pearsonr(g1_bulk,g2_bulk)[0]

    if method=='avg.spearman':
        from scipy.stats import spearmanr
        value=1-spearmanr(g1_bulk,g2_bulk)[0]
        
    if method=='avg.L1':
        from scipy.spatial.distance import cityblock
        value=cityblock(g1_bulk,g2_bulk)

    return(value)


def df_to_CI(res_df,score_cols,controls,control_col='group2',test_col='group1',ci_size=0.95):
    
    res_df.index=range(res_df.shape[0])
    groups=list(set(res_df[test_col]))
    res_cols=[]
    for score_col in score_cols:
        res_cols.append(score_col+'.mean')
        res_cols.append(score_col+'.ci.'+str(ci_size)+'.low')
        res_cols.append(score_col+'.ci.'+str(ci_size)+'.high')
        
    res_ci=pd.DataFrame(index=groups,columns=res_cols)
    for score_col in score_cols:
        for group in groups:
            rows=list(set(list(res_df.loc[res_df[test_col]==group,:].index)).intersection(set(list(res_df.loc[res_df[control_col].isin(controls),:].index))))
            group_scores=res_df.loc[rows,score_col]
            conf_intervals=get_conf_interval(group_scores,ci_size=ci_size)
            res_ci.loc[group,score_col+'.ci.'+str(ci_size)+'.low']=conf_intervals[0]
            res_ci.loc[group,score_col+'.ci.'+str(ci_size)+'.high']=conf_intervals[1]
            res_ci.loc[group,score_col+'.mean']=np.mean(group_scores)
    return(res_ci)

def get_empirical_q(res_full,res_ci,
                    controls,
                    score_cols,
                    control_col='group2',test_col='group1'):

    for score_col in score_cols:
        
        rows1=list(res_full.loc[res_full[control_col].isin(controls),:].index)
        rows2=list(res_full.loc[res_full[test_col].isin(controls),:].index)
        rows=list(set(rows1).intersection(set(rows2)))
        
        emp_scores=np.array(res_full.loc[rows,score_col].astype(float)).flatten()
        emp_scores = np.array(emp_scores[~np.isnan(emp_scores)]).flatten()
        
        #q-values
        emp_p=empirical_pvalues(res_ci[score_col+'.mean'],emp_scores)
        emp_p_ctrl=empirical_pvalues(emp_scores,emp_scores)
        emp_p_combo=list(emp_p)
        for p in emp_p_ctrl:
            emp_p_combo.append(p)
        res_ci[score_col+'.q']=qvalue(np.array(emp_p_combo))[:res_ci.shape[0]]
    return(res_ci)

def display_progress(cur_num,max_num):
    sys.stdout.write('\r'+str(int(100*(1.0*cur_num/max_num)))+' %')
    sys.stdout.flush()

def compare_groups_with_reference(data_df,
                                  labels_df,
                                 groups,
                                  controls,
                                 methods=['HotellingT2','avg.pearson'],
                                 n_bootstrap_controls=0,
                                 rng=np.random.RandomState(1234)):
    
    import copy
    groups_found=list(set(labels_df['label']))
    controls_found=list(set(groups_found).intersection(set(controls)))
    if len(controls_found)>0:
        print('Found '+str(len(controls_found))+'/'+str(len(controls_found))+' controls')
    else:
        print('ERROR: No controls found.')
        return()
    
    #setup result data
    res_col=['group1','group2']
    for c in methods:
        res_col.append(c)
    res=pd.DataFrame(columns=res_col)
    
    #score groups against all controls, and all controls against all controls
    score_mat={}
    for method in methods:
        score_mat[method]=pd.DataFrame(index=controls_found,columns=groups_found)
        
    #get group cells 
    group_cells_d={}
    for group in groups:
        if group not in groups_found:
            print('WARNING: group "'+group+'" is not in label data frame. Skipping')
            return
        group_cells=list(labels_df.loc[labels_df['label']==group].index)
        group_cells_d[group]=group_cells
        
    #get control cells
    control_cells_d={}
    controls_found_plus_bootstraps=[]
    for control in controls_found:
        controls_found_plus_bootstraps.append(control)
        control_cells=list(labels_df.loc[labels_df['label']==control].index)
        control_cells_d[control]=control_cells
        #if bootstrapping, sample with replacement from the control cells
        for b in range(n_bootstrap_controls):
            control_cells_d[control+'boot'+str(b)]=rng.choice(control_cells,size=len(control_cells))
            controls_found_plus_bootstraps.append(control+'boot'+str(b))
        
    #score groups against all controls
    counter=0
    for group in groups:
        group_cells=group_cells_d[group]
        counter+=1
        display_progress(counter,len(groups_found))
        
        for control in controls_found_plus_bootstraps:
            if control==group:
                continue
            #group is a test group and control is a bootstrapped one, skip test
            if group not in controls_found and control not in controls_found:
                continue
            control_cells=control_cells_d[control]
            res_here=pd.DataFrame({'group1':[group],
                                  'group2':[control]},index=[group+'.'+control])
            
            for method in methods:
                method_score=compare_2_groups_df(data_df.loc[list(set(group_cells).union(set(control_cells))),:],
                     group_cells,control_cells,
                     method)
                
                res_here[method]=method_score    
        
            res=pd.concat([res,res_here],axis=0,sort=True)
    res.index=list(res['group2'])
    
    #make confidence intervals
    res_ci=df_to_CI(res,methods,
                    controls_found_plus_bootstraps,
                    control_col='group2',
                    test_col='group1',
                    ci_size=0.95)
    
    #get qvalues from empirical p-values
    res_ci=get_empirical_q(res,res_ci,
                    controls=controls_found_plus_bootstraps,
                    control_col='group2',test_col='group1',
                    score_cols=methods)
    
    return(res_ci)
        
