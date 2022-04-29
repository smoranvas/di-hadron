import numpy as np
import pandas as pd
def getRatio(df_A,df_D,df_trigger_A,df_trigger_D, variable='h2_z',trig_cut = 'h1_z>0.5', pair_cut='',minz=0.05,
             maxz=0.5,nbins=9, applyweight=False):
    print ('Print Trigger Cut ' ,trig_cut)
    print ('Total Cut ', trig_cut + pair_cut)
    norm_A = df_trigger_A.query(trig_cut).shape[0] 
    norm_D = df_trigger_D.query(trig_cut).shape[0]
    bins= np.linspace(minz,maxz,nbins)
    y_A, x_conditional = np.histogram(df_A.query(trig_cut+pair_cut)[variable],bins=bins)
    y_D, x_conditional = np.histogram(df_D.query(trig_cut+pair_cut)[variable],bins=bins)
    erry_A, x_conditional = np.histogram(df_A.query(trig_cut+pair_cut)[variable],bins=bins)
    erry_D, x_conditional = np.histogram(df_D.query(trig_cut+pair_cut)[variable],bins=bins)
    x_conditional = (x_conditional[1:] + x_conditional[:-1])/2.0
    err_A = np.true_divide(np.sqrt(erry_A),y_A)
    err_D = np.true_divide(np.sqrt(erry_D),y_D)
    y_A = np.true_divide(y_A,norm_A)
    y_D = np.true_divide(y_D,norm_D)
    ratio_conditional = np.true_divide(y_A,y_D)
    error_conditional = np.multiply(ratio_conditional, np.sqrt(np.power(err_A,2.0) + np.power(err_D,2.0)))
    return ratio_conditional,error_conditional,x_conditional

def getRatio2d(df_A,df_D,df_trigger_A,df_trigger_D,trig_cut = 'h1_z>0.5', pair_cut='', varx="h2_z",bins=None, vary="h1_z", applyweight=False, verbose=False):
    if verbose:
        print ('Print Trigger Cut ' ,trig_cut)
        print ('Total Cut ', trig_cut + pair_cut)
    if varx != "h1_z" and vary != "h1_z":
        norm_A = df_trigger_A.query(trig_cut).shape[0] 
        norm_D = df_trigger_D.query(trig_cut).shape[0]
        #when binning in z1, the normalizations, N1pi(z1), are differential as well.  
    elif varx == "h1_z":
        norm_A = np.histogram(df_trigger_A.query(trig_cut).eval("h1_z"), bins=bins[0])
        norm_D = np.histogram(df_trigger_D.query(trig_cut).eval("h1_z"), bins=bins[0])
        norm_A = np.tile(norm_A,(len(bins[1])-1,1)).transpose()
        norm_D = np.tile(norm_D,(len(bins[1])-1,1)).transpose()
    elif vary == "h1_z":
        norm_A = np.histogram(df_trigger_A.query(trig_cut).eval("h1_z"), bins=bins[1])[0]
        norm_D = np.histogram(df_trigger_D.query(trig_cut).eval("h1_z"), bins=bins[1])[0]
        norm_A = np.tile(norm_A,(len(bins[0])-1,1))
        norm_D = np.tile(norm_D,(len(bins[0])-1,1))
    df_A = df_A.query(trig_cut+pair_cut)
    df_D = df_D.query(trig_cut+pair_cut)
    
    n_A, x,y = np.histogram2d(df_A.eval(varx), df_A.eval(vary),bins)
    n_D, x,y= np.histogram2d(df_D.eval(varx), df_D.eval(vary),bins)
    err_A = np.sqrt(n_A)
    err_D = np.sqrt(n_A)
    #x = (x[1:] + x[:-1])/2.0
    #y = (y[1:] + y[:-1])/2.0
    z_A = np.true_divide(n_A,norm_A)
    z_D = np.true_divide(n_D,norm_D)
    ratio_conditional = np.true_divide(z_A,z_D)
    #error_conditional = np.multiply(ratio_conditional, np.sqrt(np.power(err_A,2.0) + np.power(err_D,2.0)))
    error_conditional=ratio_conditional*np.sqrt(1/n_A+1/n_D)
    return ratio_conditional,error_conditional,x, y

def getMultiplicity(df,df_trigger, variable='h2_z',trig_cut = 'h1_z>0.5', pair_cut='',minz=0.05,maxz=0.5,nbins=9, applyweight=False):
    norm = df_trigger.query(trig_cut).shape[0]
    bins= np.linspace(minz,maxz,nbins)
    y, x = np.histogram(df.query(trig_cut+pair_cut)[variable],bins=bins)
    erry, x = np.histogram(df.query(trig_cut+pair_cut)[variable],bins=bins)
    x = (x[1:] + x[:-1])/2.0
    err = np.true_divide(np.sqrt(erry),y)
    y = np.true_divide(y,norm)
    return y,err,x


def applyCut(inputDataframe, cut, text=None):
    nbeforecut = inputDataframe.shape[0]
    cutDataframe = None
    if nbeforecut>0:
        cutDataframe = inputDataframe.query(cut)
        if text:
            print (text, cutDataframe.shape[0], ' (%2.2f '%(100.0*cutDataframe.shape[0]/nbeforecut), '%)')
    return cutDataframe

# defining the limits and number of bins for each variable.
dpionMassBinsEdges=11
maxzmass=1.7
minzmass=0.3

dpiondphiBinsEdges=9
maxdphi=3.14
mindphi=0.0

dpionz2BinsEdges=9
maxz2=0.45
minz2=0.05


p_thr=2.7
#Nphe_thr=15
Nphe_thr=5
Nphe_h1_cut='h1_Nphe> (%d*( (h1_z*nu*h1_z*nu-0.13957*0.13957)>(%f*%f) and h1_pid==211 ) -10000* (not ( (h1_z*nu*h1_z*nu-0.13957*0.13957)>(%f*%f) and h1_pid==211 ) ) )'%(Nphe_thr,p_thr,p_thr,p_thr,p_thr)
Nphe_h2_cut='h2_Nphe> (%d*( (h2_p*nu*h2_z*nu-0.13957*0.13957)>(%f*%f) and h2_pid==211 ) -10000* (not ( (h2_z*nu*h2_z*nu-0.13957*0.13957)>(%f*%f) and h2_pid==211 ) ) )'%(Nphe_thr,p_thr,p_thr,p_thr,p_thr)
Chi2CC_h1_cut='h1_Chi2CC < (0.08726*( (h1_z*nu*h1_z*nu-0.13957*0.13957)>(%f*%f) and h1_pid==211 ) +10000* (not ( (h1_z*nu*h1_z*nu-0.13957*0.13957)>(%f*%f) and h1_pid==211 ) ) )'%(p_thr,p_thr,p_thr,p_thr)
Chi2CC_h2_cut='h2_Chi2CC < (0.08726*( (h2_z*nu*h2_z*nu-0.13957*0.13957)>(%f*%f) and h2_pid==211 ) +10000* (not ( (h2_z*nu*h2_z*nu-0.13957*0.13957)>(%f*%f) and h2_pid==211 ) ) )'%(p_thr,p_thr,p_thr,p_thr)
StatCC_h1_cut='h1_StatCC>((( (h1_z*nu*h1_z*nu-0.13957*0.13957)>(%f*%f) and h1_pid==211 ) -1 ))'%(p_thr,p_thr)
StatCC_h2_cut='h2_StatCC>((( (h2_z*nu*h2_z*nu-0.13957*0.13957)>(%f*%f) and h2_pid==211 ) -1 ))'%(p_thr,p_thr)

trigger_cut_nom   ='h1_z>0.5 and abs(h1_deltaZ)<3.0 and TargType!=0 and SampFracEl25==1 and h1_FidCutPiPlus==1 and %s and %s and %s  '%(Nphe_h1_cut, Chi2CC_h1_cut,StatCC_h1_cut)
pair_cut_nom  ='h1_z>0.5 and abs(h2_deltaZ)<3.0 and abs(h1_deltaZ)<3.0 and TargType!=0 and SampFracEl25==1 and %s and %s and %s and %s and %s and %s and h1_FidCutPiPlus==1 and h2_FidCutPiPlus==1 and (h1_p+h2_p)/nu<1.0'%(Nphe_h1_cut, Nphe_h2_cut,Chi2CC_h1_cut,Chi2CC_h2_cut,StatCC_h1_cut, StatCC_h2_cut)

pair_cut_nom_pi_p  ='h1_z>0.5 and abs(h2_deltaZ)<3.0 and abs(h1_deltaZ)<3.0 and TargType!=0 and SampFracEl25==1 and h1_FidCutPiPlus==1 and h2_FidCutPiPlus==1 and %s and %s and %s'%(Nphe_h1_cut,Chi2CC_h1_cut,StatCC_h1_cut)


def applyCuts(fullDataframe,name='default',isMC=False,isTrigger=True, nomCuts=False, minPt=None, verbose=False): 
    dataframe = fullDataframe
    if(dataframe.shape[0]>0) and verbose:
        print ('Entries before cut ', dataframe.shape[0])
    dataframe.eval('inelasticity = nu/5.014', inplace=True)
    dataframe.eval('h1_e = h1_z*nu', inplace=True)
    dataframe.eval('h1_p = sqrt(h1_e*h1_e-0.13957*0.13957)', inplace=True)
    if 'h1_Betta' in dataframe.columns:
        dataframe.eval('h1_mass_TOF = h1_p/h1_Betta*sqrt(1-h1_Betta**2)', inplace=True)
    dataframe = applyCut(dataframe, 'Q2>1.0 and Q2<4.0', 'Q2>1.0 and Q2<4.0' if verbose else None)
    dataframe = applyCut(dataframe, 'h1_p <5.0 and h1_p>0.2 ', '0.2<h1_p<5.0 ' if verbose else None)
    dataframe = applyCut(dataframe, 'inelasticity<0.85','inelasticity < 0.85' if verbose else None)
    dataframe = applyCut(dataframe, 'abs(h1_pid)==211', 'h1_pid = pions (trigger)'  if verbose else None)
    dataframe = applyCut(dataframe, 'nu>2.2 and nu<4.2', '2.2 < nu <4.2'  if verbose else None)
    
    #different names for the same variable: h1_cm_pt and h1_pt
    if 'h1_cm_pt' in dataframe.columns:
        #dataframe = applyCut(dataframe, 'h1_cm_pt*cos(h1_cm_ph)<0.4', 'px in gamma* frame is <0.4 GeV')
        if minPt is not None:
            dataframe = applyCut(dataframe, f'h1_cm_pt>{minPt}', f'pt>{minPt} GeV'  if verbose else None)
    elif 'h1_pt' in dataframe.columns:
        #dataframe = applyCut(dataframe, 'h1_pt*cos(h1_cm_ph)<0.4', 'px in gamma* frame is <0.4 GeV')
        if minPt is not None:
            dataframe = applyCut(dataframe, f'h1_pt>{minPt}', f'pt>{minPt} GeV'  if verbose else None)
    #dataframe = applyCut(dataframe, 'h1_cm_pt>sqrt(0.005)', 'pT^2>0.005 GeV^2 for leading hadron')
    
    #This cut is completely dumb
#     if("yes" in "is this a dumb cut? yes, it is" and not isMC):
#         h1_pz= "sqrt(h1_p**2-h1_cm_pt**2)*(2*(h1_rap>0)-1)"
#         theta_e="2*arcsin(sqrt(Q2/(4*5.014*(5.014-nu))))"
#         theta_q=f'arctan2((5.014-nu)*sin({theta_e}),(5.014-(5.014-nu)*cos({theta_e})))'
#         theta_qe=f"{theta_q}+{theta_e}"
#         e_pz = f"((5.014-nu)*cos(theta_qe))"
#         e_pt = f"((5.014-nu)*sin(theta_qe))"
#         df = dataframe
#         df['theta_e'] = df.eval(theta_e)
#         df['theta_qe'] = df.eval(theta_qe)
#         df['e_q_pz'] = df.eval(e_pz)
#         df['e_q_pt'] = df.eval(e_pt)
#         dtheta_epi1=f"""arccos((e_px*h1_px+e_py*h1_py+e_pz*h1_pz)
#                               /sqrt(e_px**2+e_py**2+e_pz**2)
#                               /sqrt(h1_px**2+h1_py**2+h1_pz**2))*{180/np.pi}""".replace("\n", "")

#         df['dtheta_epi1']  = df.eval(f"arccos((cos(theta_qe)*{h1_pz}+sin(theta_qe)*h1_cm_pt*cos(h1_cm_ph))/h1_p)*{180/np.pi}")
#         dataframe = applyCut(df, 'dtheta_epi1>0', '25 degree separation between pion and electron.  This is a dumb cut')
    

    
    
    if max(dataframe['h1_th'])<np.pi:
        dataframe['h1_th'] = dataframe['h1_th']*180/np.pi
    if (not isMC):
        #dataframe = applyCut(dataframe, 'h1_th<90 and h1_th>10', '10< h1_th<90')
        #dataframe = applyCut(dataframe, '(h1_pid>0) | (h1_pid==-211 & h1_th<90 & h1_th>25 & (h1_p>0.5 | h1_th>40))','Theta/P fiducial region selected for trigger')
        dataframe = applyCut(dataframe, 'h1_th>10', '10< h1_th'  if verbose else None)
        #dataframe = applyCut(dataframe, '(h1_pid>0) | (h1_pid==-211  & h1_th>25 & (h1_p>0.5 | h1_th>40))','Theta/P fiducial region selected for trigger')
        dataframe = applyCut(dataframe, '(h1_pid>0) | (h1_pid==-211 & h1_th>25 & h1_p>0.7) | (h1_pid==-211 & h1_p>0.5 & h1_th>30) | (h1_pid==-211 & h1_th>40)','Theta/P fiducial region selected for trigger hadron'  if verbose else None)
        #dataframe = applyCut(dataframe, '(h1_pid==211) | (h1_pid==-211 & h1_th>25 & h1_th<90) | (h1_pid==-211 & h1_th<40 & h1_th>25 & h1_p>0.5)','Theta/P fiducial region selected for trigger')
        if (nomCuts):    dataframe = applyCut(dataframe,trigger_cut_nom, 'Nom cuts for the trigger applied' if verbose else None)
        
    return dataframe

def applyCutsPair(fullDataframe,name='default',isMC=False,nomCuts=False,h2Proton=False, minPt=None, verbose=False):
    print ('Starting election on dipion variables')
    if (isMC):
        print ('This is MC')
    else: 
        print ('This is Data')
    dataframe = fullDataframe
    dataframe.eval('z_tot = h1_z+ h2_z', inplace=True)
    dataframe.eval('h1_e = h1_z*nu', inplace=True)
    dataframe.eval('h1_p = sqrt(h1_e*h1_e-0.13957*0.13957)', inplace=True)
    dataframe.eval('h2_e = h2_z*nu', inplace=True)
    
   
    
    
    
    
    m_h2 = '(0.13957*(abs(h2_pid)==211) + .93827208816*(abs(h2_pid)==2212))'
    dataframe.eval(f'h2_p = sqrt(h2_e*h2_e-{m_h2}**2)', inplace=True)
    if 'h1_Betta' in dataframe.columns:
        dataframe.eval('h1_mass2_TOF = h1_p**2/h1_Betta**2*(1-h1_Betta**2)', inplace=True)
        dataframe.eval('h2_mass2_TOF = h2_p**2/h2_Betta**2*(1-h2_Betta**2)', inplace=True)
    dataframe = applyCut(dataframe, 'Q2>1.0 and Q2<4.0', '1.0< Q2 <4.0')
    dataframe = applyCut(dataframe, 'nu>2.2 and nu<4.2', '2.2 < nu < 4.2')

    if 'pair_pt' not in dataframe.columns:
        dataframe.eval('pair_pt = sqrt(h1_cm_pt**2+h2_cm_pt**2+2*h2_cm_pt*h1_cm_pt*cos(h1_cm_ph-h2_cm_ph))', inplace=True)
    
    dataframe.eval('pair_pt2 = pair_pt*pair_pt', inplace=True)
    if h2Proton :
        dataframe = applyCut(dataframe, 'h2_pid==2212', 'secondary hadrons are protons' if verbose else None) 
    else :
        dataframe = applyCut(dataframe, 'abs(h2_pid)==211', 'secondary hadrons are pions' if verbose else None) 
        #dataframe = applyCut(dataframe, '(h1_p+h2_p)/nu<1.0', '(h1_p+h2_p)/nu<1.0') #already in NomCuts
        dataframe = applyCut(dataframe, '(h1_pid*h2_pid)<0.', '(h1_pid*h2_pid)<0, opposite charge pions' if verbose else None)
        dataframe = applyCut(dataframe, 'h2_z > %f and h2_z< %f'%(minz2,maxz2), 'kinematic limits for h2_z' if verbose else None)
        dataframe = applyCut(dataframe, 'dphi > %f and dphi< %f'%(mindphi,maxdphi), 'kinematic limits for dphi' if verbose else None)
        dataframe = applyCut(dataframe, 'pair_mass > %f and pair_mass< %f'%(minzmass,maxzmass), 'kinematic limits for pair_mass' if verbose else None)
    dataframe = applyCut(dataframe, 'abs(h1_pid)==211', 'leading hadrons are pions' if verbose else None)    
    dataframe = applyCut(dataframe, 'h2_p>0.2 and h2_p<5.0', '0.2<h2_p<5.0' if verbose else None)
    dataframe = applyCut(dataframe, 'h1_p>0.2 and h1_p<5.0', '0.2<h1_p<5.0' if verbose else None)
    
    if 'h1_cm_pt' in dataframe.columns:
        #dataframe = applyCut(dataframe, 'h1_cm_pt*cos(h1_cm_ph)<0.4', 'px in gamma* frame is <0.4 GeV')
        if minPt is not None:
            dataframe = applyCut(dataframe, f'h1_cm_pt>{minPt} and h2_cm_pt>{minPt}', f'pt>{minPt} GeV' if verbose else None)
    elif 'h1_pt' in dataframe.columns:
        #dataframe = applyCut(dataframe, 'h1_pt*cos(h1_cm_ph)<0.4', 'px in gamma* frame is <0.4 GeV')
        if minPt is not None:
            dataframe = applyCut(dataframe, f'h1_pt>{minPt} and h2_pt>{minPt}', f'pt>{minPt} GeV' if verbose else None)
    
    #This cut is completely dumb
#     if("yes" in "is this a dumb cut? yes, it is" and not isMC):
#         h1_pz= "sqrt(h1_p**2-h1_cm_pt**2)*(2*(h1_rap>0)-1)"
#         theta_e="2*arcsin(sqrt(Q2/(4*5.014*(5.014-nu))))"
#         theta_q=f'arctan2((5.014-nu)*sin({theta_e}),(5.014-(5.014-nu)*cos({theta_e})))'
#         theta_qe=f"{theta_q}+{theta_e}"
#         e_pz = f"((5.014-nu)*cos(theta_qe))"
#         e_pt = f"((5.014-nu)*sin(theta_qe))"
#         df = dataframe
#         df['theta_e'] = df.eval(theta_e)
#         df['theta_qe'] = df.eval(theta_qe)
#         df['e_q_pz'] = df.eval(e_pz)
#         df['e_q_pt'] = df.eval(e_pt)
#         dtheta_epi1=f"""arccos((e_px*h1_px+e_py*h1_py+e_pz*h1_pz)
#                               /sqrt(e_px**2+e_py**2+e_pz**2)
#                               /sqrt(h1_px**2+h1_py**2+h1_pz**2))*{180/np.pi}""".replace("\n", "")

#         df['dtheta_epi1']  = df.eval(f"arccos((cos(theta_qe)*{h1_pz}+sin(theta_qe)*h1_cm_pt*cos(h1_cm_ph))/h1_p)*{180/np.pi}")
#         dataframe = applyCut(df, 'dtheta_epi1>0', '25 degree separation between pion and electron.  This is a dumb cut')

    
    
    #dataframe = applyCut(dataframe, 'h2_cm_pt>sqrt(0.005) and h1_cm_pt>sqrt(0.005)', 'pT^2>0.005 GeV^2 for both hadrons')
    ## Theta cuts are not applied in the GiBUU case
    ## For GiBUU case is the only time isMC=True
    if (not isMC):
        #convert to degrees if necessary
        if max(dataframe['h2_th'])<np.pi:
            dataframe['h2_th'] = dataframe['h2_th']*180/np.pi
        if max(dataframe['h1_th'])<np.pi:
            dataframe['h1_th'] = dataframe['h1_th']*180/np.pi
        #dataframe = applyCut(dataframe, 'h1_cm_pt*cos(h1_cm_ph)<0.4', 'px in gamma* frame is <0.4 GeV')


        #dataframe = applyCut(dataframe, 'h2_th<90 and h2_th>10', '10<h2_th<90')
        #dataframe = applyCut(dataframe, '(h2_pid>0) | (h2_pid==-211 & h2_th<90 & h2_th>25 & (h2_p>0.5 | h2_th>40))','Theta/P fiducial region selected for secondary hadron')
        #dataframe = applyCut(dataframe, 'h1_th<90 and h1_th>10', '10< h1_th<90')
        #dataframe = applyCut(dataframe, '(h1_pid>0) | (h1_pid==-211 & h1_th<90 & h1_th>25 & (h1_p>0.5 | h1_th>40))','Theta/P fiducial region selected for trigger hadron')
        
        dataframe = applyCut(dataframe, 'h2_th>10', '10<h2_th')
        dataframe = applyCut(dataframe, '(h2_pid>0) | (h2_pid==-211 & h2_th>25 & h2_p>0.7) |(h2_pid==-211 & h2_p>0.5 & h1_th>30) | (h2_pid==-211 & h2_th>40)','Theta/P fiducial region selected for secondary hadron' if verbose else None)
        dataframe = applyCut(dataframe, ' h1_th>10', '10< h1_th')
        dataframe = applyCut(dataframe, '(h1_pid>0) | (h1_pid==-211 & h1_th>25 & h1_p>0.7) |(h1_pid==-211 & h1_p>0.5 & h1_th>30) | (h1_pid==-211 & h1_th>40)','Theta/P fiducial region selected for trigger hadron' if verbose else None)
        
        
        if (nomCuts):
            if h2Proton:
                dataframe = applyCut(dataframe, pair_cut_nom_pi_p, 'Nom cuts for the pair applied (pi p)' if verbose else None)
            else :
                dataframe = applyCut(dataframe, pair_cut_nom, 'Nom cuts for the pair applied (pi pi)' if verbose else None)
    return dataframe

def printPairBreakdown(dataframe,h2Proton=True):
    allpairs = 1.0*dataframe.shape[0]
    print ('All pairs ', allpairs)
    print ('Pairs with Leading pi+', np.true_divide(dataframe.query('h1_pid==211').shape[0],allpairs))
    print ('Pairs with Leading pi-', dataframe.query('h1_pid==-211').shape[0]/allpairs)
    print ('Pairs with Sub-Leading pi+', dataframe.query('h2_pid==211').shape[0]/allpairs)
    print ('Pairs with Sub-Leading pi-', dataframe.query('h2_pid==-211').shape[0]/allpairs)
    print ('Pairs with Sub-Leading p', dataframe.query('h2_pid==2212').shape[0]/allpairs)
    print ('pi+ pi+ pairs',dataframe.query('h1_pid==211 and h2_pid==211').shape[0]/allpairs)
    print ('pi- pi- pairs',dataframe.query('h1_pid==-211 and h2_pid==-211').shape[0]/allpairs)
    print ('pi+ pi- pairs',dataframe.query('h1_pid==211 and h2_pid==-211').shape[0]/allpairs)
    print ('pi- pi+ pairs',dataframe.query('h1_pid==-211 and h2_pid==211').shape[0]/allpairs)
    print ('pi+ p pairs',dataframe.query('h1_pid==211 and h2_pid==2212').shape[0]/allpairs)
    print ('pi- p pairs',dataframe.query('h1_pid==-211 and h2_pid==2212').shape[0]/allpairs)
    
    print ('//////////////////////////////////////////////////////')
    return



