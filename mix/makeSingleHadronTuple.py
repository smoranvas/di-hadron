import ROOT
from ROOT import TFile
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import pandas as pd 
import time
from os import listdir

import root_pandas

class particle:
    def __init__(self, pid, fourvector, virtual_photon,ThetaPQ ):## '__init__' is the constructor of the class
        self.virtual_photon = virtual_photon
        Nu = virtual_photon.E()   ##components of a 4-vector TLorentzVector
        Q2 = -virtual_photon.M2() ## magnitud squared of a 4-vector TLorentzVector
        self.proton = ROOT.TLorentzVector()  ## proton is an attribute of the class 'particle' just created.
        self.proton.SetPxPyPzE(0,0,0, 0.938)  ## SetPxPyPzE is just a function of ROOT
        self.W = (virtual_photon + self.proton).M() ##.M() return the magnitud of a TLorentzVector [W2=(p+q)2]
        incoming_e = ROOT.TLorentzVector()
        incoming_e.SetPxPyPzE(0,0,5.014,5.014)
        part1 = virtual_photon.Vect().Cross(incoming_e.Vect()).Unit()
        part2 = virtual_photon.Vect().Cross(fourvector.Vect()).Unit()
        sign  = np.sign(part1.Dot(fourvector.Vect())) ## sign returns -1 0 or 1 if the input is negative, zero or positive.
        self.PhiPQ = sign*np.arccos(part1.Dot(part2))
        photon_pz = np.sqrt(Nu*Nu+Q2) #direction is positive by definition
        self.bcm = photon_pz/(Nu + 0.938)#photon-nucleon center-of-mass velocity 
        self.ycm = 0.5*np.log(( 1+self.bcm)/(1-self.bcm)) #photon-nucleon center-of-mass rapidity
        self.LorentzVector = fourvector #hadron four-vector. 4-vector is an input of this class
        self.PhiLab = self.LorentzVector.Phi()
        self.ThetaLab = self.LorentzVector.Theta()
        self.E = self.LorentzVector.E() #energy in lab frame
        self.vector = self.LorentzVector.Vect()
        self.Pt = self.vector.Perp(virtual_photon.Vect().Unit()) #pT with respect to photon direction
        self.Pl  = self.vector.Dot(virtual_photon.Vect().Unit()) #pL with respect to photon direction (in lab frame)
        self.y =  0.5*np.log( (self.E+self.Pl)/(self.E-self.Pl)) #rapidity in lab frame
        self.mT = np.sqrt(self.LorentzVector.M2() + self.Pt*self.Pt)
        self.y_star = self.y - self.ycm
        self.Pl_star = self.mT*np.sinh(self.y_star) # y is rapidity
        self.Xf = 2.0*self.Pl_star/self.W 
        self.pid = pid
        self.Zh = self.E/Nu
        self.ThetaPQ = np.arctan(self.Pt/self.Pl)
        self.P = np.sqrt(self.LorentzVector.Px()**2+self.LorentzVector.Py()**2+self.LorentzVector.Pz()**2)
        #self.Nphe=Nphe
        #self.deltaZ=deltaZ
        #self.fidCut=fid
        
    def redefine(self, new_virtual_photon):
        incoming_e = ROOT.TLorentzVector()
        incoming_e.SetPxPyPzE(0,0,5.014,5.014)
        part1 = new_virtual_photon.Vect().Cross(incoming_e.Vect()).Unit()
        part2 = new_virtual_photon.Vect().Cross(self.LorentzVector.Vect()).Unit()
        sign  = np.sign(part1.Dot(self.LorentzVector.Vect()))
        self.PhiPQ = sign*np.arccos(part1.Dot(part2)) 
        self.Pt = self.LorentzVector.Vect().Perp(new_virtual_photon.Vect().Unit()) #pT with respect to photon direction
        self.Pl  = self.LorentzVector.Vect().Dot(new_virtual_photon.Vect().Unit()) #pL with respect to photon direction (in lab frame)
        self.y =  0.5*np.log( (self.E+self.Pl)/(self.E-self.Pl)) #rapidity in lab frame
        self.ThetaPQ = np.arctan(self.Pt/self.Pl)
        self.virtual_photon = new_virtual_photon
        
        return
        
    def print_properties(self):
        print ('Hello, let me introduce myself, i am particle pid = ' , self.pid, ' with index ', self.index, ', from event  #', self.ievt, ' Nu and W', self.Nu, ' ' , self.W)
        print ('zh = ', self.Zh, ' phi_pq= ', self.PhiPQ, ' theta_pq=' , self.ThetaPQ, 'E = ', self.E, ' xf', self.Xf,'Pt ', self.Pt, ' Pl= ', self.Pl, ' rapidity=' ,  self.y)
        print ('pid = ' , self.pid)

Ebeam=5.014

def getDataframes(filename, Target=1,maxevents=1e9,tree_name='ntuple_data',isMC=False):
    dphi = np.array([])  ## standard way to create an array using numpy ('nampai') numerical python
    
    ParticlesFromPrevious = [] ## square brackets is for lists in python, we can use '.append' for this for example.
    
    try:
        myfile = TFile.Open('%s'%filename,'READ')
        myfile.Print()
    except:
        print("could not open file")
    mytree = myfile.Get(tree_name)
        
    print (filename, ' has ', mytree.GetEntries(), ' entries')
    
    tupla = {}  ## this is how you define a dictionary, with curly brcaes {}
    tupla['dphi'] = []  ## this is how you create a new entry to the dictionary and how you set it to be a list ([])
    tupla['dphi_lab'] = []
    tupla['drap'] = []
    tupla['h1_z'] = [] 
    tupla['h2_z'] = []
    tupla['h1_cm_pt'] = []
    tupla['h2_cm_pt'] = []
    tupla['h1_xf'] = []
    tupla['h2_xf'] = []
    tupla['h1_rap'] = []
    tupla['ycm'] = []
    tupla['h2_rap'] = []
    tupla['h1_pid'] = []
    tupla['h2_pid'] = []
    tupla['h1_cm_ph'] = []
    tupla['h2_cm_ph'] = []
    tupla['h1_cm_th'] = []  
    tupla['h2_cm_th'] = []  
    tupla['pair_mass'] = []
    tupla['pair_pt'] = []
    tupla['mx_eh1h2x'] = []
    tupla['mx_eh1x'] = []
    tupla['mx_eh2x'] = []
    tupla['t']  = []
    tupla['Q2'] = [] 
    tupla['nu'] = []
    tupla['W']  = []
    tupla['x'] = []
    tupla['u']  = []
    tupla['h1_ph'] = []
    tupla['h1_th'] = []
    tupla['h2_ph'] = []
    tupla['h2_th'] = []
    
    ## here we create another dictionary
    tupla_mix = {}
    tupla_mix['dphi'] = []
    tupla_mix['dphi_lab'] = []
    tupla_mix['drap'] = []
    tupla_mix['h1_z'] = []
    tupla_mix['h2_z'] = []
    tupla_mix['h1_cm_pt'] = []
    tupla_mix['h2_cm_pt'] = []
    tupla_mix['h1_xf'] = []
    tupla_mix['h2_xf'] = []
    tupla_mix['h1_rap'] = []
    tupla_mix['ycm'] = []
    tupla_mix['h2_rap'] = []
    tupla_mix['h1_pid'] = []
    tupla_mix['h2_pid'] = []
    tupla_mix['h1_cm_ph']   = []
    tupla_mix['h2_cm_ph'] = []
    tupla_mix['h1_cm_th'] = []
    tupla_mix['h2_cm_th'] = []
    tupla_mix['pair_mass'] = []
    tupla_mix['pair_pt'] = []
    tupla_mix['mx_eh1h2x'] = []
    tupla_mix['mx_eh1x'] = []
    tupla_mix['mx_eh2x'] = []
    tupla_mix['t']  = []
    tupla_mix['Q2'] = []
    tupla_mix['nu'] = []
    tupla_mix['W']  = []
    tupla_mix['x'] = []
    tupla_mix['u']  = []
    tupla_mix['h1_ph'] = []
    tupla_mix['h1_th'] = []
    tupla_mix['h2_ph'] = []
    tupla_mix['h2_th'] = []
    tupla_mix['dphi_norot'] = [] #save variable before rotation to new virtual photon frame
    tupla_mix['h2_cm_ph_norot'] = [] #save variable before rotation to new virtual photon frame
    tupla_mix['h2_cm_th_norot'] = [] #save variable before rotation to new virtual photon frame
    
    ## here we create another dictionary
    tupla_trigger = {}
    tupla_trigger['E'] = []
    tupla_trigger['e_p'] = []
    tupla_trigger['e_th'] = []
    tupla_trigger['e_ph'] = []
    tupla_trigger['h_pid'] = []
    tupla_trigger['h_xf'] = []
    tupla_trigger['h_xf_default'] = []
    tupla_trigger['h_z']  = []
    tupla_trigger['h_cm_pt'] = []
    tupla_trigger['h_rap']  = []
    tupla_trigger['h_cm_rap'] = []
    tupla_trigger['Q2'] = []
    tupla_trigger['x'] = []
    tupla_trigger['nu'] = []
    tupla_trigger['W'] = []
    tupla_trigger['h_cm_ph'] = []
    tupla_trigger['h_cm_th'] = []
    tupla_trigger['TargType'] = []
    tupla_trigger['missing_mass'] = []
    tupla_trigger['h_ph'] = []
    tupla_trigger['h_th'] = []
    tupla_trigger['h_p'] = []
    tupla_trigger['h_deltaZ'] = []
    tupla_trigger['h_Nphe'] = []
    tupla_trigger['h_Sector'] = []
    tupla_trigger['h_FidCut'] = []
    tupla_trigger['h_Chi2CC'] = []
    tupla_trigger['h_StatCC'] = []
    tupla_trigger['SampFracEl25'] = []
    tupla_trigger['SampFracEl20'] = []
    
    start = time.time()
    print('About to loop over ', mytree.GetEntries() , ' entries')
    for ievt  in range(mytree.GetEntries()):
        #print('evnt: ',ievt, 'len ParticlesFromPrevious ', (len(ParticlesFromPrevious)))
        
        #for kkk  in range (len(ParticlesFromPrevious)):
        #    print('evnt',ievt, 'ParticlesFromPrevious, pid:', ParticlesFromPrevious[kkk].pid, 'zh: ', ParticlesFromPrevious[kkk].Zh, 'W: ',ParticlesFromPrevious[kkk].W )

        mytree.GetEntry(ievt)   
        if mytree.W<2.0 or mytree.Q2<1.0: continue
        
        if ievt>maxevents: break        
        if(mytree.TargType==1):
            TargType=1
        elif(mytree.TargType==2):
            TargType=2
        else:
            TargType=0
        #print (TargType,  ' ' , Target)
        if not(isMC) and (TargType!=Target): continue ## 'Target' is a argument of this function
                
        #print('PASO all cuts')
        W = mytree.W
        Nu = mytree.Nu
        #get electron momentum:
        Pe = np.sqrt(mytree.Pex*mytree.Pex + mytree.Pey*mytree.Pey+ mytree.Pez*mytree.Pez)
        scattered_e = ROOT.TLorentzVector()
        scattered_e.SetPxPyPzE(mytree.Pex, mytree.Pey, mytree.Pez, Pe)
        incoming_e = ROOT.TLorentzVector()
        incoming_e.SetPxPyPzE(0,0,5.014,5.014)
        virtual_photon  = incoming_e - scattered_e 
        virtual_photon_unitvector = virtual_photon.Vect().Unit()
        proton = ROOT.TLorentzVector()
        proton.SetPxPyPzE(0,0,0, 0.938)
        ## for each event the particles are saved in the 'particles' list.        
        particles = []  ## this is how you define a list in python, this is created for each event
        #print (' Entering main loop over particles')
        for i in range(len(mytree.pid)):
            #print(mytree.pid[i])
            ## when the condition is true the 'continue' statement
            ## takes me to the next iteration of the loop
            if (abs(mytree.pid[i]) !=211 and mytree.pid[i]!=2212): continue ## with this I make sure that if I 
                                                                            ## have a pid=0 for example the 
                                                                            ## rest of the code is dismiss and 
                                                                            ## we skip to the next iteration
                            
            ## aplying the rest of cuts:
            if ( mytree.FidCheckCutPiPlus[i]==False or abs(mytree.deltaZ[i])>=3.0 or mytree.P[i]<0.2 or mytree.ThetaLab[i]<10 or mytree.ThetaLab[i]>120 ):    continue
            if ( mytree.pid[i]==-211 and (mytree.ThetaLab[i]<25) ):    continue                
            if ( mytree.pid[i]==-211 and mytree.ThetaLab[i]>40 and mytree.P[i]<0.2  ):    continue
            if ( mytree.pid[i]==-211 and mytree.ThetaLab[i]<40 and mytree.P[i]<0.5  ):    continue 
            if ( mytree.pid[i]==-211 and mytree.ThetaLab[i]<30 and mytree.P[i]<0.7  ):    continue
            if ( mytree.pid[i]==211  and mytree.P[i]>=2.7 and (mytree.Nphe[i]<15 or mytree.Chi2CC[i]>0.0872 or mytree.StatCC[i]<=0 or mytree.NRowsCC[i]==0) ):    continue     

            #print(mytree.pid[i])
            i_lv = ROOT.TLorentzVector()    ## 4-vector of the hadron
            i_lv.SetPxPyPzE(mytree.Px[i],mytree.Py[i],mytree.Pz[i],mytree.Zh[i]*Nu) ## this is the 4-vector of the hadron
            i_part = particle(mytree.pid[i], i_lv, virtual_photon, mytree.ThetaPQ[i] ) ## particle is the class defined previously
            ## in this 'particles' list there are NO cuts applied (except the pid and obvious ones)
            particles.append(i_part)   ## save that particle in the 'particles' list
            X = (virtual_photon + proton -  i_part.LorentzVector) #unobserved hadronic system
            #print('event:',ievt,'particle i:', i,' with PID:',  i_part.pid, 'and Zh:', i_part.Zh, ', W: ',i_part.W)
            if i_part.Zh > 0: #only save triggers and do correlations if they have z>0.4
            #if i_part.Pt>1.0: #only save triggers if pT>1.0 
            ## HERE WE SAVE THE VARIABLES FOR THE TRIGGER PARTICLE (THE ONE WITH Zh>0.4)
                tupla_trigger['TargType'].append(mytree.TargType)
                tupla_trigger['h_pid'].append(i_part.pid)
                tupla_trigger['h_xf'].append(i_part.Xf)
                tupla_trigger['h_xf_default'].append(-1)
                tupla_trigger['h_z'].append(i_part.Zh)
                tupla_trigger['h_cm_pt'].append(i_part.Pt)
                tupla_trigger['h_rap'].append(i_part.y_star)
                tupla_trigger['h_cm_rap'].append(i_part.ycm)
                tupla_trigger['h_cm_ph'].append(i_part.PhiPQ)
                tupla_trigger['h_cm_th'].append(i_part.ThetaPQ)
                tupla_trigger['missing_mass'].append(X.M())
                tupla_trigger['e_p'].append(mytree.Pe)
                tupla_trigger['e_th'].append(np.arctan2(np.hypot(mytree.Pex,mytree.Pey),mytree.Pez))
                tupla_trigger['e_ph'].append(np.arctan2(mytree.Pey,mytree.Pex))
                tupla_trigger['E'].append(Ebeam)
                tupla_trigger['Q2'].append(mytree.Q2)
                tupla_trigger['x'].append(mytree.Xb)
                tupla_trigger['nu'].append(mytree.Nu)
                tupla_trigger['W'].append(mytree.W)
                tupla_trigger['h_p'].append(mytree.P[i])
                tupla_trigger['h_ph'].append(mytree.PhiLab[i]*np.pi/180)
                tupla_trigger['h_th'].append(mytree.ThetaLab[i]*np.pi/180)
                tupla_trigger['h_deltaZ'].append(mytree.deltaZ[i])
                tupla_trigger['h_Nphe'].append(mytree.Nphe[i])
                tupla_trigger['h_Sector'].append(mytree.Sector[i])
                tupla_trigger['h_FidCut'].append(mytree.FidCheckCutPiPlus[i])
                tupla_trigger['h_Chi2CC'].append(mytree.Chi2CC[i])
                tupla_trigger['h_StatCC'].append(mytree.StatCC[i])
                tupla_trigger['SampFracEl25'].append(mytree.SampFractionEl25)
                tupla_trigger['SampFracEl20'].append(mytree.SampFractionEl20)
                #i_part.LorentzVector.Theta())
                #print 'mytree.Pt[i] ' , mytree.Pt[i], ' check: ' ,i_part.Vector.Perp(virtual_photon_unitvector)
                
                #print('Testing theta PQ', mytree.ThetaPQ[i],  ' '  , 180.0*i_part.ThetaPQ/np.pi)
                #print('Testing phi PQ', mytree.PhiPQ[i],  ' '  , 180.0*i_part.PhiPQ/np.pi)
#,Nphe,deltaZ                
                '''for j in range(len(mytree.pid)): 
                    if i==j: continue
                    if (abs(mytree.pid[j]) !=211 and mytree.pid[j]!=2212): continue
                    #print('evnt', ievt,' j:', j, ' , lenpid:',len(mytree.pid) )
                    ## aplying the rest of cuts:
                    if ( mytree.FidCheckCutPiPlus[j]==False or abs(mytree.deltaZ[j])>=3.0 or mytree.P[j]<0.2 or mytree.ThetaLab[j]<10 or mytree.ThetaLab[j]>120 ):    continue
                    if ( mytree.pid[j]==-211 and (mytree.ThetaLab[j]>90 or mytree.ThetaLab[j]<25) ):    continue                
                    if ( mytree.pid[j]==-211 and mytree.ThetaLab[j]>40 and mytree.P[j]<0.2  ):    continue                
                    if ( mytree.pid[j]==-211 and mytree.ThetaLab[j]<40 and mytree.P[j]<0.5  ):    continue 
                    if ( mytree.pid[j]==211  and mytree.P[j]>=2.7 and (mytree.Nphe[j]<15 or mytree.Chi2CC[j]>0.0872 or mytree.StatCC[j]<=0 or mytree.NRowsCC[j]==0) ):    continue     
                        
                    #print('inside j loop')
                    j_lv = ROOT.TLorentzVector()    
                    j_lv.SetPxPyPzE(mytree.Px[j],mytree.Py[j],mytree.Pz[j],mytree.Zh[j]*Nu)
                    j_part = particle(mytree.pid[j], j_lv, virtual_photon, mytree.ThetaPQ[j] ) ## particle is the defined class

                    
                    ## TVector2.Phi_mpi_pi is a built in function of TVector2 class (2D vectors).
                    ## Returns Phi angle in the interval [-pi,pi]. 
                    ## I used just to have the angles in the rangle between -pi and pi 
                    ## Also, we used the 'abs' to obtain the 0->pi range
                    
                    dphi = abs(ROOT.TVector2.Phi_mpi_pi(i_part.PhiPQ-j_part.PhiPQ))  
                    dphi_lab = abs(ROOT.TVector2.Phi_mpi_pi(i_part.PhiLab-j_part.PhiLab))
                    
                    dy = i_part.y-j_part.y  ## i_part and j_part are object of the class particles 
                                            # and 'y' is a methos, is the rapidity in the lab frame
                    deta = dy
                    ## Di-hadron 4-vector, is the sum P_{di-hadron} = P_{trigger-hadron} + P_{other-hadron}
                    dipion = i_part.LorentzVector+j_part.LorentzVector  ## LorentzVector is a method of the 
                                    #particle class, is the '4-vector' of the particle, 3rd argument of the class
                    
                    X = (virtual_photon + proton - dipion) #unobserved hadronic system
                    X1 = (virtual_photon + proton - i_part.LorentzVector)
                    X2 = (virtual_photon + proton - j_part.LorentzVector)

                    tupla['dphi'].append(dphi)
                    tupla['dphi_lab'].append(dphi_lab)
                    tupla['drap'].append(dy)
                    tupla['h1_z'].append(i_part.Zh)
                    tupla['h2_z'].append(j_part.Zh)
                    tupla['h1_cm_pt'].append(i_part.Pt)
                    tupla['h2_cm_pt'].append(j_part.Pt)
                    tupla['h1_xf'].append(i_part.Xf)
                    tupla['h2_xf'].append(j_part.Xf)
                    tupla['h1_rap'].append(i_part.y_star)
                    tupla['ycm'].append(i_part.ycm)
                    tupla['h2_rap'].append(j_part.y_star)
                    tupla['h1_pid'].append(i_part.pid)
                    tupla['h2_pid'].append(j_part.pid)
                    tupla['h1_cm_ph'].append(i_part.PhiPQ)
                    tupla['h2_cm_ph'].append(j_part.PhiPQ)
                    tupla['h1_cm_th'].append(i_part.ThetaPQ)
                    tupla['h2_cm_th'].append(j_part.ThetaPQ)
                    tupla['pair_mass'].append(dipion.M())
                    tupla['pair_pt'].append( dipion.Vect().Perp(virtual_photon_unitvector))
                    tupla['mx_eh1h2x'].append(X.M())
                    tupla['mx_eh1x'].append(X1.M())
                    tupla['mx_eh2x'].append(X2.M())
                    tupla['t'].append( -(virtual_photon- dipion).M2())
                    tupla['Q2'].append(mytree.Q2)
                    tupla['x'].append(mytree.Xb)
                    tupla['nu'].append(mytree.Nu)
                    tupla['W'].append(mytree.W)
                    tupla['u'].append(-(scattered_e-proton).M2())
                    tupla['h1_ph'].append(mytree.PhiLab[i])
                    tupla['h1_th'].append(mytree.ThetaLab[i])
                    tupla['h2_ph'].append(mytree.PhiLab[j])
                    tupla['h2_th'].append(mytree.ThetaLab[j])
         #end loop over secondary loop    
                #print '//////// Entering mixed event correlations with # ', len(ParticlesFromPrevious) , ' paticles in previous event'
                
                ## here we are still under the condition of Zh>0.4
                #print(ievt,i,j)
                #print('')
                for mixparticle in ParticlesFromPrevious: ## ParticlesFromPrevious is a list 
                                                          ## with 'particle' class objects
                    #print('\ninside mixparticle loop\n')
                    #print('i: ',i,'mixparticles pid:',mixparticle.pid,' zh :', mixparticle.Zh , ', W: ', mixparticle.W, ' i_part Zh:', i_part.Zh, 'i_part.W: ',i_part.W)
                    dphi = abs(ROOT.TVector2.Phi_mpi_pi(i_part.PhiPQ-mixparticle.PhiPQ))
                    tupla_mix['dphi_norot'].append(dphi)
                    tupla_mix['h2_cm_ph_norot'].append(mixparticle.PhiPQ)
                    tupla_mix['h2_cm_th_norot'].append(mixparticle.ThetaPQ)

                    mixparticle.redefine(virtual_photon) 
                    #recalculates variables in this' event photon frame (not in the previous one)
                    #dphi = abs(ROOT.TVector2.Phi_mpi_pi(i_part.PhiPQ-mixparticle.PhiPQ))
                    #print 'dphi_pq after redefinition: ', dphi , ' phi_pq ', mixparticle.PhiPQ
                    
                    dipion = i_part.LorentzVector+mixparticle.LorentzVector
                    X  = (virtual_photon + proton - dipion)
                    X1 = (virtual_photon + proton - i_part.LorentzVector)
                    X2 = (virtual_photon + proton - mixparticle.LorentzVector)

                    #recalculate the phi_pq. It has to be with respect to the photon direction
                    
                    dphi = abs(ROOT.TVector2.Phi_mpi_pi(i_part.PhiPQ-mixparticle.PhiPQ))
                    dphi_lab = abs(ROOT.TVector2.Phi_mpi_pi(i_part.PhiLab-mixparticle.PhiLab))

                    dy = i_part.y-mixparticle.y
                    deta = dy#i_part.ThetaPQ-mixparticle.ThetaPQ
                    tupla_mix['dphi'].append(dphi)
                    tupla_mix['dphi_lab'].append(dphi_lab)
                    tupla_mix['drap'].append(dy)
                    tupla_mix['h1_z'].append(i_part.Zh)
                    tupla_mix['h2_z'].append(mixparticle.Zh)
                    tupla_mix['h1_cm_pt'].append(i_part.Pt)
                    tupla_mix['h2_cm_pt'].append(mixparticle.Pt)
                    tupla_mix['h1_xf'].append(i_part.Xf)
                    tupla_mix['h2_xf'].append(mixparticle.Xf)
                    tupla_mix['h1_rap'].append(i_part.y_star)
                    tupla_mix['ycm'].append(i_part.ycm)
                    tupla_mix['h2_rap'].append(mixparticle.y_star)
                    tupla_mix['h1_pid'].append(i_part.pid)
                    tupla_mix['h2_pid'].append(mixparticle.pid)
                    tupla_mix['h1_cm_ph'].append(i_part.PhiPQ)
                    tupla_mix['h2_cm_ph'].append(mixparticle.PhiPQ)
                    tupla_mix['h1_cm_th'].append(i_part.ThetaPQ)
                    tupla_mix['h2_cm_th'].append(mixparticle.ThetaPQ)
                    tupla_mix['pair_mass'].append(dipion.M())
                    tupla_mix['pair_pt'].append( dipion.Vect().Perp(virtual_photon_unitvector))
                    tupla_mix['mx_eh1h2x'].append(X.M())
                    tupla_mix['mx_eh1x'].append(X1.M())
                    tupla_mix['mx_eh2x'].append(X2.M())
                    tupla_mix['t'].append( -(virtual_photon- dipion).M2())
                    tupla_mix['Q2'].append(mytree.Q2)
                    tupla_mix['x'].append(mytree.Xb)
                    tupla_mix['nu'].append(mytree.Nu)
                    tupla_mix['W'].append(mytree.W)
                    tupla_mix['u'].append(-(scattered_e-proton).M2())
                    tupla_mix['h1_ph'].append(i_part.LorentzVector.Phi())
                    tupla_mix['h1_th'].append(mytree.ThetaLab[i])
                    tupla_mix['h2_ph'].append(mixparticle.LorentzVector.Phi())
                    tupla_mix['h2_th'].append(mytree.ThetaLab[j])
                    #for kk  in range (len(ParticlesFromPrevious)):
                        #print('ParticlesFromPrevious, pid:', ParticlesFromPrevious[kk].pid, 'zh: ', ParticlesFromPrevious[kk].Zh, 'W: ',ParticlesFromPrevious[kk].W )
        #print (' Exiting main loop over particles (i loop, not over all entries)')
        ParticlesFromPrevious = particles
        #print ' going for next event'    
        #print ' particles in event', len(particles
        ##end loop over events correlations    
               '''
    end = time.time()
    print ('Processed in',  end-start, 'seconds')
    ##printing the 3 tuples to the output file
    #df = pd.DataFrame(tupla)
    #df_mix= pd.DataFrame(tupla_mix)
    df_trigger = pd.DataFrame(tupla_trigger)
    #print ('Number of triggers with z>0.4,  ', df.query('h1_z>0.4').shape[0])
    #print ('Number of pairs with z>0.4, '    , df_trigger.query('h1_z>0.4').shape[0]) 
    myfile.Close()
    return df_trigger


if __name__ == '__main__':
    import sys
    input_filename,input_treename = sys.argv[1].split(":")
    output_filename,output_treename = sys.argv[2].split(":")
    isMC = "--MC" in sys.argv[3:]
    Target = 1
    maxevents = 1e9
    for i in range(len(sys.argv)):
        if '-N' == sys.argv[i]:
            i+=1
            maxevents = int(sys.argv[i])
        if '-T' == sys.argv[i]:
            i+=1
            Target = int(sys.argv[i])
    df_trigger =getDataframes(input_filename, Target=Target,maxevents=maxevents,tree_name=input_treename,isMC=isMC)
    df_trigger.to_root(output_filename,output_treename)
