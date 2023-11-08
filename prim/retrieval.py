#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 18:02:30 2023

@author: idajandl@gmx.de
"""
import os
from warnings import warn

import numpy as np
import matplotlib.pyplot as plt
from numpy import ma
import pickle as pkl


from prim.helper import slit_conv, StopExecution
from prim.optics import Optics
from prim.sun import Sun
from prim.atmosphere import Atmosphere
from prim.surface import Surface
from prim.prisma import PRISMA

test= []
class PRISCO:
    def __init__(self,settings,location):
        self.location=location
        self.settings=settings
        
        self.atm = Atmosphere(settings,location)
        self.surface = Surface(settings)
        
        self.optics = Optics(settings)
        print(self.optics.species)
        self.measurement = PRISMA(settings,location)
        
        meas=self.measurement.read_netcdf4(0,0)
        self.sun_lbl= Sun(self.settings, meas['wave']).lbl
        #self.pklfile = '../tmp/optics_prop_ch4_2.pkl'
        self.pklfile = '../../tmp/optics_prop'+self.location.name+'.pkl' #TODO:
        if os.path.exists(self.pklfile):
            self.optics.prop = pkl.load(open(self.pklfile,'rb'))

        else:
            if not os.path.exists("tmp"):
                os.makedirs("tmp")
            self.optics.cal_molec_xsec(self.atm)
            pkl.dump(self.optics.prop,open(self.pklfile,'wb'))
        self.ref_CO2,self.ref_CH4,self.ref_H2O,self.ref_N2O,self.ref_NH3 = [self.atm.atmo['CO2'],self.atm.atmo['CH4'],self.atm.atmo['H2O'],self.atm.atmo['N2O'],self.atm.atmo['NH3']]
        
        

    
    def transmission(self, surface, mu0, muv):
        """
        Calculate transmission solution given
        geometry (mu0,muv) using matrix algebra
        
        Parameters 
        ----------
        optics: optic_prop object
        surface: surface_prop object
        mu0: cosine of the solar zenith angle [-]
        muv: cosine of the viewing zenith angle [-]
        
        Returns
        ------- 
        rad_trans: single scattering relative radiance [wavelength] [1/sr]
        """    
        while True:
            if 0. <= mu0 <= 1. and -1. <= muv <= 1.:
               break
            else:
               print("ERROR! transmission: input out of range.")
               raise StopExecution

        
        tautot    = self.optics.prop['tautot'][:]
        mueff     = abs(1./mu0) + abs(1./muv)
        fact      = mu0/np.pi
        exptot    = np.exp(-tautot*mueff)
        
        rad_trans = self.sun_lbl*fact*surface.alb*exptot 
        dev_tau   =-mueff*rad_trans
        dev_alb   = fact*exptot*self.sun_lbl
        
        return rad_trans, dev_tau, dev_alb    
    
    def nonscat_fwd_model(self,slit, atm, surface, mu0, muv, fwhm):
        """
        Function that convolutes the synthetic forward model with the IRSF (Gauss)
        
        """
        #species = [spec for spec in dev if spec[0:5]=='molec']
        
        # Init class with current_optics.prop dictionary
        # Read current_optics.prop dictionary from pickle file
        # self.optics.prop = pkl.load(open(self.pklfile,'rb'))
        #self.optics = Optics(self.settings)
        #self.optics.prop = pkl.load(open(self.pklfile,'rb'))
        self.optics.set_opt_depth_species(atm)
        
        rad_lbl, dev_tau_lbl, dev_alb_lbl = self.transmission(surface, mu0, muv)
       
        fwd={}
        fwd['rad']  = slit.dot(rad_lbl)
        fwd['alb0'] = slit.dot(dev_alb_lbl)
        fwd['alb1'] = slit.dot( dev_alb_lbl*self.settings.wave_alb)
        fwd['alb2'] = slit.dot(dev_alb_lbl*self.settings.wave_alb**2)
        fwd['alb3'] = slit.dot(dev_alb_lbl*self.settings.wave_alb**3)
        fwd['alb4'] = slit.dot(dev_alb_lbl*self.settings.wave_alb**4)
        fwd['alb5'] = slit.dot(dev_alb_lbl*self.settings.wave_alb**5)
        fwd['alb6'] = slit.dot(dev_alb_lbl*self.settings.wave_alb**6)
        fwd['alb7'] = slit.dot(dev_alb_lbl*self.settings.wave_alb**7)
        fwd['alb8'] = slit.dot(dev_alb_lbl*self.settings.wave_alb**8)
        fwd['alb9'] = slit.dot(dev_alb_lbl*self.settings.wave_alb**9)
        for spec in self.settings.species:
            fwd[spec] =slit.dot(np.sum(self.optics.prop[spec]['tau_per_molec'], axis=1)*dev_tau_lbl)
        
        return fwd 
    
    def run(self,i,j,x_dict0,sy,nu=0,y=None,normalization=False,method="LSQ",maxc=1E99,plot=False,save=False):
        """
        Least square solution for non-linear 
        forward model included in a Gauss Newton technique
        
        Parameters 
        ----------
        nu (int): Initial guess for wavelength shift. Default is 0 
        plot: boolean 
              if plot=True, noise from the measurement and residual
              between measurement and retrieved model are plotted in one figure, 
              the initial guess, the retrieved model and the measurement 
               are plotted in another figure. 
               default value: False 
               
         save: boolean 
               in case save=True plot should be set to True. The plot will be saved in the file where the 
                           jupyter notebook is located. 
                           default value: False
                           
        Returns
        -------             
        LSQ: dictionary for Least-Square Solution
        LSQ["x_lst"]: retrieved least square solution if konvergence is reached
        LSQ["Chi_sqrt"]: last Chi squared value
        LSQ["Sx"]: Error covariance matrix of least square solution
        LSQ["residuum"]: y-retrieved spectra 
        LSQ["Fn"]: last retrieved spectra 
        LSQ["wave_meas"]: retrieved wavelength 
        conc: dictionary, caculated concentration of molecules
        conc["CH4"]: concentration CH4  ppb
        conc["CO2"]: concentration CO2  ppm
        conc["H2O"]: concentration H2O  %
        conc_err["CH4"]: concentration error (gauss error progagation) CH4  ppb
        conc_err["CO2"]: concentration error (gauss error progagation) CO2  ppb
        conc_err["H2O"]: concentration error (gauss error progagation) H2O  ppb
                      
                        
        """
        global test
        atm = self.atm.atmo.copy()
        surf = self.surface
        meas = self.measurement.read_netcdf4(i,j)
        #nu = self.settings.nu
        
        x_dict = x_dict0.copy()
        if y is None: 
            y=meas["spec"]#/np.max(meas["spec"])
        if normalization:
            y=meas["spec"]/np.max(meas["spec"])
        
        # if test:
        #     out=[]
        #     print(atm)
        #     print(test[0])
        #     out.append(atm==test[0])
        #     out.append(all(surf==test[1]))
        #     out.append(all(meas==test[2]))
        #     out.append(all(nu==test[3]))
        #     out.append(all(x_dict==test[4]))
        print("Same as last measurement?", test==x_dict)
        test=x_dict
        #plt.plot(y)
        mu0=meas["mu0"]
        muv=meas["muv"]
        
        
        
        
        yerr=y*0.01
        # to be able to calculate the derivative with a small perturbation, it is necessary to work with 64 float self.settings.numbers. 
        # meas["wave"] is saved as a reference wavelength grid which is than shifted by self.settings.nu. 
        wave_meas_ref=np.array([meas["wave"]],dtype=np.float64)[0]
        fwhm = meas["fwhm"]
        k_dic={}   # dictionary of the state vector, later inserted in nonscattered forward model
        for key in x_dict.keys():
            if (key=='CO2'):
                k_dic[key]=key
               
            if (key=='CH4'):
                k_dic[key]=key
               
            if (key=='H2O'):
                k_dic[key]=key
                
            if (key=='N2O'):
                k_dic[key]=key
                
            if (key=='NH3'):
                k_dic[key]=key
                
        m=0
        while "alb%d"%(m)in x_dict.keys():
                k_dic["alb%d"%(m)]="alb%d"%(m)
                m=m+1 
           
        # start with iteration
        iteration=0
        konvergence=False
        Chi_sqrt=[]  
        alb_lst=[]
        for key in x_dict.keys():
            if (key=='CO2'):
                atm["CO2"]=x_dict["CO2"]*self.ref_CO2
            if (key=='CH4'):
               atm["CH4"]=x_dict["CH4"]*self.ref_CH4
            if (key=='H2O'):
               atm["H2O"]=x_dict["H2O"]*self.ref_H2O
            if (key=='N2O'):
               atm["N2O"]=x_dict["N2O"]*self.ref_N2O
        m=0    
        while "alb%d"%(m)in x_dict.keys():
                alb_lst.append(x_dict["alb%d"%(m)])
                m=m+1        
        
        #shift wave_meas_ref by initial guess of nu
        wave_meas=wave_meas_ref-nu
        
        x0_lst=[]
        for value in x_dict.values():
            x0_lst.append(value)
        x0_lst.append(nu)
        x0_lst=np.asarray(x0_lst)
        
        xapr=x0_lst
        # Calculate surface data
        surf.get_albedo_poly(alb_lst)
        
        
        wave_lbl64=np.array([self.settings.wave_lbl],dtype=np.float64)[0]
        slit=slit_conv(fwhm, wave_lbl64, wave_meas)
        
        # Calculate nonscattered forward model
        F = self.nonscat_fwd_model(slit,atm, surf, mu0, muv, fwhm)
        ytilde = y - F["rad"]  # Difference between forward model and measured spectrum
        F_0=F
        #Chi_sqrt.append(ytilde.T@ np.linalg.inv(sy)@  ytilde/(len(F["rad"])-(len(x0_lst))))
        try:
            Chi_sqrt.append(ytilde.T@ np.linalg.pinv(sy)@  ytilde/(len(F["rad"])-(len(x0_lst))))
        except: print("no chi calculation possible") 
        ###################################################################################################
        iteration=1
        while (konvergence==False):
            try: #Try Least square solution
            
                # Little perturbation of wavelength to calculate the derivative of the wavelenght shift 
                wave_dash=wave_meas+0.000001
                slit_per=slit_conv(fwhm, self.settings.wave_lbl, wave_dash)
                F_dash = self.nonscat_fwd_model(slit_per,atm, surf, mu0, muv, fwhm)#
                # Calculate convoluted Jacobians
                K= np.zeros([len(F["rad"]),len(x_dict)+1])
            
                l=0
                for value in k_dic.values():
                        K[:,l]=F[value]
                        l+=1
                # Derivative of the wavelenght shift 
                K[:,l]=(F_dash["rad"]-F["rad"])/0.000001
                
                if method=="LSQ":
                    
                    # Calculated least square solution   
                    syinv = np.linalg.pinv(sy)                    # inverse of covariance matrix of the measurement
                    #print(self.settings.syinv)
                    #print(np.dot(self.settings.syinv,K))
                    Sx = np.linalg.pinv(np.dot(K.T,np.dot(syinv,K))) # covariance matrix of the estimated least square solution
                    G  = np.dot(Sx,np.dot(K.T,syinv))               # gain matrix
                    x  = np.dot(G,y-F["rad"])                   # least square solution
                    
                if method=="SVD":
                    m = len(y) # dim. measurements
                    n = len(x_dict)+1 # dim. state vector
                    
                    # Consider measurement uncertainties by merging them into K and y,
                    # advantage: we do not need to carry sy along 
                    Kw = np.zeros([m,n])
                    yw = np.zeros([m])
                    Fw =np.zeros([m])
                    
                    for i in range(m):
                        Kw[i,:] = K[i,:]/yerr[i] 
                        yw[i]   = y[i]/yerr[i]
                        Fw[i]=F["rad"][i]/yerr[i]
                    
                    # Singular value decomposition of forward model
                    # By putting full_matrices = False, the obvious self.settings.null-space due to
                    # m > n is removed
                    U,Svec,VT = np.linalg.svd(Kw,full_matrices=False)
                    UT = U.T # U transpose
                    V = VT.T # V
                    Sinv = np.array([1./s for s in Svec]) # Inverse of singular values
                    
                    
                    # Calculate condition self.settings.number and truncate all contributions
                    # producing condition self.settings.numbers > maxc
                    CNvec = np.array([max(Svec)/s for s in Svec]) # Condition self.settings.number per singular value
                    mask = CNvec <= maxc # Mask: True if condition self.settings.number is smaller than limit, False otherwise
                    Sinv_trunc = np.where(mask == False, 0, Sinv) # Set all inverse singular values to zero if mask False
                    
                    # Moore Penrose pseudo-inverse
                    G = V.dot(np.diag(Sinv_trunc).dot(UT)) # a.k.a. gain matrix
                    # Estimates
                    x = G.dot(yw-Fw) # state estimate
                    Sx = G.dot(G.T)# a posteriori error covariance, remember: errors were merged into K 
                    A = G.dot(Kw) # averaging kernel, averaging kernel could also be calculated via VVT, 
                                  # but then, we would need to truncate V
                if method=="MAP":                  
                    n = len(x_dict)+1 # dim. state vector
                    m = len(y) # dim. measurements
                    
                    # Consider measurement uncertainties by merging them into F, K, y, calculate sy for error propagation
                    Fw = np.zeros([m])
                    Kw = np.zeros([m,n])
                    yw = np.zeros([m])
                    sy = np.eye(m)
                    syinv = np.eye(m)
                    for i in range(m):
                        Fw[i]     = F["rad"][i]/yerr[i]
                        Kw[i,:]   = K[i,:]/yerr[i] 
                        yw[i]     = y[i]/yerr[i]
                        sy[i,i]   = yerr[i]**2
                        syinv[i,i]= 1./yerr[i]**2
                    
                    Sapr = np.eye(n)
                    # Blockwise: 20 x CO2, 20 x H2O, 3 x albedo
                    # CO2 parameters (in 20 layers)
                    #for i in range(0):
                    Sapr[1,1]=Sapr[1,1]*(1E-1)**2 # 1% standard deviation
                    Sapr[2,2]=Sapr[1,1]*(1E+1)**2 # 1% standard deviation
                    # H2O parameters (in 20 layers)
                    #for i in range(1):
                    Sapr[0,0]=Sapr[0,0]*(1E+1)**2 # 5% standard deviation
                    # Albedo parameters (3 polynomial coefficients)
                    Sapr[3:n,3:n]=Sapr[3:n,3:n]*1E3**2 # infinity standard deviation
                    # Inverse of a priori covariance
                    Saprinv= np.linalg.pinv(Sapr)  
                    # Posterior error matrix with optional Levenberg-Marquardt: 
                    # (KT.K + Sapr-1 + gLM*1)-1
                    gLM=0
                    gRS=0
                    if iteration==0:
                        gRS=1
                    SLM = np.linalg.pinv(Kw.T.dot(Kw) + (1.+gLM)*Saprinv)
                    # Gain matrix / pseudo-inverse: S.KT
                    Gw = SLM.dot(Kw.T)
                    # State update dx = G.(y-F) - S.Sapr-1.(xn - xapr)
                    dx = Gw.dot(yw-Fw) - SLM.dot( Saprinv.dot(x0_lst-xapr) )
                    # New state, with optional stepsize reduction: xn+1 = xn + (1./1.+gRS) * dx
                    x = (1./(1.+gRS)) * dx
                    
                    # Posterior noise plus smoothing covariance (without LM component)
                    # (KT.K + Sapr-1)-1
                    Sx = np.linalg.pinv(K.T.dot(syinv.dot(K)) + Saprinv)
                    # Gain matrix / pseudo-inverse to be applied to y (not yw)
                    G = Sx.dot(K.T.dot(syinv))
                    # Averaging kernel: G.K
                    A = G.dot(K) 
                    # State noise covariance G.GT
                    Snoise = G.dot(sy.dot(G.T))
                
                Sx_dict={}
                m=0
                for key in x_dict.keys():
                    Sx_dict[key]=Sx[m,m]**0.5
                    m+=1
    
                ylsq = K.dot(x)
                #plt.show()
                #plt.plot(wave_meas,y-F["rad"]-ylsq)
                #plt.show()
                #ylsqw = G.dot(x)
                k=0
                x_lst=[]
                x_dict_1=x_dict.copy()
                for key in x_dict.keys():
                    if key=="CO2" or key =="CH4" or key=="H2O":
                        x_dict[key]=np.abs(x0_lst[k]+x[k])
                        x_lst.append(np.abs(x0_lst[k]+x[k]))
                        k+=1
                    else:# remember that we estimated x-x0, calculate x
                        x_dict[key]=x0_lst[k]+x[k]
                        x_lst.append(x0_lst[k]+x[k])
                        k+=1
                
                nu=x0_lst[k]-x[k] #TODO: Why "-"?
                #print(x_dict)
                # Update state vector
                alb_lst=[]
                for key in x_dict.keys():
                    if (key=='CO2'):
                        atm["CO2"]=x_dict["CO2"]*self.ref_CO2
                    if (key=='CH4'):
                       atm["CH4"]=x_dict["CH4"]*self.ref_CH4
                    if (key=='H2O'):
                       atm["H2O"]=x_dict["H2O"]*self.ref_H2O
                    if (key=='N2O'):
                       atm["N2O"]=x_dict["N2O"]*self.ref_N2O
                    
                
                #print(x_dict)  
                m=0    
                while "alb%d"%(m)in x_dict.keys():
                        alb_lst.append(x_dict["alb%d"%(m)])
                        m=m+1 
                wave_meas=wave_meas_ref-nu
                #print(x0_lst)
                x0_lst=[]
                for value in x_dict.values():
                    x0_lst.append(value)
                x0_lst.append(nu)
                
                x0_lst=np.asarray(x0_lst)
                
                # Calculate updated surface data
                surf.get_albedo_poly(alb_lst)
               
                slit=slit_conv(fwhm, self.settings.wave_lbl, wave_meas)
                # Calculate nonscattered forward model with updated state vector
                F = self.nonscat_fwd_model(slit,atm, surf, mu0, muv, fwhm)
                ytilde = y - F["rad"]  # Difference between forward model and measured spectrum
                
                # Calculate Chi squared for updated state vector
                Chi_sqrt.append(ytilde.T@ np.linalg.pinv(sy)@  ytilde/(len(F["rad"])-(len(x0_lst))))
                
                    
                # Define convergence criteria or if iterations are too large for second and higher retrieved spectrum
                if iteration!=0:
                    
                    if ma.is_masked(Chi_sqrt[iteration]):
                        warn("Point is masked",stacklevel=1)
                        break
                    
                    #if (np.abs(Chi_sqrt[iteration]-Chi_sqrt[iteration-1])<1 or iteration>=10):
                    #    if Chi_sqrt[iteration]<Chi_sqrt[iteration-1]: # Use better solution if oscillations occured
                    #        konvergence=True
                    if np.abs(x_dict_1["CH4"]-x_dict["CH4"]) < 0.04 or iteration>=10:
                        konvergence = True
                        if iteration==15: # no convergence happened -> error =True 
                            warn("No convergence in 55 iterations.",stacklevel=1)
                            break
                    #else: konvergence=False
                            
                
                iteration=iteration+1
                if np.isnan(Chi_sqrt[iteration-1]):
                    warn("Data corrupted",stacklevel=1)
                    break
            
            except:# if error in Calculation of Least Square solution (sy,x,...) # No convergence, therefore error=true
                warn("No convergence try",stacklevel=1)
                break 
                #iteration=1
            #error=True
                
        #plt.show() 
        if (plot and konvergence): 
        	
    
            #import seaborn as sns
    
            #sns.set_theme()
            plt.figure(figsize=(4,3.5))
            plt.plot(wave_meas,y,"k",label="measurement")
            #plt.plot(wave_meas,F_0['rad'],"skyblue",label='F$_0$') #+ylsq
            plt.plot(wave_meas,F["rad"],"g",label='PRIM',lw=2)#,label='F$_n$: CO$_2$ = (464$\pm$ 5)ppm \n $\chi^2_{red}$: 5.4')#label='F$_n$: ch4=(1910$\pm$120) ppb \n $\chi^2_{red}$=7.9')#,label='F$_n$: ch4 = (464$\pm$ 5)ppm \n $\chi^2_{red}$: 5.4') #+ylsq
            plt.xticks([2100,2200,2300,2400],fontsize=13)
            plt.yticks(fontsize=13)
            plt.ylabel("radiance \n [W m$^{-2}\,\mu$m$^{-1}$ sr$^{-1}$]",fontsize=16)
            plt.xlabel("wavelength [nm]",fontsize=16)
            plt.xlim(2115,2450)
            plt.title("retrieved and measured spectrum (CH$_4$)",fontsize=18)
            #plt.xlim(1975,2100)
            plt.legend(fontsize=13)
            if save==True:
                plt.savefig("spectrum.png",transparent=True,dpi=550,bbox_inches="tight")#,transparent=True
            plt.show()
            
        
            plt.plot(wave_meas,y-F["rad"],label='residual',color="red") #(ylsq+
            plt.fill_between(wave_meas,y*0.01,-y*0.01,color="blue",alpha=0.3,label="noise")
            #plt.plot(wave_meas,-y*0.01,"b-")
            plt.ylabel("W/(m$^2$ nm)",fontsize=13)
            plt.xlabel("wavelength nm",fontsize=13)
            plt.legend()
            if save==True:
                plt.savefig("noise_residual.png")
            plt.show()
        
         
        print("Did it konverge?",konvergence)  
        #print(",how much iterations?",iteration)    
        print(",Chi_sqrt:",Chi_sqrt)   
        print(x_dict,nu)                
        
        conc={}
        conc_err={}
        if not konvergence: # no convergence => np.nan (with right dimensions for x_lst and Sx)
            x0_lst=[]
            for value in x_dict.values():
                x0_lst.append(value)
            x_lst=np.asarray(x0_lst)
            Sx=np.zeros([len(x0_lst),len(x0_lst)])
            LSQ={}
            LSQ["x_lst"]=x_lst
            LSQ["Sx"]=Sx
            LSQ["Chi_sqrt"]=np.nan
            LSQ["Fn"]=F_0
            LSQ["residuum"]=y-F["rad"]
            LSQ["wave_meas"]=np.nan
            LSQ["nu"]=np.nan
            LSQ["F"]=F['rad']
            #plt.plot(wave_meas,F['rad'],label='F_ini')
            #plt.plot(wave_meas,y,label="y [Satellite]")
            #plt.show()
            #plt.plot(wave_meas,y-F["rad"],label='residual',color="red") #(ylsq+
            #plt.fill_between(wave_meas,y*0.01,-y*0.01,color="blue",alpha=0.3,label="noise")
            conc["AIR"] = np.nan
            conc["AIR_0"] = np.nan
            for key in x_dict.keys():
                if (key=='CO2'):
                    conc[key+"ppm"]=np.nan
                    conc_err[key+"ppm"]=np.nan
                    conc[key + "ppm/m^2"] = np.nan
                    conc_err[key + "ppm/m^2"] = np.nan
                if (key=='H2O'):
                    conc[key +"pr"]=np.nan
                    conc_err[key +"pr"]=np.nan
                    conc[key + "pr/m^2"]= np.nan
                    conc_err[key + "pr/m^2"] = np.nan
                if (key=='CH4'):
                    conc[key + "ppb"]= np.nan
                    conc_err[key + "ppb"]=np.nan
                    conc[key + "ppb/m^2"] = np.nan
                    conc_err[key + "ppb/m^2"] = np.nan
                if (key=='N2O'):
                    conc[key]=np.nan
                    conc_err[key]=np.nan
        if konvergence:  
            LSQ={}
            LSQ["x_lst"]=x_lst
            LSQ["Sx"]=Sx
            LSQ["Chi_sqrt"]=Chi_sqrt[iteration-1]
            LSQ["Fn"]=F['rad']#+ylsq
            LSQ["residuum"]=y-F['rad']#(ylsq+
            LSQ["wave_meas"]=wave_meas
            LSQ["nu"]=nu
            conc["AIR"] = sum(atm['AIR'])
            conc["AIR_0"]=atm['AIR'][19]
            for key in x_dict.keys():
                if (key=='CO2'):
                    atm["CO2"]=x_dict["CO2"]*self.ref_CO2
                    co2_conc= sum(atm['CO2'])/sum(atm['AIR'])*1E6 #ppm
                    co2_err= Sx_dict[key]/x_dict["CO2"]*co2_conc
    
                    #sum(Sx_lst[0]/x_dict["CO2"]*ref_CO2)/(sum(atm['AIR']))*1E6
                    conc[key+"ppm"]=co2_conc
                    conc_err[key+"ppm"]=co2_err
                    conc[key + "ppm/m^2"] = sum(atm['CO2'])*1E6
                    conc_err[key + "ppm/m^2"] = sum(atm['CO2'])*1E6*Sx_dict[key]/x_dict["CO2"]
            
                if (key=='CH4'):
                   atm["CH4"]=x_dict["CH4"]*self.ref_CH4
                   ch4_conc= sum(atm['CH4'])/sum(atm['AIR'])*1E9 #ppb
                   conc["CH4withoutAIR"]=x_dict["CH4"]*self.ref_CH4/atm['AIR']
                   ch4_err= Sx_dict[key]/x_dict["CH4"]*ch4_conc #sum(atm['CH4']/(sum(atm['AIR']))*Sx_lst[0]*1E9)
                   conc[key + "ppb"]=ch4_conc
                   conc_err[key + "ppb"]=ch4_err
                   conc[key + "ppb/m^2"] = sum(atm['CH4'])*1E9
                   conc_err[key + "ppb/m^2"] = sum(atm['CH4'])*1E9*Sx_dict[key]/x_dict["CH4"]
            
                if (key=='H2O'):
                   atm["H2O"]=x_dict["H2O"]*self.ref_H2O
                   h2o_conc= sum(atm['H2O'])/sum(atm['AIR'])*1E2 # %
                   h2o_err= Sx_dict[key]/x_dict["H2O"]*h2o_conc #sum(atm['H2O']/(sum(atm['AIR']))*Sx_lst[1]*1E2) # %
                   conc[key +"pr"]=h2o_conc
                   conc_err[key +"pr"]=h2o_err
                   conc[key + "pr/m^2"]= sum(atm['H2O'])*1E2
                   conc_err[key + "pr/m^2"] = sum(atm['H2O'])*1E9*Sx_dict[key]/x_dict["H2O"]
                
                if (key=='N2O'):
                   atm["N2O"]=x_dict["N2O"]*self.ref_N2O
                   n2o_conc= sum(atm['N2O'])/sum(atm['AIR'])*1E9 # %
                   n2o_err= Sx_dict[key]/x_dict["N2O"]*n2o_conc #sum(atm['H2O']/(sum(atm['AIR']))*Sx_lst[1]*1E2) # %
                   conc[key]=n2o_conc
                   conc_err[key]=n2o_err
        return LSQ,conc,conc_err #TODO: build dictionary
class PRIM(PRISCO):
    def __init__(self,settings,location):
        super().__init__(self,settings,location)
    def run(self,i,j,x_dict0,y=None,normalization=False,method="LSQ",maxc=1E99,plot=False,save=False):
        
            
        sy = np.eye(len(y))*(0.01*y)**2
        super().run(self,i,j,x_dict0,sy,y=y,normalization=normalization,method=method,maxc=maxc,plot=plot,save=save)
        

        
        
        
        
        
        
        
        
        
        
