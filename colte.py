def colte(sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,bprp_ex=False,pmod=False,outfile=False,MC=False,trials=False,wato=False,elogg=[],efeh=[],egg=[],ebp=[],erp=[],ej2=[],eh2=[],ek2=[],eebv=[]):
    
    '''
 PURPOSE:
     Compute stellar effective temperatures using colour-Teff relations in the
     Gaia and 2MASS photometric system

 EXPLANATION:
     The relations used to derive Teff are from Casagrande+2020.
     For each star, Teffs are computed from up to 12 different colour indices
     and results are written into a csv file.
     If the option for a MonteCarlo is set, Teff uncertainties are computed
     for each colour index, and a final weighted average for Teff along with
     its weighted standard deviation is derived. 
     The routine applies a few bare quality cuts on input photometry
     by removing BP and RP<5, G<6, J<5.0, H<4.8, K<4.2, and if available 
     also removing ej2>0.05, eh2>0.05, ek2>0.05. These cuts are mostly to
     avoid uncertainties due to saturation at bright magnitudes.  
     Further quality cuts on Gaia photometry can be set with input parameters
     bprp_ex= and pmod=    
     Also, stars with ebv<0, logg < 0 or > 5, or feh < -4 or > 0.6
     will be excluded. 

 REQUIRED INPUT PARAMETERS
 sid:     star name/ID
 logg:    surface gravity
 feh:     [Fe/H]
 gg:      Gaia G  (phot_g_mean_mag) 
 bp:      Gaia BP (phot_bp_mean_mag)
 rp:      Gaia RP (phot_rp_mean_mag)
 j2:      2MASS J
 h2:      2MASS H
 k2:      2MASS K
 ebv:     Reddening E(B-V)

 For each star, bp,rp,logg,feh,ebv are indispensable parameters needed 
 to derive Teff from at least bp-rp. Note that sid,gg,j2,h2,k2 are also 
 required inputs, but empty entries can be passed if some of these quantities
 are unavaible for a star.
  
 OPTIONAL INPUT PARAMETERS
 bprp_ex: to remove stars with bad phot_bp_rp_excess_factor (Eq 2, Arenou+18)
 pmod:    to retain only stars with phot_proc_mode=0  (Riello+18)
 outfile: output file. If not passed, then the default output file is colte.csv
 MC:      to perform a MonteCarlo of Teff uncertainties in different bands
 ej2:     2MASS J uncertainty. If not provided, 0.022 mag is assumed
 eh2:     2MASS H uncertainty. If not provided, 0.024 mag is assumed
 ek2:     2MASS K uncertainty. If not provided, 0.022 mag is assumed

 OPTIONAL INPUT PARAMETERS relevant ONLY if MC=True
 trials:  number of MC realizations for each star.If not set, default is 1000
          Default value is a good compromise between speed of execution and 
          convergence. The latter depends on the colour index, but as a rule of
          thumb, with 1000 trials, uncertainties typically converge to within 
          a few K, or ~10K in worst cases. With 100 trials, convergence is
          ~10K in most cases, and up to ~70K in worst cases. With 10000 trials
          convergence is always within a few K
 wato:    to write Weighted Averaged Teff Only in the output file 
 elogg:   logg uncertainty.    If not provided, 0.2 dex is assumed
 efeh:    [Fe/H] uncertainty.  If not provided, 0.1 dex is assumed
 egg:     Gaia G uncertainty.  If not provided, 0.005 mag is assumed
 ebp:     Gaia BP uncertainty. If not provided, 0.005 mag is assumed
 erp:     Gaia RP uncertainty. If not provided, 0.005 mag is assumed
 eebv:    Reddening uncertainty. If not provided, 10% of input ebv is assumed

 OUTPUT 
 The routine will write an output file providing for each star the adopted
 sid, logg, feh, ebv + Teffs computed from up to 12 colour indices. If Teff
 cannot be determined in a colour index, NaN is returned for that index. 
 Note that the program makes a number of basic quality cuts on input
 photometry, and requires a value for logg, feh and ebv. Hence, the output
 file might contain fewer stars than the input file.
 If MC is set, then an uncertainty is provided for each Teff, along with
 weighted averaged Teff and weighted standard deviation. If WATO is set,
 only weighted average and weighted standard deviation are written. Note that
 weighted averaged Teff and weighted standard deviation might change by a few
 Kelvin each time, because of the MC nature of the errors (more robust
 convergence can be achieved by increasing trials). 

 HISTORY
 Written Nov 5, 2020 by Luca Casagrande

    '''
    import numpy as np

    # remove warning messages arising when np.where encounters NaN 
    import warnings
    warnings.simplefilter(action = "ignore", category = RuntimeWarning)
    
    ppm=np.zeros(len(bp))

    '''
  ; First cut to retain only usable stars. Further cuts applied afterwards.
  ; -Gaia BP and RP needed for all stars. Finite function to remove NaN and Inf
  ;  in case there is any.  
  ; -Only BP and RP > 5 are considered because of uncalibrated systematics
  ;  at brighter magnitudes.
  ; -Stars without a value of reddening, logg or feh are excluded.
  ; -Stars with ebv<0, logg < 0 or > 5, feh < -4 or > 0.6 are also excluded.
  ; -If bprp_ex and/or pmod keywords are set, then stars with bad
  ;  phot_bp_rp_excess_factor and/or phot_proc_mode are also excluded  
    '''
    
    if pmod== True:
        ppm=pmod
    if bprp_ex==True:
        ok=np.where((np.isfinite(bp)==True) & (np.isfinite(rp)==True)      & \
                    (bp>=5) & (rp>=5) & (np.isfinite(ebv)==True)           & \
                    (np.isfinite(logg)==True) & (np.isfinite(feh)==True)   & \
                    (ebv>=0) & (logg>=0) & (logg<=5)&(feh>=-4)& (feh<=0.6) & \
                    (bprp_ex< (1.3+0.06*(bp-rp)*(bp-rp)) )                 & \
                    ( bprp_ex>(1.+0.015*(bp-rp)*(bp-rp)))                  & \
                    (ppm==0));ntot=len(ok[0])
    else:
        ok=np.where((np.isfinite(bp)==True) & (np.isfinite(rp)==True)    & \
                    (bp>=5) & (rp>=5) & (np.isfinite(ebv)==True)         & \
                    (np.isfinite(logg)==True) & (np.isfinite(feh)==True) & \
                    (ebv>=0) & (logg>=0)&(logg<=5)&(feh>=-4)&(feh<=0.6)  & \
                    (ppm==0));ntot=len(ok[0])


    if ntot >=1:
        sid    = [str(i) for i in sid[ok]]
        gg     = gg[ok]; bp     = bp[ok]; rp     = rp[ok]
        j2     = j2[ok]; h2     = h2[ok]; k2     = k2[ok]
        ebv    = ebv[ok];logg   = logg[ok];feh   = feh[ok]

        fatto=0

        #write output file with results
        if type(outfile)==str:
            f=open(outfile,'w')
        else: f=open('colte.csv','w')

        #things done if MC is set
        if MC==True:
            if wato==True: f.write('#star_name,logg,feh,ebv,wTeff,ewTeff  \n')
            else: f.write('#star_name,logg,feh,ebv,T_bprp,eT_bprp,T_bpj,eT_bpj, T_bph,eT_bph,T_bpk,eT_bpk,T_rpj,eT_rpj,T_rph,eT_rph,T_rpk,eT_rp,T_gj,eT_gj,T_gh,eT_gh,T_gk,eT_gk,T_gbp,eT_gbp,T_grp,eT_grp,wTeff,ewTeff \n')

            #number of MC realizations
            if trials==False: trials=1000
            
            dcol =np.zeros([12,trials])
            dteff=np.zeros([12,trials])
            eteff=np.zeros(12)

            #if errors are given then retain them. If errors are not given or
            #set to zero, then reasonable default values are assigned. 

            if len(elogg)!=0: #errors only arrays or false
                nbad   = 0
                elogg  = elogg[ok]
                bad=np.where((elogg<=0) | (np.isfinite(elogg)==False));nbad=len(bad[0])
                if nbad>=1: elogg[bad]=0.2
            else:
                elogg=np.zeros(ntot)+0.2

            if len(efeh)!=0:
                nbad   = 0
                efeh  = efeh[ok]
                bad=np.where((efeh<=0) | (np.isfinite(efeh)==False));nbad=len(bad[0])
                if nbad>=1: efeh[bad]=0.1
            else:
                efeh=np.zeros(ntot)+0.1

            if len(egg)!=0:
                nbad   = 0
                egg  = egg[ok]
                bad=np.where((egg<=0) | (np.isfinite(egg)==False));nbad=len(bad[0])
                if nbad>=1: egg[bad]=0.005
            else:
                egg=np.zeros(ntot)+0.005

            if len(ebp)!=0:
                nbad   = 0
                ebp  = ebp[ok]
                bad=np.where((ebp<=0) | (np.isfinite(ebp)==False));nbad=len(bad[0])
                if nbad>=1: ebp[bad]=0.005
            else:
                ebp=np.zeros(ntot)+0.005

            if len(erp)!=0:
                nbad   = 0
                ebp  = erp[ok]
                bad=np.where((erp<=0) | (np.isfinite(erp)==False));nbad=len(bad[0])
                if nbad>=1: erp[bad]=0.005
            else:
                erp=np.zeros(ntot)+0.005

            #default error for reddening is 10%
            if len(eebv)!=0:
                nbad   = 0
                eebv  = eebv[ok]
                bad=np.where((eebv<=0) | (np.isfinite(eebv)==False));nbad=len(bad[0])
                if nbad>=1: eebv[bad]=0.1*ebv[bad]
            else:
                eebv=0.1*ebv

        else:
            f.write('#star_name,logg,feh,ebv,T_bprp,T_bpj,T_bph,T_bpk,T_rpj,T_rph,T_rpk,T_gj,T_gh,T_gk,T_gbp,T_grp \n')

        #as above, retain or assign resonable 2MASS errors
        if len(ej2)!=0:
            nbad   = 0
            ej2  = ej2[ok]
            bad=np.where((ej2<=0) | (np.isfinite(ej2)==False));nbad=len(bad[0])
            if nbad>=1: ej2[bad]=0.022
        else:
            ej2=np.zeros(ntot)+0.022

        if len(eh2)!=0:
            nbad   = 0
            eh2  = eh2[ok]
            bad=np.where((eh2<=0) | (np.isfinite(eh2)==False));nbad=len(bad[0])
            if nbad>=1: eh2[bad]=0.024
        else:
            eh2=np.zeros(ntot)+0.024

        if len(ek2)!=0:
            nbad   = 0
            ek2  = ek2[ok]
            bad=np.where((ek2<=0) | (np.isfinite(ek2)==False));nbad=len(bad[0])
            if nbad>=1: ek2[bad]=0.022
        else:
            ek2=np.zeros(ntot)+0.022
        
        #remove Gaia G<6 because of uncalibrated CCD saturation
        nog=np.where(gg<6); nnog=len(nog[0])
        if nnog>=1: gg[nog]=float("nan")

        #correct Gaia G as per Maiz Apellaniz & Weiler (2018)
        gp=gg-0.0032*(gg-6)
        
        #more house-cleaning to remove 2MASS photometry likely to be bad
        noj=np.where((j2<5)  |(ej2>0.05)); nnoj=len(noj[0])
        noh=np.where((h2<4.8)|(eh2>0.05)); nnoh=len(noh[0])
        nok=np.where((k2<4.2)|(ek2>0.05)); nnok=len(nok[0])
        if nnoj>=1: j2[noj]=float("nan")
        if nnoj>=1: h2[noh]=float("nan")
        if nnok>=1: k2[nok]=float("nan")
        
        #compute colour dependent extinction coefficients
        bprp0 = (bp-rp)- 1.339*ebv
        R_gg  = 3.068  - 0.505*bprp0 + 0.053*bprp0*bprp0
        R_bp  = 3.533  - 0.114*bprp0 - 0.219*bprp0*bprp0 + 0.070*bprp0*bprp0*bprp0
        R_rp  = 2.078  - 0.073*bprp0
        R_j2  = np.zeros(ntot) + 0.899
        R_h2  = np.zeros(ntot) + 0.567
        R_k2  = np.zeros(ntot) + 0.366
        
        #colour-teff coefficients for Gaia DR2
        cpol=np.zeros([12,15])
        cpol[0]  = np.array([7928.0505, -3663.1140,  803.3017,    -9.3727,        0., 325.1324, -500.1160,  279.4832,  -53.5062,      0.,   -2.4205, -128.0354, 49.4933,  5.9146,  41.3650]) #B P-RP
        cpol[1]  = np.array([8217.8748, -2526.8430,  458.1827,   -28.4540,        0., 234.0113, -205.3084,   63.4781,   -7.2083,      0.,  -85.7048,  -50.1557, 32.3428, -2.3553,  20.0671]) # BP-J
        cpol[2]  = np.array([8462.0737, -2570.3684,  537.5968,   -44.3644,        0., 189.1198, -106.7584,   31.1720,   -4.9137,      0.,   -9.2587, -189.8600, 75.8619, -6.8592,  16.7226]) # BP-H
        cpol[3]  = np.array([8404.4760, -2265.1355,  403.4693,   -27.9056,        0., 193.5820, -145.3724,   47.7998,   -6.4572,      0.,  -34.5438, -130.2559, 52.6470, -4.4777,  15.8249]) # BP-K

        cpol[4]  = np.array([9073.7917, -7670.6606, 3164.0525,         0., -126.1476,       0.,   -7.3816,  -12.5168,        0., -2.0452,        0.,   76.1144,      0.,      0., -45.8056]) # RP-J
        cpol[5]  = np.array([8924.1553, -4779.3394, 1319.8989,         0.,  -16.6676,       0.,  -23.6583,   22.4243,        0., -4.3066,        0.,   35.0102,      0.,      0., -28.7228]) # RP-H
        cpol[6]  = np.array([8940.4628, -4450.6138, 1138.6816,         0.,  -10.5749,       0.,  -42.3037,   33.3365,        0., -3.2535,        0.,   41.0402,      0.,      0., -21.9922]) # RP-K

        cpol[7]  = np.array([8369.9905, -3559.7710,  895.8869,   -86.7011,        0., 180.7568, -164.9264,   24.4263,    4.2318,      0., -127.9640,   72.1449,      0.,      0.,  13.7683]) #  G-J
        cpol[8]  = np.array([8185.8827, -2536.7671,  503.2762,   -42.7871,        0., 230.4871, -254.5291,  104.6258,  -17.4859,      0., -122.0732,   45.0572,      0.,      0.,   6.9992]) #  G-H
        cpol[9]  = np.array([8103.2039, -1857.7194,        0.,    73.1834,   -1.7576, 236.0335, -345.9070,  170.4915,  -28.8549,      0., -131.4548,   49.6232,      0.,      0.,  10.0777]) #  G-K

        cpol[10] = np.array([7555.3516,  5803.7715,        0., -2441.7124,  437.7314, 455.0997, 2243.1333, 3669.4924, 1872.7035,      0.,   19.1085,   75.2198,      0.,      0., -83.9777]) # G-BP 
        cpol[11] = np.array([7971.3782, -5737.5049,        0.,  1619.9946, -203.8234, 255.7408, -492.8268,  160.1957,  103.1114,      0.,  -64.3289,   34.3339,      0.,      0.,  54.7224]) # G-RP

        # colour range for dwarfs
        d_r=np.array([2.00,3.00,4.00,4.20,1.05,1.60,1.85,2.10,2.60,2.80,-0.10,0.85])
        d_b=np.array([0.15,0.25,0.40,0.30,0.20,0.20,0.20,0.15,0.25,0.20,-1.00,0.15])

        # colour range for giants
        g_r=np.array([2.55,4.20,4.90,5.30,1.55,2.45,2.70,2.80,3.70,3.90,-0.10,1.15])
        g_b=np.array([0.15,0.90,0.40,0.30,0.60,0.20,0.20,1.00,0.25,0.20,-1.40,0.15])

        clr0       = np.zeros([12,ntot])
        teff_cal   = np.zeros([12,ntot])

        clr0[0]  = bp-rp - (R_bp-R_rp)*ebv
        clr0[1]  = bp-j2 - (R_bp-R_j2)*ebv
        clr0[2]  = bp-h2 - (R_bp-R_h2)*ebv
        clr0[3]  = bp-k2 - (R_bp-R_k2)*ebv

        clr0[4]  = rp-j2 - (R_rp-R_j2)*ebv
        clr0[5]  = rp-h2 - (R_rp-R_h2)*ebv
        clr0[6]  = rp-k2 - (R_rp-R_k2)*ebv

        clr0[7]  = gp-j2 - (R_gg-R_j2)*ebv
        clr0[8]  = gp-h2 - (R_gg-R_h2)*ebv
        clr0[9]  = gp-k2 - (R_gg-R_k2)*ebv

        clr0[10] = gp-bp - (R_gg-R_bp)*ebv
        clr0[11] = gp-rp - (R_gg-R_rp)*ebv
        

        #derive Teff in all colour indices
        for j in range(0,12):
            teff_cal[j] = cpol[j,0] + cpol[j,1]*clr0[j] + cpol[j,2]*clr0[j]*clr0[j] + cpol[j,3]*clr0[j]*clr0[j]*clr0[j] + cpol[j,4]*clr0[j]*clr0[j]*clr0[j]*clr0[j]*clr0[j] + cpol[j,5]*logg + cpol[j,6]*logg*clr0[j] + cpol[j,7]*logg*clr0[j]*clr0[j] + cpol[j,8]*logg*clr0[j]*clr0[j]*clr0[j] + cpol[j,9]*logg*clr0[j]*clr0[j]*clr0[j]*clr0[j]*clr0[j] + cpol[j,10]*feh + cpol[j,11]*feh*clr0[j] + cpol[j,12]*feh*clr0[j]*clr0[j] + cpol[j,13]*feh*clr0[j]*clr0[j]*clr0[j] + cpol[j,14]*feh*logg*clr0[j]
            
        #go through each star (i) and each colour index (j)
        for i in range(0,ntot):
            if MC==True:
               dfeh    = feh [i]   + efeh[i] *np.random.normal(size=trials)
               dlog    = logg[i]   + elogg[i]*np.random.normal(size=trials)

               dcol[0]  = clr0[0,i]  + ebp[i]*np.random.normal(size=trials)-erp[i]*np.random.normal(size=trials) - (R_bp[i]-R_rp[i])*eebv[i]*np.random.normal(size=trials)
               dcol[1]  = clr0[1,i]  + ebp[i]*np.random.normal(size=trials)-ej2[i]*np.random.normal(size=trials) - (R_bp[i]-R_j2[i])*eebv[i]*np.random.normal(size=trials)
               dcol[2]  = clr0[2,i]  + ebp[i]*np.random.normal(size=trials)-eh2[i]*np.random.normal(size=trials) - (R_bp[i]-R_h2[i])*eebv[i]*np.random.normal(size=trials)
               dcol[3]  = clr0[3,i]  + ebp[i]*np.random.normal(size=trials)-ek2[i]*np.random.normal(size=trials) - (R_bp[i]-R_k2[i])*eebv[i]*np.random.normal(size=trials)

               dcol[4]  = clr0[4,i]  + erp[i]*np.random.normal(size=trials)-ej2[i]*np.random.normal(size=trials) - (R_rp[i]-R_j2[i])*eebv[i]*np.random.normal(size=trials)
               dcol[5]  = clr0[5,i]  + erp[i]*np.random.normal(size=trials)-eh2[i]*np.random.normal(size=trials) - (R_rp[i]-R_h2[i])*eebv[i]*np.random.normal(size=trials)
               dcol[6]  = clr0[6,i]  + erp[i]*np.random.normal(size=trials)-ek2[i]*np.random.normal(size=trials) - (R_rp[i]-R_k2[i])*eebv[i]*np.random.normal(size=trials)

               dcol[7]  = clr0[7,i]  + egg[i]*np.random.normal(size=trials)-ej2[i]*np.random.normal(size=trials) - (R_gg[i]-R_j2[i])*eebv[i]*np.random.normal(size=trials)
               dcol[8]  = clr0[8,i]  + egg[i]*np.random.normal(size=trials)-eh2[i]*np.random.normal(size=trials) - (R_gg[i]-R_h2[i])*eebv[i]*np.random.normal(size=trials)
               dcol[9]  = clr0[9,i]  + egg[i]*np.random.normal(size=trials)-ek2[i]*np.random.normal(size=trials) - (R_gg[i]-R_k2[i])*eebv[i]*np.random.normal(size=trials)

               dcol[10] = clr0[10,i] + egg[i]*np.random.normal(size=trials)-ebp[i]*np.random.normal(size=trials) - (R_gg[i]-R_bp[i])*eebv[i]*np.random.normal(size=trials)
               dcol[11] = clr0[11,i] + egg[i]*np.random.normal(size=trials)-erp[i]*np.random.normal(size=trials) - (R_gg[i]-R_rp[i])*eebv[i]*np.random.normal(size=trials)

            for j in range(0,12):
                dump=-1

                #remove Teffs outside of colour cuts
                if (logg[i]>3.2 and (clr0[j,i]>d_r[j] or clr0[j,i] < d_b[j])) \
                   or (np.isfinite(teff_cal[j,i])==False):
                    teff_cal[j,i]=float("nan")
                    dump=j
                if (logg[i]<=3.2 and (clr0[j,i]> g_r[j] or clr0[j,i] < g_b[j])) \
                  or (np.isfinite(teff_cal[j,i])==False):
                    teff_cal[j,i]=float("nan")
                    dump=j

                if MC==True:
                    if dump>=0: dteff[dump] = float("nan")
                    else: dteff[j] = cpol[j,0] + cpol[j,1]*dcol[j] + cpol[j,2]*dcol[j]*dcol[j] + cpol[j,3]*dcol[j]*dcol[j]*dcol[j] + cpol[j,4]*dcol[j]*dcol[j]*dcol[j]*dcol[j]*dcol[j] + cpol[j,5]*dlog + cpol[j,6]*dlog*dcol[j] + cpol[j,7]*dlog*dcol[j]*dcol[j] + cpol[j,8]*dlog*dcol[j]*dcol[j]*dcol[j] + cpol[j,9]*dlog*dcol[j]*dcol[j]*dcol[j]*dcol[j]*dcol[j] + cpol[j,10]*dfeh + cpol[j,11]*dfeh*dcol[j] + cpol[j,12]*dfeh*dcol[j]*dcol[j] + cpol[j,13]*dfeh*dcol[j]*dcol[j]*dcol[j] + cpol[j,14]*dfeh*dlog*dcol[j]
                    eteff[j]=np.std(dteff[j])
                    
            if MC==True:
                nte=np.where(np.isfinite(eteff)==True);nav=len(nte[0])
                if nav>=1:
                    fatto = 1
                    wei   = 1./(eteff[nte]*eteff[nte])
                    tca   = teff_cal[nte,i]

                    #weighted average from Teffs in all avaiable colour indices 
                    wteff = np.sum(tca*wei)/np.sum(wei)
                    #se   = 1./sqrt(np.sum(wei)) standard error of weighted mean

                    #weighted sample standard error is a better measure of 
                    #the overdispersion of the data wrt 1./sqrt(np.sum(wei))
                    #add 20K to account for zero-point uncertainty of Teff scale
                    wse   = np.sqrt(np.sum(wei*(tca-wteff)*(tca-wteff))/np.sum(wei)) + 20.
                    
                    if wato==False:
                        f.write('{},{:.2f},{:.2f},{:.3f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f}\n'.format(sid[i],logg[i],feh[i],ebv[i],\
                        teff_cal[0,i] ,eteff[0] +20.,teff_cal[1,i] ,eteff[1] +20.,\
                        teff_cal[2,i] ,eteff[2] +20.,teff_cal[3,i] ,eteff[3] +20.,\
                        teff_cal[4,i] ,eteff[4] +20.,teff_cal[5,i] ,eteff[5] +20.,\
                        teff_cal[6,i] ,eteff[6] +20.,teff_cal[7,i] ,eteff[7] +20.,\
                        teff_cal[8,i] ,eteff[8] +20.,teff_cal[9,i] ,eteff[9] +20.,\
                        teff_cal[10,i],eteff[10]+20.,teff_cal[11,i],eteff[11]+20.,\
                        wteff,wse))
                    else:
                        f.write('{},{:.2f},{:.2f},{:.3f},{:.0f},{:.0f}\n'.format(sid[i],logg[i],feh[i],ebv[i],wteff,wse))
            else:
                nte   = np.where(np.isfinite(teff_cal[:,i])==True);nav=len(nte[0])
                if nav>=1:
                    fatto=1
                    f.write('{},{:.2f},{:.2f},{:.3f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f},{:.0f}\n'.\
                        format(sid[i],logg[i],feh[i],ebv[i],\
                        teff_cal[0,i],teff_cal[1,i],teff_cal[2,i],teff_cal[3,i],\
                        teff_cal[4,i],teff_cal[5,i],teff_cal[6,i],teff_cal[7,i],\
                        teff_cal[8,i],teff_cal[9,i],teff_cal[10,i],teff_cal[11,i]))

        f.close()
        if fatto==0: print('*** No star passes quality requirements ***')
    else: print('*** No star satisfies basic input requirements ***')
