pro colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,DR2=DR2,DR3=DR3,bprp_ex=bprp_ex,$
          pmod=pmod,COD=COD,outfile=outfile,csv=csv,MC=MC,trials=trials,$
          wato=wato,elogg=elogg,efeh=efeh,egg=egg,ebp=ebp,erp=erp,ej2=ej2,$
          eh2=eh2,ek2=ek2,eebv=eebv
;+
; NAME:
;     COLTE
;
; PURPOSE:
;     Compute stellar effective temperatures using colour-Teff relations for the
;     Gaia and 2MASS photometric systems. User has to choose either Gaia DR2 or
;     DR3 photometry: no mixing of the two! The default extinction law is that
;     of Fitzpatrick (1999, renormalized as per Schlafly & Finkbeiner 2011 -
;     FSF). The option to use the extinction law of Cardelli, Clayton & Mathis
;     (1989, with optical from O'Donnell 1994 - COD) is available.
;
; EXPLANATION:
;     The relations used to derive Teff are from Casagrande et al. (2021).
;     For each star, Teffs are computed from up to 12 different colour indices
;     and results are written into a file.
;     If the option for a MonteCarlo is set, Teff uncertainties are computed
;     for each colour index, and a final weighted average Teff along with its
;     weighted standard deviation is derived. Teff from weighted average will 
;     likely have the best accuracy. However, in the pursue of precision, one 
;     might be better off by choosing Teff from colour indices with small
;     intrinsic scatter (see discussion in Section 4 of Casagrande+21).  
;     The routine applies a few bare quality cuts on input photometry by 
;     removing BP and RP<5, G<6, J<5.0, H<4.8, K<4.2, and if uncertainties are
;     passed in, also removing ej2>0.05, eh2>0.05, ek2>0.05. These cuts are 
;     mainly to avoid issues due to saturation at bright magnitudes (and to
;     some extent large photometric errors for faint 2MASS magnitudes).
;     Further quality cuts on Gaia photometry can be set with input parameters
;     bprp_ex= and pmod=    
;     Also, stars with ebv<0, logg < 0 or > 5, or feh > 0.6 will be excluded. 
;     Due to the decreased sensitivity of colours to low metallicities, stars 
;     with feh < -4 are assigned constant feh = - 4. Stars with feh < -8 are
;     assumed to not have a valid feh measurement, and are excluded.
;
; CALLING SEQUENCE:
; COLTE, SID, LOGG, FEH, GG, BP, RP, J2, H2, K2, EBV, /DR? [ BPRP_EX=, PMOD=, 
;        /COD, OUTFILE=, /CSV, /MC, TRIALS=, /WATO, ELOGG=, EFEH=, EGG=, EBP=,
;        ERP=, EJ2=, EH2=, EK2=, EEBV= ] 
;
; REQUIRED INPUT PARAMETERS
; sid:     star name/ID (strings and/or numbers can be passed) 
; logg:    surface gravity
; feh:     [Fe/H]
; gg:      Gaia G  (phot_g_mean_mag) 
; bp:      Gaia BP (phot_bp_mean_mag)
; rp:      Gaia RP (phot_rp_mean_mag)
; j2:      2MASS J
; h2:      2MASS H
; k2:      2MASS K
; ebv:     Reddening E(B-V)
; DR?:     Gaia DR2 or DR3 needs to be specified
;
; For each star, bp,rp,logg,feh,ebv are indispensable parameters needed 
; to derive Teff from at least bp-rp. Note that sid,gg,j2,h2,k2 are also 
; required inputs, but NaN entries can be passed if some of these quantities
; are unavaible for a star. It must also be specified whether photometry from
; Gaia DR2 or DR3 is passed as input.
;  
; CORRECTIONS TO INPUT PHOTOMETRY. DOS & DON'TS
; Gaia DR2: 6<G<16 are corrected by colte following Maiz Apellaniz & Weiler
;           (2018, A&A, 619, 180), with a constant zero-point offset for G>16
;           G<6 are excluded by colte to avoid any issue with saturation
; Gaia DR3: G magnitudes from DR3 (gaiadr3.gaia_source) should be passed as they
;           are. Note that G magnitudes from EDR3 (gaiaedr3.gaia_source) for
;           sources with 2 or 6-parameter astrometric solutions should be 
;           corrected as per Riello+21, A&A, 649, 3 (see
;           github.com/agabrown/gaiaedr3-6p-gband-correction) *before* passing
;           them to colte 
;           G<8 are corrected for saturation by colte following Riello+21, A&A,
;           649, 3 (Eq. C.1). G<6 are excluded by colte to avoid any issue with 
;           saturation            
;  
; OPTIONAL INPUT PARAMETERS
; bprp_ex: to remove stars with bad phot_bp_rp_excess_factor (For DR2 see Eq. 2,
;          Arenou+18, A&A, 616, 17. For DR3 see Eq. 2, Gaia Collaboration+21,
;          A&A, 649, 8). Note that for DR3 phot_bp_rp_excess_factor should be
;          used as given in the Gaia catalog, without applying the correction of
;          Riello+21, A&A, 649, 3, see:
;          github.com/agabrown/gaiaedr3-flux-excess-correction
;          If this correction is applied, user is in charge of changing the 
;          range of tolerance for the corrected excess factor (see suggested
;          values in Appendix A of Casagrande+21).           
;          If bprp_ex option is called, and excess is passed as NaN because
;          not available, the star will be removed
; pmod:    to retain only stars with phot_proc_mode=0  (Riello+18,+21)
;          If pmod option is called, and NaN is passed because not available,
;          the star will be removed
; COD:     to use extinction coefficients computed from the extinction law of
;          Cardelli, Clayton & Mathis (1989, with optical from O'Donnell 1994). 
;          If COD is not chosen, default extinction coefficients are from the
;          law of Fitzpatrick (1999, renormalized as per Schlafly & 
;          Finkbeiner 2011 - FSF)
; outfile: output file. If not passed, then the default output file is colte.dat
; csv:     to write output file in CSV format
; MC:      to perform a MonteCarlo for Teff uncertainties in different bands
; ej2:     2MASS J uncertainty. If not provided, 0.022 mag is assumed
; eh2:     2MASS H uncertainty. If not provided, 0.024 mag is assumed
; ek2:     2MASS K uncertainty. If not provided, 0.022 mag is assumed
;
; OPTIONAL INPUT PARAMETERS working ONLY if MC keyword is set
; trials:  number of MC realizations for each star. If not set, default is 1000
;          Default value is a good compromise between speed of execution and 
;          convergence. The latter depends on the colour index, and input
;          uncertainties. As a rule of thumb, with 1000 trials, uncertainties
;          typically converge to within a few K, or ~10K in worst cases. With
;          100 trials, convergence is ~10K in most cases, and up to ~70K in
;          worst cases. With 10000 trials convergence is always within a few K
; wato:    to write Weighted Averaged Teff Only in the output file 
; elogg:   logg uncertainty.    If not provided, 0.2 dex is assumed
; efeh:    [Fe/H] uncertainty.  If not provided, 0.1 dex is assumed
; egg:     Gaia G uncertainty.  If not provided, 0.005 mag is assumed
; ebp:     Gaia BP uncertainty. If not provided, 0.005 mag is assumed
; erp:     Gaia RP uncertainty. If not provided, 0.005 mag is assumed
; eebv:    Reddening uncertainty. If not provided, or a negative eebv is passed,
;          10% of input ebv is assumed
;
; OUTPUT 
; The routine will write an output file providing for each star the adopted
; sid, logg, feh, ebv + Teffs computed from up to 12 colour indices. If Teff
; cannot be determined in a colour index, NaN is returned for that index. 
; Note that the program makes a number of basic quality cuts on input
; photometry, and requires a value for logg, feh and ebv. Hence, the output
; file might contain fewer stars than the input file.
; If MC is set, then an uncertainty is provided for each Teff, along with
; weighted averaged Teff and weighted standard deviation. Note that all
; uncertainties are increased by 20K to account for the uncertainty on the
; zero-point of the Teff scale. However, when computing weighted averaged Teff, 
; weights do not factor this 20K increase (in case you wonder). If WATO is set,
; only weighted average and weighted standard deviation are written. Note that
; weighted averaged Teff and weighted standard deviation might change by a few
; Kelvin each time, because of the MC nature of the errors (more robust
; convergence can be achieved by increasing trials). 
; If CSV is set, output file is in csv format.
;
; EXAMPLES
; (1) For each star, compute Teffs with MonteCarlo uncertainties based on
;     known input errors. Colour-Teff relations for Gaia DR3 and default
;     extinction law (FSF) are used. Results from each colour index and weighted
;     average are written into filename set1.dat
; colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/DR3,/MC,ej2=ej,eh2=eh,ek2=ek,$
;       eebv=ered,elogg=elogg,efeh=efeh,outfile='set1.dat'
;
; (2) For each star, compute Teffs with MonteCarlo uncertainties based on
;     default errors assumed by the routine. Colour-Teff relations for Gaia DR2
;     and Cardelli/O'Donnel extinction law (COD) are used. Results from each
;     colour index and weighted average are written into filename set2.dat
; colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/DR2,/COD,/MC,outfile='set2.dat'
;
; (3) For each star, compute Teffs with MonteCarlo uncertainties based on
;     default errors assumed by the routine. Colour-Teff relations for Gaia DR3
;     and Cardelli/O'Donnel extinction law (COD) are used. Only weighted 
;     averaged Teff and its uncertainty are written into the default output
;     file colte.dat
; colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/DR3,/COD,/MC,/wato
;
; (4) For each star, compute Teffs in all available colour indices and dump
;     results into the default csv output colte.csv. Colour-Teff relations
;     for Gaia DR2 and default extinction law (FSF) are used.
; colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/DR2,/csv
; 
; HISTORY
; -November 2020 - Written by Luca Casagrande
; -July     2021 - Updated to include Gaia DR3 photometry and option to choose
;                  between COD and FSF extinction law
;- Sept     2022 - Weighted sample standard error if sample size > 1, otherwise
;                  standard deviation  
  
  ppm = intarr(n_elements(bp))
  if keyword_set(pmod) then ppm=pmod

  if not keyword_set(MC) and keyword_set(trials) then $
     print,'WARNING: trials is ignored since keyword MC is not set'
  
  if not keyword_set(MC) and keyword_set(wato)   then $
     print,'WARNING: keyword WATO is ignored since keyword MC is not set'
  
  if not keyword_set(MC) and keyword_set(elogg)  then $
     print,'WARNING: elogg is ignored since keyword MC is not set'
  
  if not keyword_set(MC) and keyword_set(efeh)  then $
     print,'WARNING: efeh is ignored since keyword MC is not set'
  
  if not keyword_set(MC) and keyword_set(egg)  then $
     print,'WARNING: egg is ignored since keyword MC is not set'

  if not keyword_set(MC) and keyword_set(ebp)  then $
     print,'WARNING: ebp is ignored since keyword MC is not set'
  
  if not keyword_set(MC) and keyword_set(erp)  then $
     print,'WARNING: erp is ignored since keyword MC is not set'
  
  if not keyword_set(MC) and keyword_set(eebv)  then $
     print,'WARNING: eebv is ignored since keyword MC is not set'
  
  if not keyword_set(MC) and keyword_set(ej2)  then $
     print,'WARNING: ej2 is used for quality cuts, but not to estimate Teff uncertainties since MC is not True'
  
  if not keyword_set(MC) and keyword_set(eh2)  then $
     print,'WARNING: eh2 is used for quality cuts, but not to estimate Teff uncertainties since MC is not True'
 
  if not keyword_set(MC) and keyword_set(ek2)  then $
     print,'WARNING: ek2 is used for quality cuts, but not to estimate Teff uncertainties since MC is not True'
  
  ; write output file with results or logs
  ; Note: any existing output file will be updated
  if keyword_set(outfile) then begin
     openw,lun,outfile,/get_lun
  endif else begin      
     if keyword_set(csv) then openw,lun,'colte.csv',/get_lun else $
        openw,lun,'colte.dat',/get_lun
  endelse

  ; Load colour-Teff polynomials, extinction coefficients and quality cuts for
  ; Gaia DR2 or DR3.
  ; First cut to retain only usable stars. Further cuts applied afterwards.
  ; -Gaia BP and RP needed for all stars. Finite function to remove NaN and Inf
  ;  in case there is any.  
  ; -Only BP and RP > 5 are considered because of uncalibrated systematics
  ;  at brighter magnitudes.
  ; -Stars without a value of reddening, logg or feh are excluded.
  ; -Stars with ebv<0, logg < 0 or > 5, feh < -8 or > 0.6 are also excluded.
  ; -If bprp_ex and/or pmod keywords are set, then stars with bad
  ;  phot_bp_rp_excess_factor and/or phot_proc_mode != 0 are also excluded.  

  data_2     = 0
  data_3     = 0
  data_track = 0
  ntot       = 0
  
  if keyword_set(DR2) then begin
     data_2 = 1

     ; colour-Teff coefficients 
     cpol   = fltarr(12,15)
     
     cpol[0,*]  = [7928.0505, -3663.1140,  803.3017,    -9.3727,        0., 325.1324, -500.1160,  279.4832,  -53.5062,      0.,   -2.4205, -128.0354, 49.4933,  5.9146,  41.3650] ; BP-RP
     cpol[1,*]  = [8217.8748, -2526.8430,  458.1827,   -28.4540,        0., 234.0113, -205.3084,   63.4781,   -7.2083,      0.,  -85.7048,  -50.1557, 32.3428, -2.3553,  20.0671] ; BP-J
     cpol[2,*]  = [8462.0737, -2570.3684,  537.5968,   -44.3644,        0., 189.1198, -106.7584,   31.1720,   -4.9137,      0.,   -9.2587, -189.8600, 75.8619, -6.8592,  16.7226] ; BP-H
     cpol[3,*]  = [8404.4760, -2265.1355,  403.4693,   -27.9056,        0., 193.5820, -145.3724,   47.7998,   -6.4572,      0.,  -34.5438, -130.2559, 52.6470, -4.4777,  15.8249] ; BP-K
     
     cpol[4,*]  = [9073.7917, -7670.6606, 3164.0525,         0., -126.1476,       0.,   -7.3816,  -12.5168,        0., -2.0452,        0.,   76.1144,      0.,      0., -45.8056] ; RP-J
     cpol[5,*]  = [8924.1553, -4779.3394, 1319.8989,         0.,  -16.6676,       0.,  -23.6583,   22.4243,        0., -4.3066,        0.,   35.0102,      0.,      0., -28.7228] ; RP-H
     cpol[6,*]  = [8940.4628, -4450.6138, 1138.6816,         0.,  -10.5749,       0.,  -42.3037,   33.3365,        0., -3.2535,        0.,   41.0402,      0.,      0., -21.9922] ; RP-K
     
     cpol[7,*]  = [8369.9905, -3559.7710,  895.8869,   -86.7011,        0., 180.7568, -164.9264,   24.4263,    4.2318,      0., -127.9640,   72.1449,      0.,      0.,  13.7683] ;  G-J
     cpol[8,*]  = [8185.8827, -2536.7671,  503.2762,   -42.7871,        0., 230.4871, -254.5291,  104.6258,  -17.4859,      0., -122.0732,   45.0572,      0.,      0.,   6.9992] ;  G-H
     cpol[9,*]  = [8103.2039, -1857.7194,        0.,    73.1834,   -1.7576, 236.0335, -345.9070,  170.4915,  -28.8549,      0., -131.4548,   49.6232,      0.,      0.,  10.0777] ;  G-K
     
     cpol[10,*] = [7555.3516,  5803.7715,        0., -2441.7124,  437.7314, 455.0997, 2243.1333, 3669.4924, 1872.7035,      0.,   19.1085,   75.2198,      0.,      0., -83.9777] ; G-BP 
     cpol[11,*] = [7971.3782, -5737.5049,        0.,  1619.9946, -203.8234, 255.7408, -492.8268,  160.1957,  103.1114,      0.,  -64.3289,   34.3339,      0.,      0.,  54.7224] ; G-RP
     
     ; Fitzpatrick/Schlafly extinction coefficients
     itbr = 0.8        
     cRg  = [2.608,-0.468, 0.048]
     cRb  = [3.007,-0.099,-0.212,0.069]
     cRr  = [1.702,-0.060]
     cRj  =  0.719
     cRh  =  0.455
     cRk  =  0.306
     
     if keyword_set(COD) then begin

        ; Cardelli/O'Donnel extinction coefficients
        itbr = 1.        
        cRg  = [3.068,-0.504, 0.053]        
        cRb  = [3.533,-0.114,-0.219,0.070]
        cRr  = [2.078,-0.073]
        cRj  =  0.899 
        cRh  =  0.567
        cRk  =  0.366
        
     endif
     
     if keyword_set(bprp_ex) then begin
     
        ok= where(finite(bp) eq 1 and finite(rp) eq 1 and bp ge 5.             $
                  and rp ge 5. and finite(ebv) eq 1 and finite(logg) eq 1      $
                  and finite(feh) eq 1 and ebv ge 0. and logg ge 0.            $
                  and logg le 5. and feh ge -8. and feh le 0.6 and             $
                  finite(bprp_ex) eq 1 and bprp_ex lt 1.3+0.06*(bp-rp)*(bp-rp) $
                  and bprp_ex gt 1.+0.015*(bp-rp)*(bp-rp) and finite(ppm) eq 1 $
                  and ppm eq 0,ntot)
     
     endif else begin
  
        ok= where(finite(bp) eq 1 and finite(rp) eq 1 and bp ge 5.        $
                  and rp ge 5. and finite(ebv) eq 1 and finite(logg) eq 1 $
                  and finite(feh) eq 1 and ebv ge 0. and logg ge 0.       $
                  and logg le 5. and feh ge -8. and feh le 0.6 and        $
                  finite(ppm) eq 1 and ppm eq 0,ntot)
     
     endelse

  endif  
  
  if keyword_set(DR3) then begin
     data_3 = 1

     ; colour-Teff coefficients 
     cpol   = fltarr(12,15)

     cpol[0,*]  = [7980.8845,  -4138.3457,  1264.9366,   -130.4388,         0.,   285.8393,   -324.2196,   106.8511,    -4.9825,        0.,     4.5138,  -203.7774, 126.6981, -14.7442,    40.7376] ; BP-RP
     cpol[1,*]  = [8172.2439,  -2508.6436,   442.6771,    -25.3120,         0.,   251.5862,   -240.7094,    86.0579,   -11.2705,        0.,   -45.9166,  -137.4645,  75.3191,  -8.7175,    21.5739] ; BP-J
     cpol[2,*]  = [8158.9380,  -2146.1221,   368.1630,    -24.4624,         0.,   231.8680,   -170.8788,    52.9164,    -6.8455,        0.,   -45.5554,  -142.9127,  55.2465,  -4.1694,    17.6593] ; BP-H
     cpol[3,*]  = [8265.6045,  -2124.5574,   355.5051,    -23.1719,         0.,   209.9927,   -161.4505,    50.5904,    -6.3337,        0.,   -27.2653,  -160.3595,  67.9016,  -6.5232,    16.5137] ; BP-K
     
     cpol[4,*]  = [9046.6493,  -7392.3789,  2841.5464,          0.,   -85.7060,         0.,    -88.8397,    80.2959,         0.,  -15.3872,         0.,    54.6816,       0.,       0.,   -32.9499] ; RP-J
     cpol[5,*]  = [8870.9090,  -4702.5469,  1282.3384,          0.,   -15.8164,         0.,    -30.1373,    27.9228,         0.,   -4.8012,         0.,    25.1870,       0.,       0.,   -22.3020] ; RP-H
     cpol[6,*]  = [8910.6966,  -4305.9927,  1051.8759,          0.,    -8.6045,         0.,    -76.7984,    55.5861,         0.,   -3.9681,         0.,    35.4718,       0.,       0.,   -16.4448] ; RP-K
     
     cpol[7,*]  = [8142.3539,  -3003.2988,   499.1325,     -4.8473,         0.,   244.5030,   -303.1783,   125.8628,   -18.2917,        0.,  -125.8444,    59.5183,       0.,       0.,    16.8172] ;  G-J
     cpol[8,*]  = [8133.8090,  -2573.4998,   554.7657,    -54.0710,         0.,   229.2455,   -206.8658,    68.6489,   -10.5528,        0.,  -124.5804,    41.9630,       0.,       0.,     7.9258] ;  G-H
     cpol[9,*]  = [8031.7804,  -1815.3523,         0.,     70.7201,    -1.7309,   252.9647,   -342.0817,   161.3031,   -26.7714,        0.,  -120.1133,    42.6723,       0.,       0.,    10.0433] ;  G-K
     
     cpol[10,*] = [7346.2000,   5810.6636,         0.,  -2880.3823,   669.3810,   415.3961,   2084.4883,  3509.2200,  1849.0223,        0.,   -49.0748,     6.8032,       0.,       0.,  -100.3419] ; G-BP 
     cpol[11,*] = [8027.1190,  -5796.4277,         0.,   1747.7036,  -308.7685,   248.1828,   -323.9569,  -120.2658,   225.9584,        0.,   -35.8856,   -16.5715,       0.,       0.,    48.5619] ; G-RP
     
     ; Fitzpatrick/Schlafly extinction coefficients
     itbr = 0.8        
     cRg  = [2.609,-0.475, 0.053]
     cRb  = [2.998,-0.140,-0.175,0.062]
     cRr  = [1.689,-0.059]
     cRj  =  0.719
     cRh  =  0.455
     cRk  =  0.306
     
     if keyword_set(COD) then begin

        ; Cardelli/O'Donnel extinction coefficients
        itbr = 1.        
        cRg  = [3.071,-0.511, 0.058]
        cRb  = [3.526,-0.168,-0.170,0.060]
        cRr  = [2.062,-0.072]
        cRj  =  0.899 
        cRh  =  0.567
        cRk  =  0.366
        
     endif
     
     if keyword_set(bprp_ex) then begin
     
        ok= where(finite(bp) eq 1 and finite(rp) eq 1 and bp ge 5.        $
                  and rp ge 5. and finite(ebv) eq 1 and finite(logg) eq 1 $
                  and finite(feh) eq 1 and ebv ge 0. and logg ge 0.       $
                  and logg le 5. and feh ge -8. and feh le 0.6            $
                  and finite(bprp_ex) eq 1 and bprp_ex gt 0.              $ 
                  and alog10(bprp_ex) lt 0.12 +0.039*(bp-rp)              $
                  and alog10(bprp_ex) gt 0.001+0.039*(bp-rp) and          $
                  finite(ppm) eq 1 and ppm eq 0,ntot)
     
     endif else begin
  
        ok= where(finite(bp) eq 1 and finite(rp) eq 1 and bp ge 5.        $
                  and rp ge 5. and finite(ebv) eq 1 and finite(logg) eq 1 $
                  and finite(feh) eq 1 and ebv ge 0. and logg ge 0.       $
                  and logg le 5. and feh ge -8. and feh le 0.6 and        $
                  finite(ppm) eq 1 and ppm eq 0,ntot)
     
     endelse

  endif
  
  if data_2 + data_3 eq 0 or data_2 + data_3 eq 2 then begin
     data_track = 1
     goto,iamsodone
  endif
  
  if ntot ge 1 then begin
     
     sid0   = string(sid(ok))
     gg0    = gg(ok)
     bp0    = bp(ok)
     rp0    = rp(ok)
     j20    = j2(ok)
     h20    = h2(ok)
     k20    = k2(ok)
     ebv0   = ebv(ok)
     logg0  = logg(ok)
     feh0   = feh(ok)

     fatto  = 0

     ; output file format
     outfor = '(A25,2(2X,F5.2),2X,F6.3,12(2X,F6.0))'
     if keyword_set(csv) then outfor = '(A25,%",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f")'
     
     ; things done if MC is set
     if keyword_set(MC) then begin
        
        outfor = '(A25,2(2X,F5.2),2X,F6.3,26(2X,F6.0))'
        
        if keyword_set(wato) then begin
           
           outfor = '(A25,2(2X,F5.2),2X,F6.3,2(2X,F6.0))'
           
           if keyword_set(csv) then begin
              printf,lun,'# star_name,logg,feh,ebv,wTeff,ewTeff'
              outfor = '(A25,%",%f,%f,%f,%f,%f")'
           endif else begin
              printf,lun,'#          star_name        logg    feh    ebv    wTeff  ewTeff  '
           endelse
              
        endif else begin
           
           if keyword_set(csv) then begin
              
              printf,lun,'# star_name,logg,feh,ebv,T_bprp,eT_bprp,T_bpj,eT_bpj,T_bph,eT_bph,T_bpk,eT_bpk,T_rpj,eT_rpj,T_rph,eT_rph,T_rpk,eT_rpk,T_gj,eT_gj,T_gh,eT_gh,T_gk,eT_gk,T_gbp,eT_gbp,T_grp,eT_grp,wTeff,ewTeff'
              outfor = '(A25,%",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f")'
              
           endif else begin
              
              printf,lun,'#          star_name        logg    feh    ebv   T_bprp  eT_bprp  T_bpj  eT_bpj   T_bph  eT_bph   T_bpk  eT_bpk   T_rpj  eT_rpj   T_rph  eT_rph   T_rpk  eT_rpk   T_gj   eT_gj    T_gh   eT_gh    T_gk   eT_gk    T_gbp  eT_gbp   T_grp  eT_grp   wTeff  ewTeff  '
              
           endelse
              
        endelse

        ; number of MC realizations
        if not keyword_set(trials) then trials = 1000

        dcol   = fltarr(12,trials) ; colour distribution
        dteff  = fltarr(12,trials) ; Teff distribution
        eteff  = fltarr(12)        ; Teff error for each colour index

       ; if errors are given then retain them. If errors are not given or
       ; set to zero, then reasonable default values are assigned. 
        if keyword_set(elogg) then begin
           nbad   = 0
           elogg0 = elogg(ok)
           bad    = where(elogg0 le 0. or finite(elogg0) eq 0,nbad)
           if nbad ge 1 then elogg0(bad) = 0.2
        endif else begin
           elogg0=fltarr(ntot)+0.2
        endelse
     
        if keyword_set(efeh) then begin
           nbad  = 0
           efeh0 = efeh(ok)
           bad   = where(efeh0 le 0. or finite(efeh0) eq 0,nbad)
           if nbad ge 1 then efeh0(bad) = 0.1
        endif else begin
           efeh0=fltarr(ntot)+0.1
        endelse

        if keyword_set(egg) then begin
           nbad = 0        
           egg0 = egg(ok)
           bad  = where(egg0 le 0. or finite(egg0) eq 0,nbad)
           if nbad ge 1 then egg0(bad) = 0.005
        endif else begin
           egg0=fltarr(ntot)+0.005
        endelse
        
        if keyword_set(ebp) then begin
           nbad = 0        
           ebp0 = ebp(ok)
           bad  = where(ebp0 le 0. or finite(ebp0) eq 0,nbad)
           if nbad ge 1 then ebp0(bad) = 0.005
        endif else begin
           ebp0=fltarr(ntot)+0.005
        endelse
        
        if keyword_set(erp) then begin
           nbad = 0        
           erp0 = erp(ok)
           bad  = where(erp0 le 0. or finite(erp0) eq 0,nbad)
           if nbad ge 1 then erp0(bad) = 0.005
        endif else begin
           erp0=fltarr(ntot)+0.005
        endelse
                
        ; default error for reddening is 10%
        if keyword_set(eebv) then begin
           nbad  = 0        
           eebv0 = eebv(ok)
           bad   = where(eebv0 le 0. or finite(eebv0) eq 0,nbad)
           if nbad ge 1 then eebv0(bad) = 0.1*ebv0(bad)
        endif else begin
           eebv0=0.1*ebv0
        endelse
        
        endif else begin
        
        if keyword_set(csv) then $
           printf,lun,'# star_name,logg,feh,ebv,T_bprp,T_bpj,T_bph,T_bpk,T_rpj,T_rph,T_rpk,T_gj,T_gh,T_gk,T_gbp,T_grp' else $
           printf,lun,'#          star_name        logg    feh    ebv   T_bprp   T_bpj   T_bph   T_bpk   T_rpj   T_rph   T_rpk   T_gj    T_gh    T_gk    T_gbp   T_grp'
        
     endelse

     ; as above, retain or assign resonable 2MASS errors
     if keyword_set(ej2) then begin
        nbad = 0        
        ej20 = ej2(ok)
        bad  = where(ej20 le 0. or finite(ej20) eq 0,nbad)
        if nbad ge 1 then ej20(bad) = 0.022
     endif else begin
        ej20=fltarr(ntot)+0.022
     endelse
     
     if keyword_set(eh2) then begin
        nbad = 0        
        eh20 = eh2(ok)
        bad  = where(eh20 le 0. or finite(eh20) eq 0,nbad)
        if nbad ge 1 then eh20(bad) = 0.024
     endif else begin
        eh20=fltarr(ntot)+0.024
     endelse
     
     if keyword_set(ek2) then begin
        nbad = 0        
        ek20 = ek2(ok)
        bad  = where(ek20 le 0. or finite(ek20) eq 0,nbad)
        if nbad ge 1 then ek20(bad) = 0.022
     endif else begin
        ek20=fltarr(ntot)+0.022
     endelse

     ; stars with -8 <= feh < -4 are assigned feh = -4.
     bad = where(feh0 lt -4.,nbad)
     if nbad ge 1 then feh0(bad) = -4.
     
     ; remove Gaia G<6 because of uncalibrated CCD saturation. Also remove Inf
     bad = where(gg0 lt 6. or finite(gg0) eq 0,nbad)
     if nbad ge 1 then gg0(bad) = !values.f_nan

     ; correct G from DR2 as per Maiz Apellaniz & Weiler (2018)
     if data_2 eq 1 then begin
        
        maw = where(gg0 lt 16.,nmaw)
        zpo = where(gg0 ge 16.,nzpo)
        
        if nmaw ge 1 then gg0(maw) = gg0(maw) - 0.0032*(gg0(maw)-6.)
        if nzpo ge 1 then gg0(zpo) = gg0(zpo) - 0.032
        
     endif

     ; if G<8 from DR3, then correct for saturation as per Riello et al. (2021)
     ; Note that G<6 are excluded anyway
     if data_3 eq 1 then begin
        
        sat = where(gg0 lt 8.,nsat)
        
        if nsat ge 1 then gg0(sat) = gg0(sat) - 0.09892 + 0.059*gg0(sat) - $
           0.009775*gg0(sat)*gg0(sat) + 0.0004934*gg0(sat)*gg0(sat)*gg0(sat)

     endif
     
     ; more house-cleaning to remove 2MASS photometry likely to be bad, or Inf
     noj = where(j20 lt 5.0 or ej20 gt 0.05 or finite(j20) eq 0,nnoj)
     noh = where(h20 lt 4.8 or eh20 gt 0.05 or finite(h20) eq 0,nnoh)
     nok = where(k20 lt 4.2 or ek20 gt 0.05 or finite(k20) eq 0,nnok)
     if nnoj ge 1 then j20(noj) = !values.f_nan
     if nnoh ge 1 then h20(noh) = !values.f_nan
     if nnok ge 1 then k20(nok) = !values.f_nan

     ; compute colour dependent extinction coefficients
     bprp0 = (bp0-rp0) - itbr*ebv0
     R_gg  = cRg[0] + cRg[1]*bprp0 + cRg[2]*bprp0*bprp0
     R_bp  = cRb[0] + cRb[1]*bprp0 + cRb[2]*bprp0*bprp0+cRb[3]*bprp0*bprp0*bprp0
     R_rp  = cRr[0] + cRr[1]*bprp0
     R_j2  = fltarr(ntot) + cRj
     R_h2  = fltarr(ntot) + cRh
     R_k2  = fltarr(ntot) + cRk

     ; colour range for dwarfs
     d_r=[2.00,3.00,4.00,4.20,1.05,1.60,1.85,2.10,2.60,2.80,-0.15,0.85]
     d_b=[0.20,0.25,0.40,0.30,0.20,0.20,0.20,0.15,0.25,0.20,-1.00,0.15]
     
     ; colour range for giants
     g_r=[2.55,4.20,4.90,5.30,1.55,2.45,2.70,2.80,3.70,3.90,-0.15,1.15]
     g_b=[0.20,0.90,0.40,0.30,0.60,0.20,0.20,1.00,0.25,0.20,-1.40,0.15]
     
     clr0       = fltarr(12,ntot)
     teff_cal   = fltarr(12,ntot)
     
     clr0[0,*]  = bp0-rp0 - (R_bp-R_rp)*ebv0        
     clr0[1,*]  = bp0-j20 - (R_bp-R_j2)*ebv0
     clr0[2,*]  = bp0-h20 - (R_bp-R_h2)*ebv0
     clr0[3,*]  = bp0-k20 - (R_bp-R_k2)*ebv0
     
     clr0[4,*]  = rp0-j20 - (R_rp-R_j2)*ebv0
     clr0[5,*]  = rp0-h20 - (R_rp-R_h2)*ebv0
     clr0[6,*]  = rp0-k20 - (R_rp-R_k2)*ebv0
     
     clr0[7,*]  = gg0-j20 - (R_gg-R_j2)*ebv0
     clr0[8,*]  = gg0-h20 - (R_gg-R_h2)*ebv0
     clr0[9,*]  = gg0-k20 - (R_gg-R_k2)*ebv0
     
     clr0[10,*] = gg0-bp0 - (R_gg-R_bp)*ebv0
     clr0[11,*] = gg0-rp0 - (R_gg-R_rp)*ebv0

     ; derive Teff in all colour indices
     for j=0,11 do begin

        teff_cal[j,*] = cpol[j,0] + cpol[j,1]*clr0[j,*] + cpol[j,2]*clr0[j,*]*clr0[j,*] + cpol[j,3]*clr0[j,*]*clr0[j,*]*clr0[j,*] + cpol[j,4]*clr0[j,*]*clr0[j,*]*clr0[j,*]*clr0[j,*]*clr0[j,*] + cpol[j,5]*logg0 + cpol[j,6]*logg0*clr0[j,*] + cpol[j,7]*logg0*clr0[j,*]*clr0[j,*] + cpol[j,8]*logg0*clr0[j,*]*clr0[j,*]*clr0[j,*] + cpol[j,9]*logg0*clr0[j,*]*clr0[j,*]*clr0[j,*]*clr0[j,*]*clr0[j,*] + cpol[j,10]*feh0 + cpol[j,11]*feh0*clr0[j,*] + cpol[j,12]*feh0*clr0[j,*]*clr0[j,*] + cpol[j,13]*feh0*clr0[j,*]*clr0[j,*]*clr0[j,*] + cpol[j,14]*feh0*logg0*clr0[j,*]

     endfor

     ; go through each star (i) and each colour index (j)
     for i=long(0),ntot-1 do begin
        
        if keyword_set(MC) then begin

           dfeh    = feh0 [i]  + efeh0[i] *randomn(seed,trials)
           dlog    = logg0[i]  + elogg0[i]*randomn(seed,trials)
           
           dcol[0,*]  = clr0[0,i]  + ebp0[i]*randomn(seed,trials)-erp0[i]*randomn(seed,trials) - (R_bp[i]-R_rp[i])*eebv0[i]*randomn(seed,trials)
           dcol[1,*]  = clr0[1,i]  + ebp0[i]*randomn(seed,trials)-ej20[i]*randomn(seed,trials) - (R_bp[i]-R_j2[i])*eebv0[i]*randomn(seed,trials)
           dcol[2,*]  = clr0[2,i]  + ebp0[i]*randomn(seed,trials)-eh20[i]*randomn(seed,trials) - (R_bp[i]-R_h2[i])*eebv0[i]*randomn(seed,trials)
           dcol[3,*]  = clr0[3,i]  + ebp0[i]*randomn(seed,trials)-ek20[i]*randomn(seed,trials) - (R_bp[i]-R_k2[i])*eebv0[i]*randomn(seed,trials)
           
           dcol[4,*]  = clr0[4,i]  + erp0[i]*randomn(seed,trials)-ej20[i]*randomn(seed,trials) - (R_rp[i]-R_j2[i])*eebv0[i]*randomn(seed,trials)
           dcol[5,*]  = clr0[5,i]  + erp0[i]*randomn(seed,trials)-eh20[i]*randomn(seed,trials) - (R_rp[i]-R_h2[i])*eebv0[i]*randomn(seed,trials)
           dcol[6,*]  = clr0[6,i]  + erp0[i]*randomn(seed,trials)-ek20[i]*randomn(seed,trials) - (R_rp[i]-R_k2[i])*eebv0[i]*randomn(seed,trials)
           
           dcol[7,*]  = clr0[7,i]  + egg0[i]*randomn(seed,trials)-ej20[i]*randomn(seed,trials) - (R_gg[i]-R_j2[i])*eebv0[i]*randomn(seed,trials)
           dcol[8,*]  = clr0[8,i]  + egg0[i]*randomn(seed,trials)-eh20[i]*randomn(seed,trials) - (R_gg[i]-R_h2[i])*eebv0[i]*randomn(seed,trials)
           dcol[9,*]  = clr0[9,i]  + egg0[i]*randomn(seed,trials)-ek20[i]*randomn(seed,trials) - (R_gg[i]-R_k2[i])*eebv0[i]*randomn(seed,trials)
           
           dcol[10,*] = clr0[10,i] + egg0[i]*randomn(seed,trials)-ebp0[i]*randomn(seed,trials) - (R_gg[i]-R_bp[i])*eebv0[i]*randomn(seed,trials)
           dcol[11,*] = clr0[11,i] + egg0[i]*randomn(seed,trials)-erp0[i]*randomn(seed,trials) - (R_gg[i]-R_rp[i])*eebv0[i]*randomn(seed,trials)
           
        endif        
        
        for j=0,11 do begin

           dump = -1
        
           ; remove Teffs outside of colour cuts
           if (logg0[i] gt 3.2 and (clr0[j,i] gt d_r[j] or clr0[j,i] lt d_b[j])) or $
              (finite(teff_cal[j,i]) eq 0) then begin
              teff_cal[j,i] = !values.f_nan
              dump = j
           endif
           
           if (logg0[i] le 3.2 and (clr0[j,i] gt g_r[j] or clr0[j,i] lt g_b[j])) or $
              (finite(teff_cal[j,i]) eq 0) then begin
              teff_cal[j,i] = !values.f_nan
              dump = j              
           endif
              
           if keyword_set(MC) then begin

              if dump ge 0 then begin
                 dteff[dump,*] = !values.f_nan
              endif else begin
                 dteff[j,*] = cpol[j,0] + cpol[j,1]*dcol[j,*] + cpol[j,2]*dcol[j,*]*dcol[j,*] + cpol[j,3]*dcol[j,*]*dcol[j,*]*dcol[j,*] + cpol[j,4]*dcol[j,*]*dcol[j,*]*dcol[j,*]*dcol[j,*]*dcol[j,*] + cpol[j,5]*dlog + cpol[j,6]*dlog*dcol[j,*] + cpol[j,7]*dlog*dcol[j,*]*dcol[j,*] + cpol[j,8]*dlog*dcol[j,*]*dcol[j,*]*dcol[j,*] + cpol[j,9]*dlog*dcol[j,*]*dcol[j,*]*dcol[j,*]*dcol[j,*]*dcol[j,*] + cpol[j,10]*dfeh + cpol[j,11]*dfeh*dcol[j,*] + cpol[j,12]*dfeh*dcol[j,*]*dcol[j,*] + cpol[j,13]*dfeh*dcol[j,*]*dcol[j,*]*dcol[j,*] + cpol[j,14]*dfeh*dlog*dcol[j,*]
              endelse

              eteff[j]=stddev(dteff[j,*])
              
           endif
           
        endfor
         
        if keyword_set(MC) then begin

           nte   = where(finite(eteff) eq 1,nav)
           
           if nav ge 1 then begin
              
              fatto = 1 
              wei   = 1./(eteff(nte)*eteff(nte))
              tca   = teff_cal[nte,i]

              ; weighted average from Teffs in all available colour indices
              wteff = total(tca*wei)/total(wei)
              ; se  = 1./SQRT(total(wei)) standard error of the weighted mean

              ; If sample size >1, weighted sample standard error is a better  
              ; measure of the overdispersion of the data wrt 1/SQRT(total(wei))
              ; Add 20K to account for zero-point uncertainty of the Teff scale
              if nav eq 1 then wse = eteff(nte)+20. else $
                 wse = SQRT(total(wei*(tca-wteff)*(tca-wteff))/total(wei)) + 20.

              if not keyword_set(wato) then $
              printf,lun,sid0[i],logg0[i],feh0[i],ebv0[i],$
                     teff_cal[0,i] ,eteff[0] +20.,teff_cal[1,i] ,eteff[1] +20.,$
                     teff_cal[2,i] ,eteff[2] +20.,teff_cal[3,i] ,eteff[3] +20.,$
                     teff_cal[4,i] ,eteff[4] +20.,teff_cal[5,i] ,eteff[5] +20.,$
                     teff_cal[6,i] ,eteff[6] +20.,teff_cal[7,i] ,eteff[7] +20.,$
                     teff_cal[8,i] ,eteff[8] +20.,teff_cal[9,i] ,eteff[9] +20.,$
                     teff_cal[10,i],eteff[10]+20.,teff_cal[11,i],eteff[11]+20.,$
                     wteff,wse,$
                     format=outfor else $
              printf,lun,sid0[i],logg0[i],feh0[i],ebv0[i],wteff,wse,$
                     format=outfor
              
           endif           
                      
        endif else begin

           nte   = where(finite(teff_cal[*,i]) eq 1,nav)
           
           if nav ge 1 then begin

              fatto = 1
                            
              printf,lun,sid0[i],logg0[i],feh0[i],ebv0[i],teff_cal[0,i],teff_cal[1,i],teff_cal[2,i],teff_cal[3,i],teff_cal[4,i],teff_cal[5,i],teff_cal[6,i],teff_cal[7,i],teff_cal[8,i],teff_cal[9,i],teff_cal[10,i],teff_cal[11,i],format=outfor 
           endif
              
        endelse

     endfor     
          
     if fatto eq 0 then begin

        print     ,'*** No star passes quality requirements ***'
        printf,lun,'*** No star passes quality requirements ***'

     endif
        
     endif else begin

        print     ,'*** No star satisfies basic input requirements to derive its Teff ***'
        printf,lun,'*** No star satisfies basic input requirements to derive its Teff ***'

     endelse
  
     iamsodone: if data_track eq 1 then begin

        print,'*** You must choose colour-Teff relations from either Gaia DR2 or DR3 ***'
        printf,lun,'*** You must choose colour-Teff relations from either Gaia DR2 or DR3 ***'
        
     endif
     
     free_lun,lun
     close,lun

     ; Avoid 'floating illegal operand' message in some OS in case
     ; relational operations involve NaN
     error = check_math(MASK=128)
     
  end
