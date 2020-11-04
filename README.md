# Colour-T<sub>eff</sub> relations
-------------------------------------------------------
The IDL and python routines provided here are equivalent, and they allow users to derive photometric stellar effective temperatures (T<sub>eff</sub>) using the colour-T<sub>eff</sub> relations of Casagrande et al. (2020). These relations have been obtained implementing Gaia and 2MASS photometry in the InfraRed Flux Method, and applying it to some 300,000 stars across different evolutionary stages in the [GALAH DR3](https://galah-survey.org) survey. Please refer to the paper for further details on the method, and its validation against benchmark effective temperatures. Note that these relations are valid for Gaia DR2 photometry, and will be updated once Gaia EDR3 is out (a keyword will be included to choose between the two).

For each stars, up to 12 different T<sub>eff</sub> can be derived from various combinations of Gaia G, BP, RP and 2MASS J, H and K<sub>S</sub> photometry. Different colour indices have different range of applicability. Users have the choice to compute for each star as many T<sub>eff</sub> as possible, or to include a weighted averaged T<sub>eff</sub>. Weights are estimated with a MonteCarlo. Results are written in an output file, with the meaning of each column given in its header (e.g., "T_bprp" is T<sub>eff</sub> derived from the BP-RP colour index, and "eT_bprp" its MonteCarlo uncertainty. And so on for the other colour indices, "bpj", "rpk", "grp", etc ... "wTeff" is the weighted average for T<sub>eff</sub> and "ewTeff" its uncertainty).

See the routines for detailed examples of their use, required and optional input parameters. In short, to derive the effective temperature of a star, Gaia BP and RP photometry, surface gravity, metallicity, and reddening must be provided at the very least. If e.g., surface gravity is not know, one can assume log(g)=2 for a giant and log(g)=4 for a dwarf, and check the impact if different values were adopted. Same reasoning can be done for metallicity. If a star is unaffected by reddening, that should be passed as 0. The effects of surface gravity and metallicity will vary depending on colour index and T<sub>eff</sub> (stronger when moving to cooler stars). A star identifier, Gaia G and 2MASS J, H and K<sub>S</sub> photometry are also required inputs, but empty entries can be passed if some of these quantities are unavailable. The more photometric bands are provided, the more robustly a weighted averaged T<sub>eff</sub> can be derived. Optional input parameters are Gaia photometric flags ("phot_bp_rp_excess_factor" and "phot_proc_mode") and uncertainties on the adopted photometry, surface gravity, metallicity, and reddening. If uncertainties are not provided, the routines will assume reasonable uncertainties when doing a MonteCarlo. The routines will remove very bright stars for which Gaia and 2MASS photometry is likely to be unreliable. A few bare quality cuts on input photometric errors are performed, but users are ultimately responsible for what they pass in. 

Here are four different examples of calling these routines in python or IDL (note that the output file is always in csv format in the python routine, and ascii in the IDL one, unless the /csv keyword is set).

1. For each star, compute effective temperatures with MonteCarlo uncertainties based on some known input errors. T<sub>eff</sub> from each possible colour index and weighted averaged T<sub>eff</sub> are written into a filename called ``set1.dat``
```python
colte(sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,MC=True,ej2=ej2,eh2=eh2,ek2=ek2,eebv=ered,elogg=elogg,efeh=efeh,outfile='set1.dat')
```
```IDL
colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/MC,ej2=ej,eh2=eh,ek2=ek,eebv=ered,elogg=elogg,efeh=efeh,outfile='set1.dat'
```
-------------------------------------------------------
2. For each star, compute effective temperatures with MonteCarlo uncertainties based on default errors assumed by the routine. T<sub>eff</sub> from each possible colour index and weighted averaged T<sub>eff</sub> are written into a filename called ``set2.dat``
```python
colte(sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,MC=True,outfile='set2.dat')
```
```IDL
colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/MC,outfile='set2.dat'
```
-------------------------------------------------------
3. For each star, compute effective temperatures with MonteCarlo uncertainties based on default errors assumed by the routine. Only the weighted averaged T<sub>eff</sub> and its uncertainty are written into the default output file (``colte.csv`` in python, ``colte.dat`` in IDL)
```python
colte(sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,MC=True,wato=True')
```
```IDL
colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/MC,/wato
```
-------------------------------------------------------
4. For each star, compute effective temperatures in all available colour indices and dump the results into the default csv output file. Note that weighted averaged T<sub>eff</sub> is not computed since MC keyword is not set. 
```python
colte(sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv)
```
```IDL
colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/csv
```
