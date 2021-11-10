# Gaia DR2, DR3 and 2MASS colour-T<sub>eff</sub> relations
-------------------------------------------------------
The IDL and python routines provided here are equivalent, and they allow users to derive photometric stellar effective temperatures (T<sub>eff</sub>) using the colour-T<sub>eff</sub> relations of [Casagrande et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2684C/abstract). These relations have been obtained implementing Gaia and 2MASS photometry in the InfraRed Flux Method, and applying it to over 360,000 stars across different evolutionary stages in the [GALAH+ DR3](https://docs.datacentral.org.au/galah/dr3/overview/) survey. Please refer to the [paper](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2684C/abstract) for further details on the method, and its validation against benchmark effective temperatures. Users must specify whether they are using Gaia DR2 or EDR3 photometry.
 - Gaia DR2 photometry should be _passed as is_. For G>6 the routine will correct the [Gaia DR2 G magnitude dependent offset](https://ui.adsabs.harvard.edu/abs/2018MNRAS.479L.102C/abstract) according to [this paper](https://ui.adsabs.harvard.edu/abs/2018A%26A...619A.180M/abstract). G<6 will be removed to avoid uncalibrated systematics. 
 - Gaia EDR3 G-band photometry for sources with 2 or 6-parameter astrometric solutions is **not** corrected by the routine. It is responsibility of the user to do so as per [here](https://github.com/agabrown/gaiaedr3-6p-gband-correction). _No further correction to input EDR3 data is needed by the routine_. The routine will [corrected for saturation](https://ui.adsabs.harvard.edu/abs/2021A%26A...649A...3R/abstract) G<8. In abundance of caution, note that G<6 will still be removed to avoid uncalibrated systematics. If the option to apply a quality cut based on "phot_bp_rp_excess_factor" is used, this excess factor must be passed as given in the Gaia catalog (i.e., do **not** use the corrected factor defined [here](https://github.com/agabrown/gaiaedr3-flux-excess-correction)).

For each stars, up to 12 different T<sub>eff</sub> can be derived from various combinations of Gaia G, BP, RP and 2MASS J, H and K<sub>S</sub> photometry. Different colour indices have different range of applicability. Users have the choice to compute for each star as many T<sub>eff</sub> as possible. If a MonteCarlo is called, a weighted averaged T<sub>eff</sub> is also provided. Effective temperatures from weighted average will likely have the best accuracy. However, in the pursue of precision, one might be better off by choosing effective temperatures from colour indices with small intrinsic scatter. Results are written in an output file, with the meaning of each column given in its header (e.g., "T_bprp" is T<sub>eff</sub> derived from the BP-RP colour index, and "eT_bprp" its MonteCarlo uncertainty, and so on for the other colour indices, "bpj", "rpk", "grp", etc... . "wTeff" is the weighted average for T<sub>eff</sub> and "ewTeff" its uncertainty). 

See the routine for a detailed discussion on the required and optional input parameters, and practical examples. In short, to derive the effective temperature of a star, Gaia BP and RP photometry, surface gravity, metallicity, and reddening must be provided at the very least. It must also be specified whether Gaia DR2 or EDR3 photometry is used. The routine will use colour dependent extinction coefficients derived from an assumed [extinction law](https://ui.adsabs.harvard.edu/abs/1999PASP..111...63F/abstract) to deredden magnitudes. The option to use extinction coefficients from a [different law](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C/abstract) is available. If e.g., surface gravity is not know, one can assume log(g)=2 for a giant and log(g)=4 for a dwarf, and check the impact if different values were adopted (using the MonteCarlo option and allowing for large uncertainties in log(g) is also advisable, although uncertainties will be centred around assumed input values). Same reasoning can be done for metallicity. If a star is unaffected by reddening, that should be passed as 0. The effects of surface gravity and metallicity will vary depending on colour index and T<sub>eff</sub> (stronger when moving to cooler stars). A star identifier, Gaia G and 2MASS J, H and K<sub>S</sub> photometry are also required inputs, but empty entries can be passed if some of these quantities are unavailable. The more photometric bands are provided, the more robustly a weighted averaged T<sub>eff</sub> can be derived. Optional input parameters are Gaia photometric flags ("phot_bp_rp_excess_factor" and "phot_proc_mode") and uncertainties on the adopted photometry, surface gravity, metallicity, and reddening. If uncertainties are not provided, the routines will assume reasonable uncertainties when doing a MonteCarlo. The routines will remove very bright stars for which Gaia and 2MASS photometry is likely to be unreliable, but users are ultimately responsible for what they pass in. 

Here are four different examples of calling these routines in python or IDL. Note that the output file is always in csv format for the python routine, and ascii for IDL (unless the /csv keyword is set).

1. For each star, compute effective temperatures with MonteCarlo uncertainties based on some known input errors. Colour-T<sub>eff</sub> relations for Gaia EDR3 and default extinction law are used. Effective temperatures from each possible colour index and weighted averaged T<sub>eff</sub> are written into a filename called ``set1.dat``
```python
colte(sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,DR3=True,MC=True,ej2=ej2,eh2=eh2,ek2=ek2,eebv=ered,elogg=elogg,efeh=efeh,outfile='set1.dat')
```
```IDL
colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/DR3,/MC,ej2=ej,eh2=eh,ek2=ek,eebv=ered,elogg=elogg,efeh=efeh,outfile='set1.dat'
```
-------------------------------------------------------
2. For each star, compute effective temperatures with MonteCarlo uncertainties based on default errors assumed by the routine. Colour-T<sub>eff</sub> relations for Gaia DR2 and Cardelli/O'Donnel (COD) extinction law are used. Effective temperatures from each possible colour index and weighted averaged T<sub>eff</sub> are written into a filename called ``set2.dat``
```python
colte(sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,DR2=True,COD=True,MC=True,outfile='set2.dat')
```
```IDL
colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/DR2,/COD,/MC,outfile='set2.dat'
```
-------------------------------------------------------
3. For each star, compute effective temperatures with MonteCarlo uncertainties based on default errors assumed by the routine. Colour-T<sub>eff</sub> relations for Gaia EDR3 and Cardelli/O'Donnel (COD) extinction law are used. Only the weighted averaged T<sub>eff</sub> and its uncertainty are written into the default output file (``colte.csv`` in python, ``colte.dat`` in IDL)
```python
colte(sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,DR3=True,COD=True,MC=True,wato=True')
```
```IDL
colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/DR3,/COD,/MC,/wato
```
-------------------------------------------------------
4. For each star, compute effective temperatures in all available colour indices and dump the results into the default csv output file. Colour-T<sub>eff</sub> relations for Gaia DR2 and default extinction law are used. Note that weighted averaged T<sub>eff</sub> is not computed since MC keyword is not set
```python
colte(sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,DR2=True)
```
```IDL
colte,sid,logg,feh,gg,bp,rp,j2,h2,k2,ebv,/DR2,/csv
```
