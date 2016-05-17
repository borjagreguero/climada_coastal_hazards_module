{\rtf1\ansi\ansicpg1252\cocoartf1265\cocoasubrtf210
{\fonttbl\f0\fmodern\fcharset0 Courier;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720

\f0\fs24 \cf0 17-March-2014\
\
Here are the files which contain the ocean and ice components, sums and uncertainties as used in preparing the IPCC AR5 report (2014), with some slight modifications.  One small choice I made here was to combine the ocean and inverse barometer effect into one field, both for the mean and the uncertainty.  I also have provided better smoothed maps for the *time series* (the 20-mean-difference maps are the same as in the report).  This actually shouldn't be the cause for any difference in the report figures, as I didn't use this field for anything but the coastal stations in Fig. 13.23, and they have the same interpolation scheme at the coast now, just a better interpolation scheme in the open ocean (bilinear; not shown in any figure in the report).\
\
One thing to note: I made a choice as to how to provide the 5% and 95% (upper and lower 90% confidence interval) uncertainty estimates for the various fields.  I have provided the maps of these similar to the way Jonathan Gregory provided the time series to me, as the individual component upper and lower bounds. However, to combine these errors in the same way as in the report, you will need to take the difference between the upper bound and the middle value (for combining the upper uncertainty total estimate) or the lower bound and middle value (for combining the lower uncertainty total estimate), and use the formula shown in the Supplementary Material for Chapter 13 of the AR5 (SM13) - minus the IBE which is combined with the ocean field here.  Like this for the time series components datasets (for python code):\
\
\
import numpy as np\
# for the total sea level rise signal, plus gia\
totslr = ocn_m + antdyn_m + greendyn_m + glac_m + antsmb_m + greensmb_m + grw_m + gia_ts[:-1]\
\
# lower uncertainty bound, \
# adding the uncertainty back to the total so as to plot the lower total bound\
loerr = totslr - np.sqrt((np.abs(ocn_l-ocn_m) + np.abs(antsmb_l-antsmb_m) + np.abs(greensmb_l-greensmb_m))**2 + (glac_m-glac_l)**2 + (1.645*giasd_ts[:-1])**2 + (grw_m-grw_l)**2 + (greendyn_m-greendyn_l)**2 + (antdyn_m-antdyn_l)**2)\
\
# upper uncertainty bound, \
# adding the uncertainty back to the total so as to plot the upper total bound\
hierr =  totslr + np.sqrt((np.abs(ocn_h-ocn_m) + np.abs(antsmb_h-antsmb_m) + np.abs(greensmb_h-greensmb_m))**2 + (glac_m-glac_h)**2 + (1.645*giasd_ts[:-1])**2 + (grw_m-grw_h)**2 + (greendyn_m-greendyn_h)**2 + (antdyn_m-antdyn_h)**2)\
\
-----\
The GIA fields are also provided here:\
\
gia_peltier.nc - Peltier et al. 2004 GIA field\
\
gia_lambeck.nc - GIA field based on Lambeck ice history, using modified ANU model (see SM13 for details)\
\
gia_mean.nc - the mean of these two GIA fields; this is the GIA field used in the AR5 report\
\
gia_diff.nc - the absolute half-difference between these two GIA fields (i.e., the absolute difference between the GIA mean field and one of the two input GIA fields); this is the uncertainty used for GIA in the AR5 report\
\
\
The GIA time series, for the calculations shown above, were calculated as follows (again, python code):\
\
a1 = arange(105.).reshape(-1,1,1)\
\
dgia = nc4.Dataset('gia_mean.nc')\
gia = dgia.variables['rsl'][:]\
gia1 = (gia/95.)*a1\
gia_ts = gia1[10:]\
\
dgiasd = nc4.Dataset('gia_diff.nc')\
giasd = dgiasd.variables['giasd'][:]\
giasd1 = (giasd / 95.) * a1\
giasd_ts = giasd1[10:]\
\
-----\
\
Let me know if you have any questions, but please, after reading through the notes in the files and SM13.\
\
Thank you,\
\
Mark Carson\
mark.carson@zmaw.de\
\
\
\
References\
\
Chapter 13 paper:\
Church, J. A., P. Clark, A. Cazenave, J. Gregory, S. Jevrejeva, A. Levermann, M. Merrifield, G. Milne, R.S.Nerem, P. Nunn, A. Payne, W. Pfeffer, D. Stammer, and A. Unnikrishnan (2013), Sea level change, in Climate Change 2013: The Physical Science Basis, edited by T. F. Stocker, D. Qin, G.-K. Plattner, M. Tignor, S. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex, and P. Midgley, Cambridge University Press, Cambridge, UK and New York, NY. USA.\
\
Available here as of 18Mar2014: http://www.climatechange2013.org/images/report/WG1AR5_Chapter13_FINAL.pdf\
\
Supplementary Material:\
Church, J.A., P.U. Clark, A. Cazenave, J.M. Gregory, S. Jevrejeva, A. Levermann, M.A. Merrifield, G.A. Milne, R.S.\
Nerem, P.D. Nunn, A.J. Payne, W.T. Pfeffer, D. Stammer and A.S. Unnikrishnan, 2013: Sea Level Change Supplementary Material. In:\
Climate Change 2013: The Physical Science Basis. Contribution of Working Group I to the Fifth Assessment Report of the Intergovernmental Panel on Climate Change [Stocker, T.F., D. Qin, G.-K. Plattner, M. Tignor, S.K. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex and P.M. Midgley (eds.)]. Available from www.climatechange2013.org and www.ipcc.ch.\
\
(Specifically available from here, as of 18Mar2014)\
http://www.climatechange2013.org/images/report/WG1AR5_Ch13SM_FINAL.pdf\
\
\
}