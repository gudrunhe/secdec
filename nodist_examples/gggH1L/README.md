# PSD_gggH1L
Numerical evaluation of the amplitude for gg->gH process at 1-Loop LO using pySecDec

The expressions are based on the "HiggsTools Handbook", https://arxiv.org/abs/1711.09875, Section II.5.2 and II.5.3.

gggH1L:  contains all the files necessary for the numerical evaluation of the required amplitude with pySecDec

# Calculation & Verification

In the gggH1L folder, execute the following commands to numerically evaluate the required normalized amplitude and to also verify the numerical results with the analytic results in the HTL-Heavy Top Limit

$ python3 generate_F212.py

$ make -C F212/

$ python3 generate_F312.py

$ make -C F312/

$ python3 integrate_gggH1L.py

# Sample Terminal Output

----------------------------------------------------------------------------------------------------------------------

Numerical Results for all normalized Form Factors up to a phase factor :

Normalized F212 :  2228.05719411905375 - 8.28197625818457052e-7* I  +/-  1.87131842528517779e-6 + 2.94978429315952658e-6* I

Normalized F311 :  -1481.48161841203819 + 1.00900601941613428e-6* I  +/-  5.03777341256770348e-6 + 5.44341369208350936e-6* I

Normalized F332 :  4444.44484616036971 - 1.8636243843809954e-6* I  +/-  3.8668284992959365e-6 + 5.56577270100101655e-6* I 

Normalized F312 :  5191.0537518748215 + 0.0000153041235741549214* I  +/-  3.94613453985244828e-6 + 6.20950451773513458e-6* I

----------------------------------------------------------------------------------------------------------------------

Analytic results of all normalized Form Factors with heavy top quark limit

Analytic Result Normalized F212 : -2228.0569919594127

Analytic Result Normalized F311 : 1481.4814814814813

Analytic Result Normalized F332 : -4444.444444444444

Analytic Result Normalized F312 : -5191.019954922376

----------------------------------------------------------------------------------------------------------------------

Numerical result of the helicity amplitudes up to a phase factor, with the above Form Factors : 

----------------------------------------------------------------------------------------------------------------------

M++-  :  0.0300393616026186366 +I*( -1.11660185501779895e-11 )
 
 M+-+  :  -0.00333770685039121934 +I*( 2.27324204447506315e-12 )
 
 M-++  :  0.0132809934301902414 +I*( -5.56892571784054173e-12 )
 
 M+++  :  2.07980125993722556e-7 +I*( 1.56564829313294109e-10 )
                                                                                                                                                                                             
 The Other 4 Helicity Amplitudes are Identical to the above ones, upto a phase factor 

----------------------------------------------------------------------------------------------------------------------

Numerical result of the normalized amplitude squared with the above Form Factors : 
                                                                                                                                                                                             
| M_(g g -> g H) |^2 :  0.00217977663809607881  : Here we have used Eq. (II.48a)  of  1711.09875
                                                                                                                                                                                             
| M_(g g -> g H) |^2 :  0.00217977663809607898  : Here we have added all the mod squared helicity amplitudes 
                                                                                                                                                                                             
Multiply above number with ( ( Nc*(Nc^2 -1) )* (topmass^4)*(alpha_S^3) / pi*(v^2) )  To get the final ampltiude squared, which can then be averaged over colors and gluon polarisations 

----------------------------------------------------------------------------------------------------------------------

Analytic result of the normalized amplitude squared in the heavy top quark limit
                                                                                                                                                                                             
| M_(g g -> g H) |^2 :  0.0021797762426380867
                                                                                                                                                                                             
Multiply above number with ( ( Nc*(Nc^2 -1) )* (topmass^4)*(alpha_S^3) / pi*(v^2) )  to get the final ampltiude squared, which can then be averaged over colors and gluon polarisations 

----------------------------------------------------------------------------------------------------------------------

Both numerical and analytic amplitudes have been evaluated with the following kinematics : 

----------------------------------------------------------------------------------------------------------------------

 [ s12, s13, s23, topmass^2, hmass^2 ]  = [  0.0009 ,  -0.0003 ,  -0.00059842873775 ,  1.0 ,  1.5712622500000002e-06  ]

----------------------------------------------------------------------------------------------------------------------

 topmass =  1.0   |  hmass =  0.0012535  |  s1/(topmass^2) =  0.0009  |  s2/(topmass^2) =  -0.0003  |  s3/(topmass^2) =  -0.00059842873775

----------------------------------------------------------------------------------------------------------------------
