# Nunn-Qian-2014
Replication and critique of Nunn &amp; Qian (2014).  Joint Project with Guy Ishai.

Nunn & Qian (2014) find that wheat aid shipped from the US to developing countries increases conflicts in these countries on average. To deal with endogeneity, they use exogenous variation of US wheat production as an instrument. We find that using better controls for long-run time trends, the effect becomes statistically insignificant. In addition, using Monte Carlo simulations and arbitrary placebo data we find that a supposedly statistically significant effect on conflict patterns is not rare to find, implying that significance levels should be interpreted with caution.

*Final project in the course "Econometrics B" (Tel Aviv University)


This repo includes the following files:
1. final_gg_do.do - run this to reproduce our results in STATA
2. FAid_Final.dta - The original data file from Nunn & Qian (2014), downloaded from the article's AER page.
3. placebos.dta - A dataset with 3 placebo variables that we use
4. MC_Exports.dta - The dataset constructed by the Monte Carlo Simulation. Not needed in the do file but may be useful to browse the simulation results, as the simulation lasts around 12 hours.
5. FAid_GG.log - the log file of our final run
6. Readme.txt

To replicate results, make sure that files 1-3 are in the same directory, set the correct folder in the beginning of the do file, and run the do file.
