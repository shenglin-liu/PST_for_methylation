# PST_for_methylation
R functions for calculating PST between two groups based on methylation levels.

## Usage
1. Load the above two functions by sourcing this file in R (or simply copy and paste).

2. Convert your methylation data into a matrix in R like this one:
`     AA_cold AA_warm AR_cold AR_warm RA_cold RA_warm RR_cold RR_warm
[1,]   0.710   0.800   0.704   0.704   0.909   0.706   0.789   0.724
[2,]   0.727   0.710   0.586   0.792   0.750   0.586   0.815   0.846
[3,]   0.650   0.667   0.403   0.474   0.421   0.476   0.660   0.583
[4,]   0.852   0.827   0.830   0.917   0.922   0.902   0.830   0.805
[5,]   0.689   0.756   0.636   0.750   0.622   0.771   0.712   0.700
[6,]   0.867   0.700   0.864   0.846   0.783   0.909   0.893   0.852`

3. Suppose that you are comparing between the two temperature treatments. Run a code like this:
`group<-c(1,2,1,2,1,2,1,2)
pst_methyl<-PST_for_methylation(methyl,group)`

4. pst_methyl is the output. It is a data frame looking like this:
`          PST   meth1   meth2
1 0.261692638 0.77800 0.73350
2 0.017424287 0.71950 0.73350
3 0.018792170 0.53350 0.55000
4 0.007271477 0.85850 0.86275
5 0.820083908 0.66475 0.74425
6 0.109013212 0.85175 0.82675`

## NOTE
The sample size must be > 1 for each group.

## To cite
