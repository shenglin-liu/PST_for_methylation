# Functions for calculating between-group PST for methylation data.

#-----------Functions.
# Calculate PST for a quantitative trait among populations or groups.
# Y is a list of numeric vectors for a quantitative trait.
# Each vector is a population or a group.
# The element of each vector is an individual or a sample.
myPST<-function(Y)
{
	K<-length(Y)
	N<-sum(!is.na(unlist(Y)))
	ni<-sapply(Y,function(x)sum(!is.na(x)))
	mi<-sapply(Y,function(x)mean(x,na.rm=T))
	ma<-mean(unlist(Y),na.rm=T)
	vb<-sum((mi-ma)^2*ni)/(K-1)
	vw<-sum(sapply(Y,function(x)var(x,na.rm=T)*(sum(!is.na(x))-1)))/(N-K)
	vb/(vb+vw*2)
}
# Calculate PST of methylation between two groups.
# methyl: a numeric matrix of methylation levels;
#	rows correspond to methylation sites;
#	columns correspond to samples.
# group: an integer vector of 1's and 2's specifying the group attribution for the columns of methyl.
# Returns a 3-column data frame with the same number of rows as methyl. The columns represent:
#	PST;
#	mean methylation of Group 1;
#	mean methylation of Group 2.
PST_for_methylation<-function(methyl,group)
{
	meth1<-rowMeans(methyl[,group==1],na.rm=TRUE)
	meth2<-rowMeans(methyl[,group==2],na.rm=TRUE)
	PST<-sapply(1:nrow(methyl),function(x)myPST(split(methyl[x,],group)))
	data.frame(PST,meth1,meth2)
}

#-----------How to run.
# 1. Load the above two functions by sourcing this file in R (or simply copy and paste).

# 2. Convert your methylation data into a matrix in R like this one:
#     AA_cold AA_warm AR_cold AR_warm RA_cold RA_warm RR_cold RR_warm
#[1,]   0.710   0.800   0.704   0.704   0.909   0.706   0.789   0.724
#[2,]   0.727   0.710   0.586   0.792   0.750   0.586   0.815   0.846
#[3,]   0.650   0.667   0.403   0.474   0.421   0.476   0.660   0.583
#[4,]   0.852   0.827   0.830   0.917   0.922   0.902   0.830   0.805
#[5,]   0.689   0.756   0.636   0.750   0.622   0.771   0.712   0.700
#[6,]   0.867   0.700   0.864   0.846   0.783   0.909   0.893   0.852

# 3. Suppose that you are comparing between the two temperature treatments. Run a code like this:
#group<-c(1,2,1,2,1,2,1,2)
#pst_methyl<-PST_for_methylation(methyl,group)

# 4. pst_methyl is the output. It is a data frame looking like this:
#          PST   meth1   meth2
#1 0.261692638 0.77800 0.73350
#2 0.017424287 0.71950 0.73350
#3 0.018792170 0.53350 0.55000
#4 0.007271477 0.85850 0.86275
#5 0.820083908 0.66475 0.74425
#6 0.109013212 0.85175 0.82675

# NOTE: The sample size must be > 1 for each group.

