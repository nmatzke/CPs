#######################################################
# Example Ancestral Character Estimation ("ace")
# on Carnivorous Plant Trap traits
#######################################################

library(ape)
library(phytools)


wd = "/GitHub/CPs/examples/all_v1"
setwd(wd)


trfn = "gbotb_tr14_wSister_genera_edited.newick"
states_fn = "gbotb_tr14_wSister_genera_edited_states.txt"

tr = read.tree(trfn)
length(tr$tip.label)

ddf = read.table(states_fn, header=TRUE, sep="\t")
head(ddf)
dim(ddf)

# Ensure appropriate sorting
states = ddf$state
states = c(as.numeric(states) + 1)
names(states) = ddf$species
head(states)

summary(as.factor(states))
#    1    2    3    4    5    6    7    8    9   10   11 
# 2215   17  128    2    1    6   15  137   46    1   80 
# Ancestral Character Estimation ("ace")
# ER = equal rates


# See "rates matrix" Excel spreadsheet for how this works!
Matzke_2002_model = matrix(data=c(0,2,2,0,0,0,0,2,0,0,0,
1,0,3,0,0,6,0,0,0,0,0,
1,3,0,4,0,0,6,0,0,0,0,
1,0,0,0,5,0,0,0,0,0,0,
1,0,0,0,0,0,0,0,0,0,0,
1,6,0,0,0,0,3,6,0,0,0,
1,0,6,0,0,3,0,0,6,0,0,
1,0,0,0,0,6,0,0,3,0,0,
1,0,0,0,0,0,6,3,0,4,0,
1,0,0,0,0,0,0,0,0,0,5,
1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

Matzke_2002_model


equal_rates_model = matrix(data=1, nrow=11, ncol=11, byrow=TRUE)
diag(equal_rates_model) = 0
equal_rates_model


resER = fitMk.parallel(tree=tr, x=states, model="ER", fixedQ=NULL); save(resER, file="resER_noParallel.Rdata")

resER
# Log-likelihood: -581.909837

resER = fitMk.parallel(tree=tr, x=states, model="ER", fixedQ=NULL, ncores=11)
save(resER, file="resER.Rdata")

resM02 = fitMk.parallel(tree=tr, x=states, model=Matzke_2002_model, fixedQ=NULL, ncores=11)
save(resER, file="resM02.Rdata")
