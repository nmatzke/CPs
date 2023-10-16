
library(ape)
library(phytools)
library(BioGeoBEARS)
library(openxlsx)

wd = "/GitHub/CPs/examples/Droseraceae_v1/"
setwd(wd)

trfn = "Sen_2020_Fig1_Droceraceae_digitized_v1.newick"

tr = read.tree(trfn)
plot(tr, cex=0.7)
axisPhylo()
title("Droseraceae phylogeny from Sen et al (2020), Fig. 1")

# Print tip labels to screen, for pasting into Excel
cat(tr$tip.label, sep="\n")


# Read in Excel file
xls_fn = "Droseraceae_simple_chars_v1.xlsx"
xlsx = openxlsx::read.xlsx(xls_fn)

# Format tip data for phytools
head(xlsx)
tipdata = xlsx$char
names(tipdata) = xlsx$names
state_names = c("pitcher_air", "pitcher_ground", "sticky_air", "sticky_ground", "snap_amphibious", "snap_aquatic")

# ER = equal rates
# pi = root state frequencies (default equal)
resER = phytools::fitMk(tree=tr, x=tipdata, model="ER", pi="equal")

# SYM = symmetric rates
resSYM = phytools::fitMk(tree=tr, x=tipdata, model="SYM", pi="equal")

# ARD = All Rates Different (will be slow)
# resARD = phytools::fitMk(tree=tr, x=tipdata, model="SYM", pi="equal")

# compare the models
AIC(resER, resSYM)




# SIMMAP: stochastic mapping of possible histories
resER_smap = simmap(resER)

# print a summary of the stochastic mapping
summary(resER_smap)


# plot a posterior probabilities of ancestral state

pdffn = "resER_smap.pdf"
pdf(file=pdffn, width=12, height=12)

cols = setNames(c("violet", "blue","green3", "yellow", "orange2", "red"),sort(unique(tipdata)))
cols
plot(summary(resER_smap), colors=cols, ftype="i")

legend("bottomleft", state_names, pch=21, pt.bg=cols, pt.cex=1)


dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)



pdffn = "resER_smap_density.pdf"
pdf(file=pdffn, width=12, height=12)

cols = setNames(c("violet", "blue","green3", "yellow", "orange2", "red"),sort(unique(tipdata)))
cols

# plot posterior density on the number of changes
par(mar=c(5.1,4.1,4.1,2.1),las=1)
plot(density(resER_smap), bty="l")

title(main="Posterior distribution of changes of each type", font.main=3)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)
