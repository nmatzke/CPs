#######################################################
# Example Ancestral Character Estimation ("ace")
# on Carnivorous Plant Trap traits
#######################################################

library(ape)
library(phytools)


wd = "/GitHub/CPs/examples/drop.tip_v1/"
setwd(wd)

trfn = "gbotb_tr14_wSister_genera_edited2.newick"
tr = read.tree(trfn)

# Plot this massive tree, but WITHOUT tipnames,
# because there are sooooo many tipnames the graphic
# stalls or crashes.
plot.phylo(tr, show.tip.label=FALSE)



# To reduce the tree by cutting out monocots, we have to:
# 1. Identify the MRCA node (most recent common ancestor) of
#    the clade we want to drop from the tree
#    a. This will require finding "tip specifier" names
# 2. Identify its node number
# 3. Put this into ape::drop.tip command

tip_specifiers = c("Typha_minima", "Spathanthus_bicolor")

nodenum = getMRCA(phy=tr, tip=tip_specifiers)
nodenum

# We need a list of ALL of the tips to drop
# This will be all of the species that descend from node
# 4528

# One way to get this list, is with extract.clade

monocots_clade_subtree = extract.clade(phy=tr, node=nodenum)
monocots_clade_subtree

monocots_clade_tipnames = monocots_clade_subtree$tip.label
monocots_clade_tipnames

tr2 = drop.tip(phy=tr, tip=monocots_clade_tipnames)
tr2 

# Compare number of tips
length(tr$tip.label)
length(tr2$tip.label)

# Output tree
outfn = "gbotb_tr14_wSister_genera_edited2_minusMonocots.newick"
write.tree(tr2, file=outfn)




#######################################################
# Drop the same tips from the states file
#######################################################
states_fn = "gbotb_tr14_wSister_genera_edited_states.txt"

ddf = read.table(states_fn, header=TRUE, sep="\t")
head(ddf)
dim(ddf)

# Ensure appropriate sorting
states = ddf$state
states = c(as.numeric(states) + 1)
names(states) = ddf$species
head(states)

summary(as.factor(states))


# Subset the states to remove monocots_clade_tipnames

# Calculate TRUE/FALSE (TF) on whether a states name is in 
# the monocots list

TF = ddf$species %in% monocots_clade_tipnames
TF
length(TF)
sum(TF)

# Subset the states table
# Square brackets is always [rows, columns]
keepTF = TF == FALSE # keep the non-monocots
ddf2 = ddf[keepTF,]

dim(ddf)
dim(ddf2)

# Save to output file
outdata_fn = "gbotb_tr14_wSister_genera_edited_states_noMonocots.txt"
write.table(x=ddf2, file=outdata_fn, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)





