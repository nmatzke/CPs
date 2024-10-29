library(ape)
library(maps)
library(phytools)

setwd("/GitHub/CPs/examples/CP_ace_pie/")

trfn = "gbotb_tr14_wSister_genera_edited2_minusMonocots.newick"
states_fn = "gbotb_tr14_wSister_genera_edited_states_noMonocots_v2.txt"

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





###################################################################################
#7 rates (aerial: both-pitcher or reverse) complex trap evolution model#
###################################################################################
seven_rates_abpr_complex_trap_evolution_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
                                                              1,0,3,0,0,11,0,0,0,0,0,
                                                              1,3,0,4,0,0,9,0,0,0,0,
                                                              1,0,0,0,5,0,0,0,0,0,0,
                                                              1,0,0,0,0,0,0,0,0,0,0,
                                                              1,0,0,0,0,0,0,12,0,0,0,
                                                              1,0,13,0,0,0,0,0,10,0,0,
                                                              1,0,0,0,0,15,0,0,6,0,0,
                                                              1,0,0,0,0,0,14,6,0,7,0,
                                                              1,0,0,0,0,0,0,0,0,0,8,
                                                              1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

seven_rates_abpr_complex_trap_evolution_model











##################
###7abprCTE map###
##################

runslow = TRUE
if (runslow == TRUE)
	{
	set.seed(34321)
	stochastic_maps_7abprCTE = make.simmap(tree=tr, x=states, model=seven_rates_abpr_complex_trap_evolution_model, nsim=100)
	# Save the full stochastic maps
	save(stochastic_maps_7abprCTE, file="stochastic_maps_7abprCTE.Rdata")
	} else {
	# Loads to: stochastic_maps_7abprCTE
	load(file="stochastic7abprCTE.Rdata")	
	} # END runslow

stochastic_maps_7abprCTE

# Make and save a summary 
summary_stochastic_maps_7abprCTE = summary(stochastic_maps_7abprCTE)
save(summary_stochastic_maps_7abprCTE, file="summary_stochastic_maps_7abprCTE.Rdata")
load(file="summary_stochastic_maps_7abprCTE.Rdata")	



names(summary_stochastic_maps_7abprCTE)

summary_stochastic_maps_7abprCTE$count


# ace = ancestral character estimation
ancstates = summary_stochastic_maps_7abprCTE$ace

pdf(file="plot_ACE_7abprCTE_v2.pdf", width=12, height=40)

cols = c("white","lightblue","blue", "yellow","orange",
         "orange3","red","lightgrey", "darkgrey", "green3", "darkgreen")

# get the internal node numbers
tipnode_nums = 1:length(tr$tip.label)
internal_nodenums = (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)

# Cut to just carnivorous nodes
carnivorous_TF = ancstates[,"1"] < 0.50
#internal_nodenums_CPs = internal_nodenums[carnivorous_TF]
#ancstates_CPs = ancstates[carnivorous_TF,]



# Let's also get the node ancestor of the LCA of each clade
list_of_ancnodes_of_lca = NULL
# Lentibulariaceae
lca_of_Lentibulariaceae = getMRCA(phy=tr, tip=c("Pinguicula_alpina","Utricularia_albiflora"))
lca_of_Lentibulariaceae
edgeTF = tr$edge[,2] == lca_of_Lentibulariaceae
edgenum = (1:(length(edgeTF)))[edgeTF]
edgenum
ancnode_of_lca_Lentibulariaceae = tr$edge[edgeTF,1]
ancnode_of_lca_Lentibulariaceae
tr$edge[edgeTF,]


# Nepentheaceae-Droseraceae (to see the transtion from non-carnivorous to sticky traps)
lca_of_NeDro = getMRCA(phy=tr, tip = c("Drosera_lanata", "Nepenthes_gracilis"))
lca_of_NeDro

edgeTF_NeDro = tr$edge[,2] == lca_of_NeDro
edgeTF_NeDro
edgenum_NeDro = (1:(length(edgeTF_NeDro)))[edgeTF_NeDro]
edgenum_NeDro
ancnode_of_lca_NeDro = tr$edge[edgeTF_NeDro,1]
ancnode_of_lca_NeDro
tr$edge[edgeTF_NeDro,]


# Nepentheaceae (see transition from sticky to pitcher throught ancestral transitional state)
#trtable$tipnames[1583:1732]
lca_of_Nepentheaceae = getMRCA(phy=tr, tip = c("Nepenthes_madagascariensis", "Nepenthes_pervillei"))
lca_of_Nepentheaceae

edgeTF_Nepentheaceae = tr$edge[,2] == lca_of_Nepentheaceae
edgeTF_Nepentheaceae
edgenum_Nepentheaceae = (1:(length(edgeTF_Nepentheaceae)))[edgeTF_Nepentheaceae]
edgenum_Nepentheaceae
ancnode_of_lca_Nepentheaceae = tr$edge[edgeTF_Nepentheaceae,1]
ancnode_of_lca_Nepentheaceae
tr$edge[edgeTF_Nepentheaceae,]


#Sarraceniaceae
lca_of_Sarrac = getMRCA(phy=tr, tip = c("Sarracenia_rosea", "Darlingtonia_californica"))
lca_of_Sarrac

edgeTF_Sarrac = tr$edge[,2] == lca_of_Sarrac
edgeTF_Sarrac
edgenum_Sarrac = (1:(length(edgeTF_Sarrac)))[edgeTF_Sarrac]
edgenum_Sarrac
ancnode_of_lca_Sarrac = tr$edge[edgeTF_Sarrac,1]
ancnode_of_lca_Sarrac
tr$edge[edgeTF_Sarrac,]


#Droseraceae
lca_of_Droseraceae = getMRCA(phy=tr, tip = c("Drosera_arcturi", "Aldrovanda_vesiculosa"))
lca_of_Droseraceae

edgeTF_Droseraceae = tr$edge[,2] == lca_of_Droseraceae
edgeTF_Droseraceae
edgenum_Droseraceae = (1:(length(edgeTF_Droseraceae)))[edgeTF_Droseraceae]
edgenum_Droseraceae
ancnode_of_lca_Droseraceae = tr$edge[edgeTF_Droseraceae,1]
ancnode_of_lca_Droseraceae
tr$edge[edgeTF_Droseraceae,]


# Lentibulariaceae-transitional-pitcher-eel
lca_of_Lentibulariaceae_2 = getMRCA(phy=tr, tip=c("Genlisea_uncinata","Utricularia_uniflora"))
lca_of_Lentibulariaceae_2
edgeTF_2 = tr$edge[,2] == lca_of_Lentibulariaceae_2
edgenum_2 = (1:(length(edgeTF_2)))[edgeTF_2]
edgenum_2
ancnode_of_lca_Lentibulariaceae_2 = tr$edge[edgeTF_2,1]
ancnode_of_lca_Lentibulariaceae_2
tr$edge[edgeTF_2,]


list_of_ancnodes_of_lca = c(ancnode_of_lca_Lentibulariaceae, ancnode_of_lca_NeDro, ancnode_of_lca_Nepentheaceae, ancnode_of_lca_Sarrac, ancnode_of_lca_Droseraceae, ancnode_of_lca_Lentibulariaceae_2)

print(list_of_ancnodes_of_lca)

internal_nodenums_CPs = c(internal_nodenums[carnivorous_TF], list_of_ancnodes_of_lca)
#internal_nodenums_CPs = internal_nodenums[carnivorous_TF]


ancstates_CPs = rbind(ancstates[carnivorous_TF,], ancstates[list_of_ancnodes_of_lca,])
#ancstates_CPs = ancstates[carnivorous_TF,]

plot.phylo(tr, show.tip.label=FALSE)

nodelabels(node=internal_nodenums_CPs, pie=ancstates_CPs, piecol=cols, cex=0.5)

dev.off()
system("open plot_ACE_7abprCTE_v2.pdf")