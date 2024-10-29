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

# NJM
runslow = FALSE
if (runslow == TRUE)
	{
	set.seed(34321)
	
	# 1 minute to run
	trtable = BioGeoBEARS::prt(tr)
	save(trtable, file="gbotb_tr14_wSister_genera_edited2_minusMonocots_trtable.Rdata")
	
	# Run time:
	ML_7abprCTE = ape::ace(x=states, phy=tr, type="discrete", method="ML", model=seven_rates_abpr_complex_trap_evolution_model, marginal=FALSE, use.expm = TRUE, use.eigen=FALSE)

# ML_7abprCTE = phytools::fitMK(tree=tr, x=states, model=seven_rates_abpr_complex_trap_evolution_model, pi="equal")


	
	# ~1 hour to run
	stochastic_maps_7abprCTE = make.simmap(tree=tr, x=states, model=seven_rates_abpr_complex_trap_evolution_model, nsim=100)
	# Save the full stochastic maps
	save(stochastic_maps_7abprCTE, file="stochastic_maps_7abprCTE.Rdata")

	# Make and save a summary as well 
	# several minutes
	summary_stochastic_maps_7abprCTE = summary(stochastic_maps_7abprCTE)
	save(summary_stochastic_maps_7abprCTE, file="summary_stochastic_maps_7abprCTE.Rdata")
	} else {
	# Loads to: trtable
	load(file="gbotb_tr14_wSister_genera_edited2_minusMonocots_trtable.Rdata")


	# Loads to: stochastic_maps_7abprCTE
	load(file="stochastic_maps_7abprCTE.Rdata")	
	# Loads to: stochastic_maps_7abprCTE
	load(file="summary_stochastic_maps_7abprCTE.Rdata")	
	} # END runslow


summary_stochastic_maps_7abprCTE$count
names(summary_stochastic_maps_7abprCTE)


# ace = ancestral character estimation
ancstates = summary_stochastic_maps_7abprCTE$ace


# lca_of_Lentibulariaceae
2205

# ancnode_of_lca_Lentibulariaceae
2114

ancstates[2205,]
ancstates[2114,]

# Row of prt table
trtable[2205,]
#    node ord_ndname node_lvl node.type parent_br edge.length ancestor daughter_nds node_ht  time_bp fossils label
#2205 2205       2205       12  internal       640    10.53411     2114   2206, 2285 90.8475 32.88674      NA      

# edge/branch number:
640

stochastic_maps_7abprCTE[[1]]$maps[640]
stochastic_maps_7abprCTE[[2]]$maps[640]
stochastic_maps_7abprCTE[[3]]$maps[640]
stochastic_maps_7abprCTE[[4]]$maps[640]
stochastic_maps_7abprCTE[[5]]$maps[640]




pdf(file="plot_ACE_7abprCTE_v2.pdf", width=12, height=40)

cols = c("white","lightblue","blue", "yellow","orange",
         "orange3","red","lightgrey", "darkgrey", "green3", "darkgreen")

# get the internal node numbers
tipnode_nums = 1:length(tr$tip.label)
internal_nodenums = (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)

# Cut to just carnivorous nodes
carnivorous_TF = ancstates[,"1"] < 0.50
CPnums = (1:length(carnivorous_TF))[carnivorous_TF]
CPnums
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


# NJM:
# Let's make sure the secondarily noncarnivorous relatives
# of Droseraceae and Nepenthes are in there:
TF = grepl(tr$tip.label, pattern="Triphyophyllum")
num = (1:length(TF))[TF]
Triphyophyllum_num = num
tr$tip.label[1725:1775]
# Ancistrocladus_likoko
# Habropetalum_dawei
# Triphyophyllum_peltatum
getMRCA(phy=tr, tip=c("Ancistrocladus_likoko", "Habropetalum_dawei"))
getMRCA(phy=tr, tip=c("Ancistrocladus_likoko", "Triphyophyllum_peltatum"))
getMRCA(phy=tr, tip=c("Habropetalum_dawei", "Triphyophyllum_peltatum"))

nodes_of_nonCarnivorous_Droseraceae_relatives = BioGeoBEARS::get_daughter_nodes(nodenum=3615, tr=tr, nodes=NULL)
nodes_of_nonCarnivorous_Droseraceae_relatives

# Remove Triphyophyllum_num & put last, so it plots better
nodes_of_nonCarnivorous_Droseraceae_relatives = nodes_of_nonCarnivorous_Droseraceae_relatives[nodes_of_nonCarnivorous_Droseraceae_relatives != Triphyophyllum_num]
nodes_of_nonCarnivorous_Droseraceae_relatives = c(nodes_of_nonCarnivorous_Droseraceae_relatives, Triphyophyllum_num)


list_of_ancnodes_of_lca = c(ancnode_of_lca_Lentibulariaceae, ancnode_of_lca_NeDro, ancnode_of_lca_Nepentheaceae, ancnode_of_lca_Sarrac, ancnode_of_lca_Droseraceae, ancnode_of_lca_Lentibulariaceae_2, nodes_of_nonCarnivorous_Droseraceae_relatives)

# NJM: These ancnode numbers are for when the whole tree is labeled.
# However, the ancstates_CPs table is just internal nodes
# The way APE works is:
# tip_nodes = 1:length(tr$tip.label)
#  = 1-1879
# internal_nodes = (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)
#  = 1880-3757

# NJM compare:
length(internal_nodenums)
length(carnivorous_TF)

# Also weirdly, "ancstates" from summary() on a simmaps object
# ancstates = summary_stochastic_maps_7abprCTE$ace
# has internal nodes first, THEN tip nodes
ancstates[1:50,]
ancstates[1850:1900,] # <- changes to tips here
ancstates[3700:3717,]

# Same when you subset by carnivorous_TF
ancstates[carnivorous_TF,][1:50,]
ancstates[carnivorous_TF,][410:440,] # <- changes to tips here
ancstates[carnivorous_TF,][830:860,]


# These issues lead to this having NA at end:
internal_nodenums[carnivorous_TF]


# NJM SOLUTION: re-order ancstates, makes ancstates2 with
# standard APE node order
tip_nodes = 1:length(tr$tip.label)
internal_nodes = (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)
all_nodes = c(tip_nodes, internal_nodes)
ancstates_internal = summary_stochastic_maps_7abprCTE$ace[1:length(internal_nodes),]
ancstates_tips = summary_stochastic_maps_7abprCTE$ace[(length(internal_nodes)+1):nrow(summary_stochastic_maps_7abprCTE$ace),]
# These should have APE order
ancstates2 = rbind(ancstates_tips, ancstates_internal)

# Check ancstates2
ancstates2[1:50,]
ancstates2[1850:1900,] # <- changes to INTERNAL here
ancstates2[3700:3717,]


# It should work with ancstates2 and all_nodes
carnivorous_TF2 = ancstates2[,"1"] < 0.50

internal_nodenums_CPs = c(all_nodes[carnivorous_TF2], list_of_ancnodes_of_lca)

ancstates_CPs2 = rbind(ancstates2[carnivorous_TF2,], ancstates2[list_of_ancnodes_of_lca,])

ancstates2[carnivorous_TF,][1:100,]
ancstates2[carnivorous_TF,][101:200,]
ancstates2[carnivorous_TF,][201:300,]
ancstates2[carnivorous_TF,][301:400,]
ancstates2[carnivorous_TF,][401:500,]
ancstates2[carnivorous_TF,][501:600,]
ancstates2[carnivorous_TF,][601:700,]
ancstates2[carnivorous_TF,][701:800,]
ancstates2[carnivorous_TF,][750:860,]



#ancstates_CPs = ancstates[carnivorous_TF,]

plot.phylo(tr, show.tip.label=FALSE)

nodelabels(node=internal_nodenums_CPs, pie=ancstates_CPs2, piecol=cols, cex=0.5)




# NJM 
# double check LCA of Lentibuariaceae
# 2205
# nodelabels(text=lca_of_Lentibulariaceae, node=lca_of_Lentibulariaceae, cex=0.3)

# Check ancstates
# ancstates[lca_of_Lentibulariaceae,]



dev.off()
system("open plot_ACE_7abprCTE_v2.pdf")



