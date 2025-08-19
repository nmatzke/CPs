library(ape)
library(maps)
library(phytools)

setwd("~/GitHub/CPs/examples/CP_ace_pie/")

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
	
	# Run time for ML: FAILED (with these settings; make.simmap works though)
	# ML_7abprCTE = ape::ace(x=states, phy=tr, type="discrete", method="ML", model=seven_rates_abpr_complex_trap_evolution_model, marginal=FALSE, use.expm = TRUE, use.eigen=FALSE)

	
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
dim(ancstates)



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




pdf(file="plot_ACE_7abprCTE_v2.pdf", width=12, height=80)

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


# Nepentheaceae-Droseraceae (to see the transition from non-carnivorous to sticky traps)
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



#######################################################
# Version 3 graphic:
# Ladderize tree
# Extract key clades, with ancestral states
# plot those in reasonable fashion
# ??
#######################################################
library(BioGeoBEARS)   # for prt()

# Ladderizing a tree just makes it look nicer
ltr = ladderize(tr, right=FALSE)
# Make sure the nodes are relabeled correctly be re-loading tree from Newick
ltr = read.tree(file="", text=write.tree(ltr, file=""))

plot(tr, show.tip.label=FALSE, main="original pruned tree focusing on CPs")
plot(ltr, show.tip.label=FALSE, main="ladderized-leftwards tree")

# Make tree tables to convert the nodes
trtable = prt(tr)
ltrtable = prt(ltr)

# Key thing about trtables:
# the row number = the node number
# row 52 encodes the information for node #52 in the tr (according to the APE node numbering in R)

# Table to convert between node numbers
matchnodes_trNode_to_ltrNode = match(x=ltrtable$tipnames, table=trtable$tipnames)
matchnodes_trNode_to_ltrNode

matchnodes_ltrNode_to_trNode = match(x=trtable$tipnames, table=ltrtable$tipnames)
matchnodes_ltrNode_to_trNode

# Conversion table
convert_tr_to_ltr = cbind(1:length(matchnodes_convert_trNode_to_ltrNode),matchnodes_ltrNode_to_trNode)
convert_tr_to_ltr_df = as.data.frame(convert_tr_to_ltr, stringsAsFactors=FALSE)
names(convert_tr_to_ltr_df) = c("trNode","ltrNode")
head(convert_tr_to_ltr_df)
tail(convert_tr_to_ltr_df)

# Sort on ltr nodes to convert the other way
convert_ltr_to_tr_df = convert_tr_to_ltr_df[order(convert_tr_to_ltr_df$ltrNode),]
head(convert_ltr_to_tr_df)
tail(convert_ltr_to_tr_df)


# Check if it did what you thought (converting tr nodes to ltr nodes)
head(trtable[convert_ltr_to_tr_df$trNode,])
head(ltrtable)

tail(trtable[convert_ltr_to_tr_df$trNode,])
tail(ltrtable)

# This works, meaning:
# * The rows/nodes of trtable can be converted to the equivalent positions in ltr and ltrtable, using
#   matchnodes_convert_trNode_to_ltrNode

# The three big clades are:
ancnode_of_lca_Lentibulariaceae
ancnode_of_lca_NeDro
ancnode_of_lca_Sarrac

# Let's extract each one as a separate tree, with separate Ancestral States numbers
ancnode_of_lca_Lentibulariaceae # 2114

# This says, for the first node, place in the matchnodes_convert_trNode_to_ltrNode[1] position
# trtable[matchnodes_convert_trNode_to_ltrNode]

# So this node, should go in the matchnodes_convert_trNode_to_ltrNode[ancnode_of_lca_Lentibulariaceae] position to match ltrtable
ancnode_of_lca_Lentibulariaceae

ancnode_of_lca_Lentibulariaceae_in_ltr = convert_tr_to_ltr_df$ltrNode[ancnode_of_lca_Lentibulariaceae]

# Double-check, yes this looks correct
ltrtable[ancnode_of_lca_Lentibulariaceae_in_ltr, ]

# Get the Lenti subtree
lenti_subtree = extract.clade(phy=ltr, node=ancnode_of_lca_Lentibulariaceae_in_ltr)
lenti_subtree = read.tree(file="", text=write.tree(lenti_subtree, file=""))
sort(lenti_subtree$tip.label)


# OK, clades have been extracted to subtrees

# Let's subset the ancestral state results
ancnode_of_lca_Lentibulariaceae

# descendent node numbers in tr
ancnode_of_lca_Lentibulariaceae_daughters = BioGeoBEARS::get_daughter_nodes(nodenum=ancnode_of_lca_Lentibulariaceae, tr=tr, nodes=NULL)
ancnode_of_lca_Lentibulariaceae_daughters

# Ancstates in APE node order
ancstates2
dim(ancstates2)

# extract ancstates2 for Lentibulariaceae
# put in order of lenti_subtree

ancstates2_for_lenti_tr = ancstates2[ancnode_of_lca_Lentibulariaceae_daughters,]
dim(ancstates2_for_lenti_tr)

# Get the taxa lists descending from every node
nodes_to_keep_labels = trtable$tipnames[ancnode_of_lca_Lentibulariaceae_daughters]
nodes_to_keep_labels

# match to subtree
lenti_subtree_table = prt(lenti_subtree)
get_subtree_tipnames_in_lenti_order = match(lenti_subtree_table$tipnames, table=nodes_to_keep_labels)
ancstates2_for_lenti_subtree_ltr_order = ancstates2_for_lenti_tr[get_subtree_tipnames_in_lenti_order,]
ancstates2_for_lenti_subtree_ltr_order
dim(ancstates2_for_lenti_subtree_ltr_order)

# Test plot

pdffn = "ancstates2_for_lenti_subtree_ltr_order.pdf"
pdf(file=pdffn, width=20, height=35)

plot(lenti_subtree)
cols = c("white","lightblue","blue", "yellow","orange",
         "orange3","red","lightgrey", "darkgrey", "green3", "darkgreen")

# get the internal node numbers
tipnode_nums = 1:length(lenti_subtree$tip.label)
internal_nodenums = (length(lenti_subtree$tip.label)+1):(length(lenti_subtree$tip.label)+lenti_subtree$Nnode)

nodelabels(node=internal_nodenums, pie=ancstates2_for_lenti_subtree_ltr_order[internal_nodenums,], piecol=cols, cex=0.5)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)






#######################################################
# For Sarracenia
#######################################################
ancnode_of_lca_Sarraceniaceae_in_ltr = convert_tr_to_ltr_df$ltrNode[ancnode_of_lca_Sarrac]

# Double-check, yes this looks correct
ltrtable[ancnode_of_lca_Sarraceniaceae_in_ltr, ]

# Get the sarra subtree
sarra_subtree = extract.clade(phy=ltr, node=ancnode_of_lca_Sarraceniaceae_in_ltr)
sarra_subtree = read.tree(file="", text=write.tree(sarra_subtree, file=""))
sort(sarra_subtree$tip.label)


# OK, clades have been extracted to subtrees

# Let's subset the ancestral state results
ancnode_of_lca_Sarraceniaceae

# descendent node numbers in tr
ancnode_of_lca_Sarraceniaceae_daughters = BioGeoBEARS::get_daughter_nodes(nodenum=ancnode_of_lca_Sarrac, tr=tr, nodes=NULL)
ancnode_of_lca_Sarraceniaceae_daughters

# Ancstates in APE node order
ancstates2
dim(ancstates2)

# extract ancstates2 for Sarraceniaceae
# put in order of sarra_subtree

ancstates2_for_sarra_tr = ancstates2[ancnode_of_lca_Sarraceniaceae_daughters,]
dim(ancstates2_for_sarra_tr)

# Get the taxa lists descending from every node
nodes_to_keep_labels = trtable$tipnames[ancnode_of_lca_Sarraceniaceae_daughters]
nodes_to_keep_labels

# match to subtree
sarra_subtree_table = prt(sarra_subtree)
get_subtree_tipnames_in_sarra_order = match(sarra_subtree_table$tipnames, table=nodes_to_keep_labels)
ancstates2_for_sarra_subtree_ltr_order = ancstates2_for_sarra_tr[get_subtree_tipnames_in_sarra_order,]
ancstates2_for_sarra_subtree_ltr_order
dim(ancstates2_for_sarra_subtree_ltr_order)

# Test plot

pdffn = "ancstates2_for_sarra_subtree_ltr_order.pdf"
pdf(file=pdffn, width=12, height=25)

plot(sarra_subtree)
cols = c("white","lightblue","blue", "yellow","orange",
         "orange3","red","lightgrey", "darkgrey", "green3", "darkgreen")

# get the internal node numbers
tipnode_nums = 1:length(sarra_subtree$tip.label)
internal_nodenums = (length(sarra_subtree$tip.label)+1):(length(sarra_subtree$tip.label)+sarra_subtree$Nnode)

nodelabels(node=internal_nodenums, pie=ancstates2_for_sarra_subtree_ltr_order[internal_nodenums,], piecol=cols, cex=0.5)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)







#######################################################
# For Nepenthenales
#######################################################
ancnode_of_lca_Nepenthenales_in_ltr = convert_tr_to_ltr_df$ltrNode[ancnode_of_lca_NeDro]

# Double-check, yes this looks correct
ltrtable[ancnode_of_lca_Nepenthenales_in_ltr, ]

# Get the NeDro subtree
NeDro_subtree = extract.clade(phy=ltr, node=ancnode_of_lca_Nepenthenales_in_ltr)
NeDro_subtree = read.tree(file="", text=write.tree(NeDro_subtree, file=""))
sort(NeDro_subtree$tip.label)

# OK, clades have been extracted to subtrees

# Let's subset the ancestral state results
ancnode_of_lca_Nepenthales

# descendent node numbers in tr
ancnode_of_lca_Nepenthales_daughters = BioGeoBEARS::get_daughter_nodes(nodenum=ancnode_of_lca_NeDro, tr=tr, nodes=NULL)
ancnode_of_lca_Nepenthales_daughters

# Ancstates in APE node order
ancstates2
dim(ancstates2)

# extract ancstates2 for Nepenthales
# put in order of NeDro_subtree

ancstates2_for_NeDro_tr = ancstates2[ancnode_of_lca_Nepenthales_daughters,]
dim(ancstates2_for_NeDro_tr)

# Get the taxa lists descending from every node
nodes_to_keep_labels = trtable$tipnames[ancnode_of_lca_Nepenthales_daughters]
nodes_to_keep_labels

# match to subtree
NeDro_subtree_table = prt(NeDro_subtree)
get_subtree_tipnames_in_NeDro_order = match(NeDro_subtree_table$tipnames, table=nodes_to_keep_labels)
ancstates2_for_NeDro_subtree_ltr_order = ancstates2_for_NeDro_tr[get_subtree_tipnames_in_NeDro_order,]
ancstates2_for_NeDro_subtree_ltr_order
dim(ancstates2_for_NeDro_subtree_ltr_order)

# Test plot

pdffn = "ancstates2_for_NeDro_subtree_ltr_order.pdf"
pdf(file=pdffn, width=25, height=55)

plot(NeDro_subtree)
cols = c("white","lightblue","blue", "yellow","orange",
         "orange3","red","lightgrey", "darkgrey", "green3", "darkgreen")

# get the internal node numbers
tipnode_nums = 1:length(NeDro_subtree$tip.label)
internal_nodenums = (length(NeDro_subtree$tip.label)+1):(length(NeDro_subtree$tip.label)+NeDro_subtree$Nnode)

nodelabels(node=internal_nodenums, pie=ancstates2_for_NeDro_subtree_ltr_order[internal_nodenums,], piecol=cols, cex=0.5)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

