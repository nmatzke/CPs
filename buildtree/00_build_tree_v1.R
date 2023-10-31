# devtools::install_github("jinyizju/V.PhyloMaker2")

library(ape)
library(BioGeoBEARS)
library(V.PhyloMaker2)

source("/GitHub/bioinfRhints/R/r8s_functions_v2.R")
source("/GitHub/bioinfRhints/R/_R_tree_functions_v2.R")

wd = "/GitHub/CPs/buildtree/"
setwd(wd)

#?at.node


#######################################################
# Accessing a large, pre-calculated angiosperm phylogeny
#######################################################

tr = GBOTB.extended.WP
sum(grepl(pattern="Nepenthes", x=tr$tip.label))
sum(grepl(pattern="Utricularia", x=tr$tip.label))
sum(grepl(pattern="Pinguicula", x=tr$tip.label))
sum(grepl(pattern="Brocchinia", x=tr$tip.label))
sum(grepl(pattern="Catopsis", x=tr$tip.label))
#sum(grepl(pattern="Brocchinia_hechtioides", x=tr$tip.label))
sum(grepl(pattern="Brocchinia_reducta", x=tr$tip.label))
#sum(grepl(pattern="Catopsis_berteroniana", x=tr$tip.label))

# Cephalotus tip
TF = grepl(pattern="Cephalotus", x=tr$tip.label)
tr$tip.label[TF]

# Tip info
TF = grepl(pattern="Cephalotus", x=tips.info.LCVP$species)
tip_table = tips.info.LCVP[TF,c("species","genus","family")]
tip_table
build.nodes.2(tree=tr, tips=tip_table)

TF = grepl(pattern="Drosera", x=tips.info.LCVP$species)
tip_table = tips.info.LCVP[TF,c("species","genus","family")]
tip_table
build.nodes.2(tree=tr, tips=tip_table)


#######################################################
# Checking all the genera
#######################################################

# CP genera:
genera = c("Brocchinia_reducta",
#"Brocchinia_hechtioides",
#"Catopsis_berteroniana",
"Drosera",
"Dionaea",
"Aldrovanda",
"Drosophyllum",
"Triphyophyllum",
"Roridula",
"Sarracenia",
"Heliamphora",
"Darlingtonia",
"Nepenthes",
"Genlisea",
"Pinguicula",
"Utricularia",
"Byblis",
"Cephalotus")

pdffn = "GBOTB.extended.WP_subtrees.pdf"
pdf(file=pdffn, width=12, height=12)


num_species = NULL
root_age = NULL
node_ages = NULL
for (i in 1:length(genera))
	{
	TF = grepl(pattern=genera[i], x=tr$tip.label)
	num_species = c(num_species, sum(TF))
	
	tips_to_keep = tr$tip.label[TF]
	node_to_keep = getMRCA(phy=tr, tip=tips_to_keep)
	if (sum(TF) > 1)
		{
		subtree = extract.clade(phy=tr, node_to_keep)
		subtr_age = get_max_height_tree(subtree)
		
		plot(subtree)
		axisPhylo()
		title(genera[i])

		
		} else {
		subtr_age = 0.0
		}
	root_age = c(root_age, subtr_age)
	
	TF = grepl(pattern=genera[i], x=tips.info.LCVP$species)
	tip_table = tips.info.LCVP[TF,c("species","genus","family")]
	tip_table
	tmptable = build.nodes.2(tree=tr, tips=tip_table)
	node_ages = rbind(node_ages, tmptable)
	}


dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


sdf = as.data.frame(cbind(genera, root_age, num_species), stringsAsFactors=FALSE)
sdf

node_ages


pdffn = "Masafumi_trees.pdf"
pdf(file=pdffn, width=12, height=12)



#######################################################
# Masafumi's trees
#######################################################
# Dated
Droseraceae_trfn = "/GitHub/CPs/rawtrees/Digitized_trees/Droseraceae/Sen_2020_Fig1_Droseraceae_digitized_v1.newick"
Droseraceae_tr = read.tree(Droseraceae_trfn)
plot(Droseraceae_tr)
axisPhylo()
title("Droseraceae")
get_treeheight(Droseraceae_tr)


Genlisea_trfn = "/GitHub/CPs/rawtrees/Digitized_trees/Lentibulaceae/Genlisea/Fleischmann_2010_Fig5_Genlisea_digitized_v1.newick"
Genlisea_tr = read.tree(Genlisea_trfn)
plot(Genlisea_tr)
axisPhylo()
title("Genlisea")
get_treeheight(Genlisea_tr)

# dated
Pinguicula1_trfn = "/GitHub/CPs/rawtrees/Digitized_trees/Lentibulaceae/Pinguicula/Silva_et_al_2018_Fig7_Pinguicula_digitized_v1.newick"
Pinguicula1_tr = read.tree(Pinguicula1_trfn)
plot(Pinguicula1_tr)
axisPhylo()
title("Pinguicula1")
get_treeheight(Pinguicula1_tr)

Pinguicula2_trfn = "/GitHub/CPs/rawtrees/Digitized_trees/Lentibulaceae/Pinguicula/Shimai_et_al_2021_S1_Fig1_Pinguicula_digitized_v1.newick"
Pinguicula2_tr = read.tree(Pinguicula2_trfn)
plot(Pinguicula2_tr)
axisPhylo()
title("Pinguicula2")
get_treeheight(Pinguicula2_tr)


# dated
Utricularia_trfn = "/GitHub/CPs/rawtrees/Digitized_trees/Lentibulaceae/Utricularia/Jobson_2017_Fig3_Utricularia_digitized_v1.newick"
Utricularia_tr = read.tree(Utricularia_trfn)
plot(Utricularia_tr)
axisPhylo()
title("Utricularia")
get_treeheight(Utricularia_tr)

Nepenthes_trfn = "/GitHub/CPs/rawtrees/Digitized_trees/Nepenthaceae/Murphy_2020_Fig4_Nepenthaceae_digitized_v1.newick"
Nepenthes_tr = read.tree(Nepenthes_trfn)
Nepenthes_tr$tip.label = gsub(pattern="N.", replacement="Nepenthes_", x=Nepenthes_tr$tip.label)
plot(Nepenthes_tr)
axisPhylo()
title("Nepenthes")
get_treeheight(Nepenthes_tr)


# dated
Heliamphora_trfn = "/GitHub/CPs/rawtrees/Digitized_trees/Sarraceniaceae/Heliamphora/Liu_Smith_2020_Fig1_Heliamphora_digitized_v1.newick"
Heliamphora_tr = read.tree(Heliamphora_trfn)
plot(Heliamphora_tr)
axisPhylo()
title("Heliamphora")
get_treeheight(Heliamphora_tr)

# dated
Sarracenia_trfn = "/GitHub/CPs/rawtrees/Digitized_trees/Sarraceniaceae/Sarracenia/Ellison_2012_Fig4_Sarracenia_digitized_v1.newick"
Sarracenia_tr = read.tree(Sarracenia_trfn)
plot(Sarracenia_tr)
axisPhylo()
title("Sarracenia")
get_treeheight(Sarracenia_tr)



dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)



#######################################################
# Query the digitized trees
#######################################################

genera2 = list()
tree2 = list()

genera2[["Drosera"]] = c("Drosera")
tree2[["Drosera"]] = Droseraceae_tr

genera2[["Droseraceae"]] = c("Drosera",
"Dionaea",
"Aldrovanda")
tree2[["Droseraceae"]] = Droseraceae_tr

genera2[["DionAldro"]] = c("Dionaea",
"Aldrovanda")
tree2[["DionAldro"]] = Droseraceae_tr

genera2[["Drosophyllaceae"]] = c("Drosera",
"Dionaea",
"Aldrovanda",
"Drosophyllum")
tree2[["Drosophyllaceae"]] = Droseraceae_tr

genera2[["Dioncophyllaceae"]] = c("Nepenthes",
"Triphyophyllum")
tree2[["DroseraDioncophyllaceaeceae"]] = GBOTB.extended.WP


genera2[["DrosoNep"]] = c("Drosera",
"Dionaea",
"Aldrovanda",
"Drosophyllum",
"Nepenthes",
"Triphyophyllum")
tree2[["DrosoNep"]] = GBOTB.extended.WP

genera2[["Nepenthes"]] = c("Nepenthes")
tree2[["Nepenthes"]] = Nepenthes_tr

genera2[["Sarracenia"]] = c("Sarracenia")
tree2[["Sarracenia"]] = Sarracenia_tr

genera2[["Sarraceniacae"]] = c("Sarracenia", "Heliamphora", "Darlingtonia", "Actinidia", "Roridula")
tree2[["Sarraceniacae"]] = Sarracenia_tr

genera2[["Heliamphora"]] = c("Heliamphora")
tree2[["Heliamphora"]] = Heliamphora_tr

genera2[["Genlisea"]] = c("Genlisea")
tree2[["Genlisea"]] = Genlisea_tr

genera2[["Pinguicula"]] = c("Pinguicula")
tree2[["Pinguicula"]] = Pinguicula1_tr

genera2[["Utricularia"]] = c("Utricularia")
tree2[["Utricularia"]] = Utricularia_tr

genera2[["Byblis"]] = c("Byblis")
tree2[["Byblis"]] = GBOTB.extended.WP


pdffn = "Masafumi_subtrees2.pdf"
pdf(file=pdffn, width=12, height=12)


names2 = NULL
root_age2 = NULL
num_species2 = NULL
for (i in 1:length(genera2))
	{
	names2 = c(names2, names(genera2)[i])

	TF = rep(FALSE, times=length(tree2[[i]]$tip.label))
	for (j in 1:length(genera2[[i]]))
		{
		tmpTF = grepl(pattern=genera2[[i]][j], x=tree2[[i]]$tip.label)
		TF = TF + tmpTF
		}
	num_species2 = c(num_species2, sum(TF))
	
	tips_to_keep = tree2[[i]]$tip.label[TF == 1]
	node_to_keep = getMRCA(phy=tree2[[i]], tip=tips_to_keep)
	if (sum(TF) > 1)
		{
		subtree = extract.clade(phy=tree2[[i]], node_to_keep)
		subtr_age = get_max_height_tree(subtree)

		plot(subtree)
		axisPhylo()
		title(names(genera2)[i])

		} else {
		subtr_age = 0.0
		}
	root_age2 = c(root_age2, subtr_age)
	}

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


#######################################################
# Query the master tree
#######################################################

pdffn = "GBOTB.extended.WP_subtrees2.pdf"
pdf(file=pdffn, width=12, height=12)

names2 = NULL
root_age3 = NULL
num_species3 = NULL
for (i in 1:length(genera2))
	{
	names2 = c(names2, names(genera2)[i])

	TF = rep(FALSE, times=length(GBOTB.extended.WP$tip.label))
	for (j in 1:length(genera2[[i]]))
		{
		tmpTF = grepl(pattern=genera2[[i]][j], x=GBOTB.extended.WP$tip.label)
		TF = TF + tmpTF
		}
	num_species3 = c(num_species3, sum(TF))
	
	tips_to_keep = GBOTB.extended.WP$tip.label[TF == 1]
	node_to_keep = getMRCA(phy=GBOTB.extended.WP, tip=tips_to_keep)
	if (sum(TF) > 1)
		{
		subtree = extract.clade(phy=GBOTB.extended.WP, node_to_keep)
		subtr_age = get_max_height_tree(subtree)
		
		plot(subtree)
		axisPhylo()
		title(names(genera2)[i])

		
		} else {
		subtr_age = 0.0
		}
	root_age3 = c(root_age3, subtr_age)
	}

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


sdf2 = as.data.frame(cbind(names2, root_age2, root_age3, num_species2, num_species3), stringsAsFactors=FALSE)
sdf2


get_treeheight(Droseraceae_tr)
get_treeheight(Genlisea_tr)
get_treeheight(Pinguicula1_tr)
get_treeheight(Pinguicula2_tr)
get_treeheight(Utricularia_tr)
get_treeheight(Nepenthes_tr)
get_treeheight(Heliamphora_tr)
get_treeheight(Sarracenia_tr)





# Manual modifications to tree:

# Droseraceae
masa_time_bp = 42.72850428
gbotb = 63.040034

Droseraceae_tr

# DrosoNep / Drosophyllaceae
masa_time_bp = 51.2285
gbotb = 82.989449

# Nepenthes
masa_time_bp = NA
gbotb = 25.701614
Nepenthes_tr









#######################################################
# Play with drop.tip
#######################################################
x=drop.tip(Pinguicula1_tr, tip=c("Pinguicula_caerulea", "Pinguicula_grandiflora", "Pinguicula_lusitanica"), trim.internal=FALSE, subtree=TRUE)
plot(x)
x$tip.label








#######################################################
# Time-scale the Nepenthes tree, and add to the main tree
#######################################################
# Requirements:
# * source() the file below 
# * have r8s installed somewhere like /Applications/r8s
#

calibration_node_tip_specifiers = c("Nepenthes_pervillei", "Nepenthes_madagascariensis")
calibration_age = 25.701614
r8s_nexus_fn = "Nepenthes_r8s_nexus_fn.nex"
r8s_logfile_fn = "Nepenthes_r8s_nexus_fn.log"
run_r8s_1calib(tr=Nepenthes_tr, calibration_node_tip_specifiers=calibration_node_tip_specifiers, r8s_method="LF", addl_cmd="", calibration_age=calibration_age, nsites=1000, tmpwd=getwd(), r8s_nexus_fn=r8s_nexus_fn, r8s_logfile_fn=r8s_logfile_fn, r8s_path="/Applications/r8s")

Nepenthes_tr_dated = extract_tree_from_r8slog(logfn=r8s_logfile_fn, delimiter=" = ", printall=TRUE)
Nepenthes_tr_rates = extract_rates_from_r8slog(logfn=r8s_logfile_fn)

pdffn = "Nepenthes_tr_dated.pdf"
pdf(file=pdffn, width=12, height=12)

plot(Nepenthes_tr_dated, cex=0.5)
axisPhylo()

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


get_treeheight(Nepenthes_tr_dated)

# Add to main tree
gbotb_tr = GBOTB.extended.WP
length(gbotb_tr$tip.label)


TF = grepl(pattern="Nepenthes", x=gbotb_tr$tip.label)
sum(TF)
tips_to_drop = gbotb_tr$tip.label[TF]
gbotb_tr2 = drop.tip(phy=gbotb_tr, tip=tips_to_drop, trim.internal=FALSE, subtree=TRUE)
length(gbotb_tr2$tip.label)
TF = grepl(pattern="\\[", x=gbotb_tr2$tip.label)
tipnum = (1:length(gbotb_tr2$tip.label))[TF]
newtip = gbotb_tr2$tip.label[TF]
newtip

# New tip is "[97_tips]"
gbotb_tr3 = bind.tree(x=gbotb_tr2, y=Nepenthes_tr_dated, where=tipnum)
length(gbotb_tr3$tip.label)

is.ultrametric(gbotb_tr)
is.ultrametric(gbotb_tr2)
is.ultrametric(gbotb_tr3)








#######################################################
# Time-scale the Droseraceae tree, and add to the main tree
#######################################################
# Just do linear scaling, as it's ultrametric
source("/GitHub/bioinfRhints/R/r8s_functions_v2.R")

# Subset to remove Nepenthes
TF = grepl(pattern="Nepenthes", x=Droseraceae_tr$tip.label)
sum(TF)
tips_to_drop = Droseraceae_tr$tip.label[TF]
Droseraceae_tr2 = drop.tip(phy=Droseraceae_tr, tip=tips_to_drop, trim.internal=TRUE, subtree=FALSE)
length(Droseraceae_tr$tip.label)
length(Droseraceae_tr2$tip.label)
is.ultrametric(Droseraceae_tr)
is.ultrametric(Droseraceae_tr2)
get_treeheight(Droseraceae_tr)
get_treeheight(Droseraceae_tr2)

# Branch length multiplier
calibration_age = 63.040034
multiplier = calibration_age / get_treeheight(Droseraceae_tr2)

Droseraceae_tr3 = Droseraceae_tr2
Droseraceae_tr3$edge.length = Droseraceae_tr2$edge.length * multiplier


# Add to main tree
gbotb_tr = gbotb_tr3
length(gbotb_tr$tip.label)


TF1 = grepl(pattern="Drosera", x=gbotb_tr$tip.label)
TF2 = grepl(pattern="Dionaea", x=gbotb_tr$tip.label)
TF3 = grepl(pattern="Aldrovanda", x=gbotb_tr$tip.label)
TF = (TF1 + TF2 + TF3) > 0
sum(TF)
tips_to_drop = gbotb_tr$tip.label[TF]
gbotb_tr2 = drop.tip(phy=gbotb_tr, tip=tips_to_drop, trim.internal=FALSE, subtree=TRUE)
length(gbotb_tr2$tip.label)
TF = grepl(pattern="\\[", x=gbotb_tr2$tip.label)
tipnum = (1:length(gbotb_tr2$tip.label))[TF]
newtip = gbotb_tr2$tip.label[TF]
newtip

# New tip is "[69_tips]"
gbotb_tr3 = bind.tree(x=gbotb_tr2, y=Droseraceae_tr3, where=tipnum)
length(gbotb_tr3$tip.label)

is.ultrametric(gbotb_tr)
is.ultrametric(gbotb_tr2)
is.ultrametric(gbotb_tr3)








#######################################################
# Time-scale the Sarracenia tree, and add to the main tree
#######################################################
# Just do linear scaling, as it's ultrametric
source("/GitHub/bioinfRhints/R/r8s_functions_v2.R")

# Subset to remove Roridula, Actinidia
TF1 = grepl(pattern="Roridula", x=Sarracenia_tr$tip.label)
TF2 = grepl(pattern="Actinidia", x=Sarracenia_tr$tip.label)
TF3 = grepl(pattern="Clethra", x=Sarracenia_tr$tip.label)
TF4 = grepl(pattern="Cyrilla", x=Sarracenia_tr$tip.label)
TF = (TF1 + TF2 + TF3 + TF4) > 0
sum(TF)
tips_to_drop = Sarracenia_tr$tip.label[TF]
Sarracenia_tr2 = drop.tip(phy=Sarracenia_tr, tip=tips_to_drop, trim.internal=TRUE, subtree=FALSE)
length(Sarracenia_tr$tip.label)
length(Sarracenia_tr2$tip.label)
is.ultrametric(Sarracenia_tr)
is.ultrametric(Sarracenia_tr2)
get_treeheight(Sarracenia_tr)
get_treeheight(Sarracenia_tr2)

# Branch length multiplier
calibration_age = 41.25056
multiplier = calibration_age / get_treeheight(Sarracenia_tr2)

Sarracenia_tr3 = Sarracenia_tr2
Sarracenia_tr3$edge.length = Sarracenia_tr2$edge.length * multiplier


# Add to main tree
gbotb_tr = gbotb_tr3
length(gbotb_tr$tip.label)


TF1 = grepl(pattern="Sarracenia", x=gbotb_tr$tip.label)
TF2 = grepl(pattern="Heliamphora", x=gbotb_tr$tip.label)
TF3 = grepl(pattern="Darlingtonia", x=gbotb_tr$tip.label)
TF = (TF1 + TF2 + TF3) > 0
sum(TF)
tips_to_drop = gbotb_tr$tip.label[TF]
gbotb_tr2 = drop.tip(phy=gbotb_tr, tip=tips_to_drop, trim.internal=FALSE, subtree=TRUE)
length(gbotb_tr2$tip.label)
TF = grepl(pattern="\\[", x=gbotb_tr2$tip.label)
tipnum = (1:length(gbotb_tr2$tip.label))[TF]
newtip = gbotb_tr2$tip.label[TF]
newtip

# New tip is "[15_tips]"
gbotb_tr3 = bind.tree(x=gbotb_tr2, y=Sarracenia_tr3, where=tipnum)
length(gbotb_tr3$tip.label)

length(gbotb_tr$tip.label)
length(gbotb_tr2$tip.label)
length(gbotb_tr3$tip.label)


is.ultrametric(gbotb_tr)
is.ultrametric(gbotb_tr2)
is.ultrametric(gbotb_tr3)



#######################################################
# Manually add Roridula_dentata as sister to Roridula_gorgonias at ~11 ma
#######################################################
Sarracenia_trtable = prt(Sarracenia_tr)
Sarracenia_trtable[2,]
Sarracenia_trtable$edge.length[2]
tip_branch_length = 10.52966

fake_trtxt = "(Roridula_dentata:10.52966,Roridula_dentata2:10.52966);"
faketr = read.tree(file="", text=fake_trtxt)
faketr
plot(faketr)
axisPhylo()

TF = grepl(pattern="Roridula_gorgonias", x=gbotb_tr3$tip.label)
tipnum = (1:length(gbotb_tr3$tip.label))[TF]
tipnum

gbotb_tr4 = bind.tree(x=gbotb_tr3, y=faketr, where=tipnum, position=tip_branch_length)

# Drop the extra fake tip
TF = grepl(pattern="Roridula_dentata2", x=gbotb_tr4$tip.label)
tips_to_drop = (1:length(gbotb_tr4$tip.label))[TF]
tips_to_drop

gbotb_tr5 = drop.tip(phy=gbotb_tr4, tip=tips_to_drop, trim.internal=TRUE, subtree=FALSE)

length(gbotb_tr3$tip.label)
length(gbotb_tr4$tip.label)
length(gbotb_tr5$tip.label)

is.ultrametric(gbotb_tr3)
is.ultrametric(gbotb_tr4)
is.ultrametric(gbotb_tr5)



#######################################################
# Double-check addition of Roridula_dentata
#######################################################
TF1 = grepl(pattern="Roridula", x=gbotb_tr5$tip.label)
TF2 = grepl(pattern="Actinidia", x=gbotb_tr5$tip.label)
TF3 = grepl(pattern="Clematoclethra", x=gbotb_tr5$tip.label)
#TF4 = grepl(pattern="Cyrilla", x=gbotb_tr5$tip.label)
TF4 = grepl(pattern="Saurauia", x=gbotb_tr5$tip.label)
TF5 = grepl(pattern="Sarracenia", x=gbotb_tr5$tip.label)
TF6 = grepl(pattern="Heliamphora", x=gbotb_tr5$tip.label)
TF7 = grepl(pattern="Darlingtonia", x=gbotb_tr5$tip.label)

TF = (TF1 + TF2 + TF3 + TF4 + TF5 + TF6 + TF7) > 0
sum(TF)
tips_to_keep = gbotb_tr5$tip.label[TF]
node_to_keep = getMRCA(phy=gbotb_tr5, tip=tips_to_keep)
Sarraceniaceae_gbotb_tr5 = extract.clade(phy=gbotb_tr5, node=node_to_keep)
plot(Sarraceniaceae_gbotb_tr5)
axisPhylo()





#######################################################
# Time-scale the Heliamphora tree, and add to the main tree
#######################################################
# Just do linear scaling, as it's ultrametric
source("/GitHub/bioinfRhints/R/r8s_functions_v2.R")

# Subset to remove Nepenthes
TF1 = grepl(pattern="Roridula", x=Heliamphora_tr$tip.label)
TF2 = grepl(pattern="Actinidia", x=Heliamphora_tr$tip.label)
TF3 = grepl(pattern="Darlingtonia", x=Heliamphora_tr$tip.label)
TF4 = grepl(pattern="Sarracenia", x=Heliamphora_tr$tip.label)
TF = (TF1 + TF2 + TF3 + TF4) > 0
sum(TF)
tips_to_drop = Heliamphora_tr$tip.label[TF]
Heliamphora_tr2 = drop.tip(phy=Heliamphora_tr, tip=tips_to_drop, trim.internal=TRUE, subtree=FALSE)
length(Heliamphora_tr$tip.label)
length(Heliamphora_tr2$tip.label)
is.ultrametric(Heliamphora_tr)
is.ultrametric(Heliamphora_tr2)
get_treeheight(Heliamphora_tr)
get_treeheight(Heliamphora_tr2)

# Branch length multiplier
calibration_age = 17.975096561  # We are just adding a bigger sister clade
                                # to Sarracenia, not fitting into a megatree clade
multiplier = calibration_age / get_treeheight(Heliamphora_tr2)

Heliamphora_tr3 = Heliamphora_tr2
Heliamphora_tr3$edge.length = Heliamphora_tr2$edge.length * multiplier


# Add to main tree
gbotb_tr6 = gbotb_tr5
length(gbotb_tr6$tip.label)

# Drop Heliamphora from main tree
TF = grepl(pattern="Heliamphora", x=gbotb_tr6$tip.label)
sum(TF)
tips_to_drop = gbotb_tr6$tip.label[TF]
gbotb_tr7 = drop.tip(phy=gbotb_tr6, tip=tips_to_drop, trim.internal=TRUE, subtree=FALSE)

# Add as sister to Sarracenia

# Get height of Sarracenia
TF = grepl(pattern="Sarracenia", x=gbotb_tr7$tip.label)
tips_to_keep = gbotb_tr7$tip.label[TF]
node_to_keep = getMRCA(phy=gbotb_tr7, tip=tips_to_keep)
Sarraceniaceae_gbotb_tr7 = extract.clade(phy=gbotb_tr7, node=node_to_keep)
height_of_Sarracenia = get_treeheight(Sarraceniaceae_gbotb_tr7)
height_of_Sarracenia

# Get height of Sarracenia + Heliamphora
TF1 = grepl(pattern="Sarracenia", x=gbotb_tr6$tip.label)
TF2 = grepl(pattern="Heliamphora", x=gbotb_tr6$tip.label)
gbotb_tr6$tip.label[TF2]
TF = (TF1 + TF2) > 0
tips_to_keep = gbotb_tr6$tip.label[TF]
node_to_keep = getMRCA(phy=gbotb_tr6, tip=tips_to_keep)
both_gbotb_tr6 = extract.clade(phy=gbotb_tr6, node=node_to_keep)
height_of_both = get_treeheight(both_gbotb_tr6)
height_of_both



plot(Sarraceniaceae_gbotb_tr7)
axisPhylo()

crown_age_Heliamphora = get_treeheight(Heliamphora_tr3)
crown_age_Sarracenia = get_treeheight(Sarraceniaceae_gbotb_tr7)
crown_age_both = get_treeheight(both_gbotb_tr6)
crown_age_Heliamphora
crown_age_Sarracenia
crown_age_both

# Add a root to Heliamphora
Heliamphora_tr3$root.edge = crown_age_both - crown_age_Heliamphora
Heliamphora_tr3$root.edge

# Distance below Sarracenia crown to side branch
distance_below_Sarracenia_crown = crown_age_both - crown_age_Sarracenia
distance_below_Sarracenia_crown

# New Heliamphora clade will be sister to Sarracenia at: node_to_keep=95286
TF = grepl(pattern="Sarracenia", x=gbotb_tr7$tip.label)
tips_to_keep = gbotb_tr7$tip.label[TF]
node_to_keep = getMRCA(phy=gbotb_tr7, tip=tips_to_keep)
Sarraceniaceae_gbotb_tr7 = extract.clade(phy=gbotb_tr7, node=node_to_keep)
gbotb_tr8 = bind.tree(x=gbotb_tr7, y=Heliamphora_tr3, where=node_to_keep, position=distance_below_Sarracenia_crown)

length(gbotb_tr6$tip.label)
length(gbotb_tr7$tip.label)
length(gbotb_tr8$tip.label)

is.ultrametric(gbotb_tr6)
is.ultrametric(gbotb_tr7)
is.ultrametric(gbotb_tr8)




#######################################################
# Double-check addition of Heliamphora
#######################################################
TF1 = grepl(pattern="Roridula", x=gbotb_tr8$tip.label)
TF2 = grepl(pattern="Actinidia", x=gbotb_tr8$tip.label)
TF3 = grepl(pattern="Clematoclethra", x=gbotb_tr8$tip.label)
#TF4 = grepl(pattern="Cyrilla", x=gbotb_tr8$tip.label)
TF4 = grepl(pattern="Saurauia", x=gbotb_tr8$tip.label)
TF5 = grepl(pattern="Sarracenia", x=gbotb_tr8$tip.label)
TF6 = grepl(pattern="Heliamphora", x=gbotb_tr8$tip.label)
TF7 = grepl(pattern="Darlingtonia", x=gbotb_tr8$tip.label)

TF = (TF1 + TF2 + TF3 + TF4 + TF5 + TF6 + TF7) > 0
sum(TF)
tips_to_keep = gbotb_tr8$tip.label[TF]
node_to_keep = getMRCA(phy=gbotb_tr8, tip=tips_to_keep)
Sarraceniaceae_gbotb_tr8 = extract.clade(phy=gbotb_tr8, node=node_to_keep)
plot(Sarraceniaceae_gbotb_tr8)
axisPhylo()





#######################################################
# What's left?  Lentibulariaceae
#######################################################

# Time-scale (very approximately) the Pinguicula2 cladogram via the dated Pinguicula1
# Requirements:
# * source() the file below 
# * have r8s installed somewhere like /Applications/r8s
#
source("/GitHub/bioinfRhints/R/r8s_functions_v2.R")


# Get height of Pinguicula from megatree
# The megatree topology seems to match Pinguicula2 better than Pinguicula1
TF = grepl(pattern="Pinguicula", x=gbotb_tr8$tip.label)
tips_to_keep = gbotb_tr8$tip.label[TF]
node_to_keep = getMRCA(phy=gbotb_tr8, tip=tips_to_keep)
Pinguicula_gbotb_tr8 = extract.clade(phy=gbotb_tr8, node=node_to_keep)
height_of_Pinguicula = get_treeheight(Pinguicula_gbotb_tr8)
height_of_Pinguicula

Pinguicula_gbotb_tr8$tip.label[Pinguicula_gbotb_tr8$tip.label %in% Pinguicula2_tr$tip.label == FALSE]

Pinguicula_gbotb_tr8_drop1 = drop.tip(Pinguicula_gbotb_tr8, tip="Pinguicula_caerulea")
Pinguicula_gbotb_tr8_drop2 = drop.tip(Pinguicula_gbotb_tr8_drop1, tip="Pinguicula_caussensis")

trtable = prt(Pinguicula_gbotb_tr8_drop2)
internal_nodes = (length(Pinguicula_gbotb_tr8_drop2$tip.label)+1):(length(Pinguicula_gbotb_tr8_drop2$tip.label)+Pinguicula_gbotb_tr8_drop2$Nnode)
internal_nodes
trtable$tipnames[internal_nodes]
trtable$time_bp[internal_nodes]



calibration_node_tip_specifiers = c("Pinguicula_alpina", "Pinguicula_vulgaris")
calibration_age = 15.839350182
r8s_nexus_fn = "Pinguicula2_r8s_nexus_fn.nex"
r8s_logfile_fn = "Pinguicula2_r8s_nexus_fn.log"
run_r8s_1calib(tr=Pinguicula2_tr, calibration_node_tip_specifiers=calibration_node_tip_specifiers, r8s_method="LF", addl_cmd="", calibration_age=calibration_age, nsites=1000, tmpwd=getwd(), r8s_nexus_fn=r8s_nexus_fn, r8s_logfile_fn=r8s_logfile_fn, r8s_path="/Applications/r8s")

file.remove(r8s_logfile_fn)


# Manually edit r8s NEXUS file for 11 calibrations
r8s_nexus_fn = "Pinguicula2_r8s_nexus_fn_v2.nex"
r8s_logfile_fn = "Pinguicula2_r8s_nexus_fn_v2.nex.log"

Pinguicula2_tr_dated = extract_tree_from_r8slog(logfn=r8s_logfile_fn, delimiter=" = ", printall=TRUE)
Pinguicula2_tr_rates = extract_rates_from_r8slog(logfn=r8s_logfile_fn)

plot(Pinguicula2_tr_dated)
axisPhylo()

# Manually edit this NEXUS file to include all timepoints
# Drop this:
# "Pinguicula_caerulea"


#######################################################
# Get the time-scaled Pinguicula tree, add to the main tree
#######################################################
# Just do linear scaling, as it's ultrametric

# Subset to remove Nepenthes
TF = grepl(pattern="Pinguicula", x=gbotb_tr8$tip.label)
sum(TF)
tips_to_drop = gbotb_tr8$tip.label[TF]
gbotb_tr9 = drop.tip(phy=gbotb_tr8, tip=tips_to_drop, trim.internal=TRUE, subtree=TRUE)
length(gbotb_tr8$tip.label)
length(gbotb_tr9$tip.label)
is.ultrametric(gbotb_tr8)
is.ultrametric(gbotb_tr9)

# Branch length multiplier
calibration_age = 1
multiplier = calibration_age / 1

gbotb_tr10 = gbotb_tr9
gbotb_tr10$edge.length = gbotb_tr9$edge.length * multiplier

# Add to main tree
TF = grepl(pattern="\\[", x=gbotb_tr9$tip.label)
tipnum = (1:length(gbotb_tr9$tip.label))[TF]
newtip = gbotb_tr9$tip.label[TF]
newtip

# New tip is "[14_tips]"
gbotb_tr10 = bind.tree(x=gbotb_tr9, y=Pinguicula2_tr_dated, where=tipnum)

length(gbotb_tr8$tip.label)
length(gbotb_tr9$tip.label)
length(gbotb_tr10$tip.label)

is.ultrametric(gbotb_tr8)
is.ultrametric(gbotb_tr9)
is.ultrametric(gbotb_tr10)



#######################################################
# Do Genlisea & Utricularia
#######################################################

plot(Genlisea_tr)
plot(Utricularia_tr)

Genlisea_tr1 = drop.tip(Genlisea_tr, "Utricularia_multifida")
Utricularia_tr1 = drop.tip(Utricularia_tr, "Outgroup")

plot(Genlisea_tr1)
plot(Utricularia_tr1)





#######################################################
# Time-scale the Genlisea tree, and add to the main tree
#######################################################
# Requirements:
# * source() the file below 
# * have r8s installed somewhere like /Applications/r8s
#

calibration_node_tip_specifiers = c("Genlisea_lobata", "Genlisea_roraimensis")
calibration_age = 26.396475
r8s_nexus_fn = "Genlisea_r8s_nexus_fn.nex"
r8s_logfile_fn = "Genlisea_r8s_nexus_fn.log"
run_r8s_1calib(tr=Genlisea_tr1, calibration_node_tip_specifiers=calibration_node_tip_specifiers, r8s_method="LF", addl_cmd="", calibration_age=calibration_age, nsites=1000, tmpwd=getwd(), r8s_nexus_fn=r8s_nexus_fn, r8s_logfile_fn=r8s_logfile_fn, r8s_path="/Applications/r8s")

Genlisea_tr_dated = extract_tree_from_r8slog(logfn=r8s_logfile_fn, delimiter=" = ", printall=TRUE)
Genlisea_tr_rates = extract_rates_from_r8slog(logfn=r8s_logfile_fn)

pdffn = "Genlisea_tr_dated.pdf"
pdf(file=pdffn, width=12, height=12)

plot(Genlisea_tr_dated, cex=0.5)
axisPhylo()

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


get_treeheight(Genlisea_tr_dated)

# Add to main tree
TF = grepl(pattern="Genlisea", x=gbotb_tr10$tip.label)
sum(TF)
tips_to_drop = gbotb_tr10$tip.label[TF]
gbotb_tr11 = drop.tip(phy=gbotb_tr10, tip=tips_to_drop, trim.internal=FALSE, subtree=TRUE)
length(gbotb_tr11$tip.label)
TF = grepl(pattern="\\[", x=gbotb_tr11$tip.label)
tipnum = (1:length(gbotb_tr11$tip.label))[TF]
newtip = gbotb_tr11$tip.label[TF]
newtip

# New tip is "[15_tips]"
gbotb_tr12 = bind.tree(x=gbotb_tr11, y=Genlisea_tr_dated, where=tipnum)
length(gbotb_tr10$tip.label)
length(gbotb_tr11$tip.label)
length(gbotb_tr12$tip.label)

is.ultrametric(gbotb_tr10)
is.ultrametric(gbotb_tr11)
is.ultrametric(gbotb_tr12)




#######################################################
# Add the Utricularia tree to the main tree
#######################################################
plot(Utricularia_tr1)
axisPhylo()

# Branch length multiplier
calibration_age = 29.853656
multiplier = calibration_age / get_tr_height(Utricularia_tr1)
Utricularia_tr2 = Utricularia_tr1
Utricularia_tr2$edge.length = Utricularia_tr1$edge.length * multiplier
get_tr_height(Utricularia_tr2)

# Add to main tree
TF = grepl(pattern="Utricularia", x=gbotb_tr12$tip.label)
sum(TF)
tips_to_drop = gbotb_tr12$tip.label[TF]
gbotb_tr13 = drop.tip(phy=gbotb_tr12, tip=tips_to_drop, trim.internal=FALSE, subtree=TRUE)
length(gbotb_tr13$tip.label)
TF = grepl(pattern="\\[", x=gbotb_tr13$tip.label)
tipnum = (1:length(gbotb_tr13$tip.label))[TF]
newtip = gbotb_tr13$tip.label[TF]
newtip

# New tip is [57_tips]
gbotb_tr14 = bind.tree(x=gbotb_tr13, y=Utricularia_tr2, where=tipnum)
length(gbotb_tr12$tip.label)
length(gbotb_tr13$tip.label)
length(gbotb_tr14$tip.label)

is.ultrametric(gbotb_tr12)
is.ultrametric(gbotb_tr13)
is.ultrametric(gbotb_tr14)


#######################################################
# Remarkable lack of overlap!
#######################################################
sort(tips_to_drop)
sort(Utricularia_tr1$tip.label)

length(tips_to_drop)
length(Utricularia_tr1$tip.label)

length(unique(c(tips_to_drop, Utricularia_tr1$tip.label)))

57 + 58
# 115


out_trfn = "gbotb_tr14.newick"
write.tree(gbotb_tr14, file=out_trfn)
out_trfn = "gbotb_wCPs_merged_v14.newick"
write.tree(gbotb_tr14, file=out_trfn)

length(gbotb_tr14$tip.label)
# 72728 tips

#######################################################
# NOTES
#######################################################

# These Pinguicula were dropped in the switch from the big tree to 
# the cladogram tree:
Pinguicula_caerulea
Pinguicula_caussensis


gbotb_tr14$tip.label[grepl(pattern="Byblis", x=gbotb_tr14$tip.label)]
# Need characters for there:
Byblis_gigantea
Byblis_lamellata
Byblis_liniflora

gbotb_tr14$tip.label[grepl(pattern="Roridula", x=gbotb_tr14$tip.label)]
Roridula_gorgonias
Roridula_dentata

The original big tree had 57 Utricularia
The digitized dated tree had 58 Utricularia

However, only 4 of these overlapped!  Merging the two trees would create 
111 Utricularia.

Given that there are hundreds of Utricularia species, perhaps this is not 
surprising; but, look for additional Utricularia phylogenies, or consider
making your own from the genetic data in the future.



Genlisea_bigtree = gbotb_tr$tip.label[grepl(pattern="Genlisea", x=gbotb_tr$tip.label)]
Genlisea_molecular = gbotb_tr14$tip.label[grepl(pattern="Genlisea", x=gbotb_tr14$tip.label)]

Genlisea_bigtree[Genlisea_bigtree %in% Genlisea_molecular == FALSE]

# This was dropped in the switch to the new tree:
Genlisea_pygmaea

# Carnivorous bromeliads:
sum(grepl(pattern="Brocchinia_hechtioides", x=gbotb_tr14$tip.label))
sum(grepl(pattern="Brocchinia_reducta", x=gbotb_tr14$tip.label))
sum(grepl(pattern="Catopsis_berteroniana", x=gbotb_tr14$tip.label))




#######################################################
# Trim the huge tree to just have the 2 or 3 sister 
# groups to each carnivorous clade
#######################################################
# Droseraceae / Nepentheaceae clade:
outgroup_tips = NULL


tip_genera = c("Drosera",
"Dionaea",
"Aldrovanda",
"Drosophyllum",
"Nepenthes",
"Triphyophyllum")

TF = rep(FALSE, times=length(gbotb_tr14$tip.label))
for (i in 1:length(tip_genera))
	{
	tmpTF = grepl(pattern=tip_genera[i], x=gbotb_tr14$tip.label)
	TF = TF + tmpTF
	}
TF = TF > 0
tips_to_ID_node = gbotb_tr14$tip.label[TF]
CPnode = getMRCA(phy=gbotb_tr14, tip=tips_to_ID_node)
length(tips_to_ID_node)
CPnode

parent_node1 = get_parent(nodenum=CPnode, t=gbotb_tr14)
parent_node2 = get_parent(nodenum=parent_node1, t=gbotb_tr14)
parent_node3 = get_parent(nodenum=parent_node2, t=gbotb_tr14)

parent_node1
parent_node2
parent_node3

outgroup_tips1 = get_daughter_tipnums(nodenum=parent_node1, tr=gbotb_tr14)
outgroup_tips2 = get_daughter_tipnums(nodenum=parent_node2, tr=gbotb_tr14)
outgroup_tips3 = get_daughter_tipnums(nodenum=parent_node3, tr=gbotb_tr14)

length(outgroup_tips1) # 775 already
length(outgroup_tips2)
length(outgroup_tips3)
outgroup_tips = c(outgroup_tips, outgroup_tips3)



tip_genera = c("Roridula",
"Actinidia",
"Clematoclethra",
"Saurauia",
"Sarracenia",
"Heliamphora",
"Darlingtonia")

TF = rep(FALSE, times=length(gbotb_tr14$tip.label))
for (i in 1:length(tip_genera))
	{
	tmpTF = grepl(pattern=tip_genera[i], x=gbotb_tr14$tip.label)
	TF = TF + tmpTF
	}
TF = TF > 0
tips_to_ID_node = gbotb_tr14$tip.label[TF]
CPnode = getMRCA(phy=gbotb_tr14, tip=tips_to_ID_node)
length(tips_to_ID_node)
CPnode

parent_node1 = get_parent(nodenum=CPnode, t=gbotb_tr14)
parent_node2 = get_parent(nodenum=parent_node1, t=gbotb_tr14)
parent_node3 = get_parent(nodenum=parent_node2, t=gbotb_tr14)

parent_node1
parent_node2
parent_node3

outgroup_tips1 = get_daughter_tipnums(nodenum=parent_node1, tr=gbotb_tr14)
outgroup_tips2 = get_daughter_tipnums(nodenum=parent_node2, tr=gbotb_tr14)
outgroup_tips3 = get_daughter_tipnums(nodenum=parent_node3, tr=gbotb_tr14)

length(outgroup_tips1) # 1606
length(outgroup_tips2)
length(outgroup_tips3)
outgroup_tips = c(outgroup_tips, outgroup_tips3)






tip_genera = c("Pinguicula",
"Genlisea",
"Utricularia")

TF = rep(FALSE, times=length(gbotb_tr14$tip.label))
for (i in 1:length(tip_genera))
	{
	tmpTF = grepl(pattern=tip_genera[i], x=gbotb_tr14$tip.label)
	TF = TF + tmpTF
	}
TF = TF > 0
tips_to_ID_node = gbotb_tr14$tip.label[TF]
CPnode = getMRCA(phy=gbotb_tr14, tip=tips_to_ID_node)
length(tips_to_ID_node)
CPnode

parent_node1 = get_parent(nodenum=CPnode, t=gbotb_tr14)
parent_node2 = get_parent(nodenum=parent_node1, t=gbotb_tr14)
parent_node3 = get_parent(nodenum=parent_node2, t=gbotb_tr14)

parent_node1
parent_node2
parent_node3

outgroup_tips1 = get_daughter_tipnums(nodenum=parent_node1, tr=gbotb_tr14)
outgroup_tips2 = get_daughter_tipnums(nodenum=parent_node2, tr=gbotb_tr14)
outgroup_tips3 = get_daughter_tipnums(nodenum=parent_node3, tr=gbotb_tr14)

length(outgroup_tips1) # 550
length(outgroup_tips2) # 560
length(outgroup_tips3) # 706
outgroup_tips = c(outgroup_tips, outgroup_tips3)



tip_genera = c("Cephalotus")

TF = rep(FALSE, times=length(gbotb_tr14$tip.label))
for (i in 1:length(tip_genera))
	{
	tmpTF = grepl(pattern=tip_genera[i], x=gbotb_tr14$tip.label)
	TF = TF + tmpTF
	}
TF = TF > 0
tips_to_ID_node = gbotb_tr14$tip.label[TF]
#CPnode = getMRCA(phy=gbotb_tr14, tip=tips_to_ID_node)

tipnode_TF = grepl(pattern="Cephalotus", x=gbotb_tr14$tip.label)
nodenum = (1:length(gbotb_tr14$tip.label))[tipnode_TF]
nodenum

parent_node1 = get_parent(nodenum=nodenum, t=gbotb_tr14)
parent_node2 = get_parent(nodenum=parent_node1, t=gbotb_tr14)
parent_node3 = get_parent(nodenum=parent_node2, t=gbotb_tr14)

parent_node1
parent_node2
parent_node3

outgroup_tips1 = get_daughter_tipnums(nodenum=parent_node1, tr=gbotb_tr14)
outgroup_tips2 = get_daughter_tipnums(nodenum=parent_node2, tr=gbotb_tr14)
outgroup_tips3 = get_daughter_tipnums(nodenum=parent_node3, tr=gbotb_tr14)

length(outgroup_tips1) # 6
length(outgroup_tips2) # 172
length(outgroup_tips3) # 273
outgroup_tips = c(outgroup_tips, outgroup_tips3)




tip_genera = c("Byblis")

TF = rep(FALSE, times=length(gbotb_tr14$tip.label))
for (i in 1:length(tip_genera))
	{
	tmpTF = grepl(pattern=tip_genera[i], x=gbotb_tr14$tip.label)
	TF = TF + tmpTF
	}
TF = TF > 0
tips_to_ID_node = gbotb_tr14$tip.label[TF]
CPnode = getMRCA(phy=gbotb_tr14, tip=tips_to_ID_node)
length(tips_to_ID_node)
CPnode

parent_node1 = get_parent(nodenum=CPnode, t=gbotb_tr14)
parent_node2 = get_parent(nodenum=parent_node1, t=gbotb_tr14)
parent_node3 = get_parent(nodenum=parent_node2, t=gbotb_tr14)

parent_node1
parent_node2
parent_node3

outgroup_tips1 = get_daughter_tipnums(nodenum=parent_node1, tr=gbotb_tr14)
outgroup_tips2 = get_daughter_tipnums(nodenum=parent_node2, tr=gbotb_tr14)
outgroup_tips3 = get_daughter_tipnums(nodenum=parent_node3, tr=gbotb_tr14)

length(outgroup_tips1) # 4
length(outgroup_tips2) # 2777
length(outgroup_tips3) # 3065
outgroup_tips = c(outgroup_tips, outgroup_tips3)




tip_genera = c("Brocchinia")

TF = rep(FALSE, times=length(gbotb_tr14$tip.label))
for (i in 1:length(tip_genera))
	{
	tmpTF = grepl(pattern=tip_genera[i], x=gbotb_tr14$tip.label)
	TF = TF + tmpTF
	}
TF = TF > 0
tips_to_ID_node = gbotb_tr14$tip.label[TF]
CPnode = getMRCA(phy=gbotb_tr14, tip=tips_to_ID_node)
length(tips_to_ID_node)
CPnode

parent_node1 = get_parent(nodenum=CPnode, t=gbotb_tr14)
parent_node2 = get_parent(nodenum=parent_node1, t=gbotb_tr14)
parent_node3 = get_parent(nodenum=parent_node2, t=gbotb_tr14)

parent_node1
parent_node2
parent_node3

outgroup_tips1 = get_daughter_tipnums(nodenum=parent_node1, tr=gbotb_tr14)
outgroup_tips2 = get_daughter_tipnums(nodenum=parent_node2, tr=gbotb_tr14)
outgroup_tips3 = get_daughter_tipnums(nodenum=parent_node3, tr=gbotb_tr14)

length(outgroup_tips1) # 789
length(outgroup_tips2) # 813
length(outgroup_tips3) # 4297
outgroup_tips = c(outgroup_tips, outgroup_tips3)






tip_genera = c("Catopsis")

TF = rep(FALSE, times=length(gbotb_tr14$tip.label))
for (i in 1:length(tip_genera))
	{
	tmpTF = grepl(pattern=tip_genera[i], x=gbotb_tr14$tip.label)
	TF = TF + tmpTF
	}
TF = TF > 0
tips_to_ID_node = gbotb_tr14$tip.label[TF]
CPnode = getMRCA(phy=gbotb_tr14, tip=tips_to_ID_node)
length(tips_to_ID_node)
CPnode

parent_node1 = get_parent(nodenum=CPnode, t=gbotb_tr14)
parent_node2 = get_parent(nodenum=parent_node1, t=gbotb_tr14)
parent_node3 = get_parent(nodenum=parent_node2, t=gbotb_tr14)

parent_node1
parent_node2
parent_node3

outgroup_tips1 = get_daughter_tipnums(nodenum=parent_node1, tr=gbotb_tr14)
outgroup_tips2 = get_daughter_tipnums(nodenum=parent_node2, tr=gbotb_tr14)
outgroup_tips3 = get_daughter_tipnums(nodenum=parent_node3, tr=gbotb_tr14)

length(outgroup_tips1) # 789
length(outgroup_tips2) # 813
length(outgroup_tips3) # 4297
outgroup_tips = c(outgroup_tips, outgroup_tips3)


uniq_outgroup_tips = unique(outgroup_tips)
length(outgroup_tips)
length(uniq_outgroup_tips)






#######################################################
# For each non-carnivorous genus, remove all but 1
#######################################################

CP_genera = c("Drosera",
"Dionaea",
"Aldrovanda",
"Drosophyllum",
"Nepenthes",
"Triphyophyllum",
"Habropetalum",
"Ancistrocladus",
"Roridula",
"Actinidia",
"Clematoclethra",
"Saurauia",
"Sarracenia",
"Heliamphora",
"Darlingtonia",
"Pinguicula",
"Genlisea",
"Utricularia",
"Cephalotus",
"Byblis",
"Brocchinia",
"Catopsis")

uniq_outgroup_tipnames = gbotb_tr14$tip.label[uniq_outgroup_tips]
uniq_outgroup_tipnames_reduced = uniq_outgroup_tipnames
for (i in 1:length(CP_genera))
	{
	TF = grepl(CP_genera[i], x=uniq_outgroup_tipnames_reduced)
	uniq_outgroup_tipnames_reduced = uniq_outgroup_tipnames_reduced[TF==FALSE]
	}
length(uniq_outgroup_tipnames)
length(uniq_outgroup_tipnames_reduced)


# Get the unique genera
list_of_genera = rep("", times=length(uniq_outgroup_tipnames_reduced))
for (i in 1:length(uniq_outgroup_tipnames_reduced))
	{
	genus = strsplit(uniq_outgroup_tipnames_reduced[i], split="_")[[1]][1]
	list_of_genera[i] = genus
	}

uniq_list_of_genera = unique(list_of_genera)
length(list_of_genera)
length(uniq_list_of_genera)


tips_to_keep = rep("", times=length(uniq_list_of_genera))
for (i in 1:length(uniq_list_of_genera))
	{
	genus_hits_TF = grepl(pattern=uniq_list_of_genera[i], x=uniq_outgroup_tipnames_reduced)
	genus_hits = uniq_outgroup_tipnames_reduced[genus_hits_TF]
	tips_to_keep[i] = genus_hits[1]
	}

tips_to_keep
length(tips_to_keep)

tips_to_drop_TF = (uniq_outgroup_tipnames_reduced %in% tips_to_keep) == FALSE
tips_to_drop = uniq_outgroup_tipnames_reduced[tips_to_drop_TF]


tips_to_drop2 = gbotb_tr14$tip.label

for (i in 1:length(CP_genera))
	{
	TF = grepl(CP_genera[i], x=tips_to_drop2)
	tips_to_drop2 = tips_to_drop2[TF==FALSE]
	}

tips_to_drop3 = tips_to_drop2[(tips_to_drop2 %in% tips_to_keep) == FALSE]

gbotb_tr14_wSister_genera = drop.tip(phy=gbotb_tr14, tip=tips_to_drop3)

length(gbotb_tr14_wSister_genera$tip.label)


#plot(gbotb_tr14_wSister_genera, tip.label=FALSE)
#axisPhylo()




