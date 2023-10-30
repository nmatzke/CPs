# devtools::install_github("jinyizju/V.PhyloMaker2")

library(ape)
library(BioGeoBEARS)
library(V.PhyloMaker2)

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

sdf = as.data.frame(cbind(genera, root_age, num_species), stringsAsFactors=FALSE)
sdf

node_ages


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
tree2[["Drosophyllaceae"]] = GBOTB.extended.WP

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

names2 = NULL
root_age2 = NULL
num_species2 = NULL
for (i in 1:length(genera2))
	{
	names2 = c(names2, names(genera2)[i])

	TF = rep(FALSE, times=length(tree2[[i]]$tip.label))
	for (j in 1:length(genera2[[i]]))
		{
		tmpTF = grepl(pattern=genera[[i]][j], x=tree2[[i]]$tip.label)
		TF = TF + tmpTF
		}
	num_species2 = c(num_species2, sum(TF))
	
	tips_to_keep = tree2[[i]]$tip.label[TF]
	node_to_keep = getMRCA(phy=tree2[[i]], tip=tips_to_keep)
	if (sum(TF) > 1)
		{
		subtree = extract.clade(phy=tree2[[i]], node_to_keep)
		subtr_age = get_max_height_tree(subtree)
		} else {
		subtr_age = 0.0
		}
	root_age2 = c(root_age2, subtr_age)
	}

sdf2 = as.data.frame(cbind(names2, root_age2, num_species2), stringsAsFactors=FALSE)
sdf2


get_treeheight(Droseraceae_tr)
get_treeheight(Genlisea_tr)
get_treeheight(Pinguicula1_tr)
get_treeheight(Pinguicula2_tr)
get_treeheight(Utricularia_tr)
get_treeheight(Nepenthes_tr)
get_treeheight(Heliamphora_tr)
get_treeheight(Sarracenia_tr)




 

