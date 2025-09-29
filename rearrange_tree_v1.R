
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
pdf(file=pdffn, width=20, height=35)

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
