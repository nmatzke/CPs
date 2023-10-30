# 2015: GraphClick now free online
# http://www.arizona-software.ch/graphclick/disclaimer.html

# 2023: GraphClick is dead, long live WebPlotDigitizer on Manual mode
# https://automeris.io/WebPlotDigitizer/download.html

#Example run command:

#===================================================
# Source the TreeRogue functions below (input the directory and file where you
# have saved the source):
#===================================================
library(BioGeoBEARS)
sourcedir = '/Users/masa/Downloads/setup_TreeRogue/'
source3 = '_genericR_v2.R'
source(paste(sourcedir, source3, sep=""))
#===================================================


#===================================================
# Run with these commands
#===================================================
library(ape)

# put the text files (at bottom) into your working directory

runjunk = '
wd = "/drives/Dropbox/_notes_install/setup_TreeRogue"
setwd(wd)

xy2 = treerogue_read_files()
xy3 = treerogue_associate_branch_bottoms_with_nodes(xy2)
tr = build_tree_using_corners(xy3)
plot(tr)

##From here undated
tr2B = tr
plot(tr2B, cex=0.7)

title("Nepenthaceae phylogeny from Murphy et al (2020), Fig. 4")
trfn = "Murphy_2020_Fig4_Nepenthaceae_digitized_v1.newick"
write.tree(tr2B, file=trfn)

pdffn = "tr2B.pdf"
pdf(file=pdffn, width=7.5, height=9)
dev.off()
cmdstr = paste0("open ", pdffn)


## From here dated trees
# Flatten the tips
prttr = prt(tr)
ntips = length(tr$tip.label)
tipdates = prttr$time_bp[1:ntips]
meantipdate = mean(tipdates)
meantipdate

tr2 = tr
how_to_change_tip_edge_length = tipdates - meantipdate

for (i in 1:ntips)
	{
	tipnum = i
	edgeTF = tr2$edge[,2] == tipnum
	edgenum = (1:nrow(tr2$edge))[edgeTF]
	edgenum
	tr2$edge.length[edgenum] = tr2$edge.length[edgenum] + how_to_change_tip_edge_length[tipnum]
	}
	
prt(tr2)
is.ultrametric(tr)
is.ultrametric(tr2)



# Scale the tree
xmin = 0.0
xmax = 87.8372100
digitized_distance = xmax - xmin
actual_distance = 88 	# million years
multiply_branches_by = actual_distance / digitized_distance
tr3 = tr2
tr3$edge.length = tr3$edge.length * multiply_branches_by
plot(tr3)
axisPhylo()

trfn = "Liu_Smith_2020_Fig1_Heliamphora_digitized_v1.newick"
write.tree(tr3, file=trfn)

# Make a PDF
pdffn = "Liu_Smith_2020_Fig1_Heliamphora_digitized_v1.pdf"
pdf(file=pdffn, width=7.5, height=9)

plot(tr3, cex=0.7)
axisPhylo()
title("Dated Heliamphora phylogeny from Liu & Smith (2020), Fig. 1")
mtext(text="Millions of years ago", side=1, line=2, cex=0.8)


dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)
# End of PDF

'


#===================================================
# Source
#===================================================


##############################################################
# TreeRogue 0.21: TreeThief-like tree grabbing using x,y 
# coordinates digitized from an image of a phylogenetic tree.
##############################################################
# GOAL: to process x, y coordinates into a Newick-format tree
##############################################################
# by Nick Matzke
# Copyright 2010-infinity
# matzke@berkeley.edu
# 10/27/2010
#
# Please link/cite if you use this, email me if you have 
#   thoughts/improvements/corrections.
#
##############################################################
#
# Free to use/redistribute under:
# Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0) 
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the above license, linked here:
# 
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# 
# Summary:
#
# You are free:
#
#   * to Share ? to copy, distribute and transmit the work
#   * to Remix ? to adapt the work
#
# Under the following conditions:
#
#   * Attribution ? You must attribute the work in the manner 
#     specified by the author or licensor (but not in any way that 
#     suggests that they endorse you or your use of the work).
#   * Noncommercial ? You may not use this work for commercial purposes. 
#
#   * Share Alike ? If you alter, transform, or build upon this work,
#     you may distribute the resulting work only under the same or
#     similar license to this one. 
#
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# 
###################################################################
#
# #===================================================
# Run with these commands (after sourcing)
# ===================================================
# library(ape)
# 
# put the text files (at bottom) into your working directory
# wd = "/Users/nick/Desktop/__projects/2010-11-01_cone_snails/"
# setwd(wd)
# 
# xy2 = treerogue_read_files()
# xy3 = treerogue_associate_branch_bottoms_with_nodes(xy2)
# tr = build_tree_using_corners(xy3)
# plot(tr)
#
####################################################################
#
# NOTES
#
# *Heavily* modified from a very limited script posted here by 
# bbolker@gmail.com:
# https://stat.ethz.ch/pipermail/r-sig-phylo/2008-November/000173.html
#
# Background: I worked up this script after getting frustrated at
# (1) the failure of the famous "TreeThief" to work on any modern machine;
# (2) my failure to get the newer "TreeSnatcher" to work on Mac OS X 10.4.11
# (some weirdness about Java and X11, as far as I can tell), and
# (3) no other good options.
#
# Summary: This script takes in x,y coordinates (with the lower left as the origin)
# of nodes, tips, and branch bottoms ("corners"), and builds a tree out of it.
# 
# It assumes:
# (1) Your tree is horizontal left-to-right, with the tips on the right
# (2) Your tree is a "square" tree (i.e. no diagonal/curved branches)
#
# I captured my x,y coordinates using GraphClick 3.0, available for $8 at:
# http://www.arizona-software.ch/graphclick/
#
# (there is a free trial download, but only lets you export 10 coordinates at a
# time, so it is pointless)
#
# REQUIRED INPUTS:
#   (for each, digitize the relevant points in GraphClick, and
#    File-->Export to a text file):
# 
# (Note: all text files should have a header line)
#
# 1. A tab-delimited text file with x and y for each internal node
#
# 2. A tab-delimited text file with x and y for each tip
#
# 2a. A text file with tip names, in order from top-to-bottom
# 
# 3. A tab-delimited text file with x and y for each tip for each "corner"
#    (i.e., the bottom of each branch).
#
# 4. For now, do NOT digitize the bottom of the root of the tree, if the
#    image you are digitizing has one.  You could add the root length later
#    manually, if you like.
#
# 5. The tree must be fully dichotomous (if the image you are digitizing is not,
#    you can "fake it" by resolving polytomies by clicking digitization points
#    to, in effect, manually resolve these polytomies with very short branches.
#    Note that you will have to add a corner for each internal node you add (!).
#
#    The R APE package can collapse short branches to polytomies later, if you like.
#
# Trees to not have to be ultrametric, and digitization does not have to be 
# exact -- the script will attempt to match the most likely branch bottoms
# to the nodes (a graphic is produced by R so that you can check the results
# and tinker if necessary).
#
# Requires the APE library.
# 
# COMMON PROBLEMS WHILE DIGITIZING
# * missing internal nodes or corners
# * double-clicking a node so there are 2 nodes where there should be 1
#
# 
#############################################################

library(ape)


# Assumes default filenames
treerogue_read_files <- function()
	{
	internal = read_table_good("internal.txt")
	tips = read_table_good("tips.txt")
	tipnames = read_table_good("tipnames.txt")
	corners = read_table_good("corners.txt")
	
	# sort the tips from top to bottom in y
	tips = tips[order(tips$y, decreasing = TRUE), ]
	
	# sort the internals from left to right in x
	internal = internal[order(internal$x, decreasing=FALSE), ]
	
	if (nrow(tips) != nrow(tipnames))
		{
		print("ERROR: the number of tips must equal the length of the tipnames!")
		print(paste("Instead, nrow(tipnames) =", nrow(tipnames), "and nrow(tips) =", nrow(tips), sep=" "))
		}
	
	if ((nrow(tips)-1) != nrow(internal))
		{
		print("ERROR: the number of tips-1 must equal the number of the internal nodes!")
		print(paste("Instead, nrow(tips) =", nrow(tips), "and nrow(internal) =", nrow(internal), sep=" "))
		}
	
	nodetypes = c(rep("tip", nrow(tipnames)), rep("internal", nrow(internal)))
	nodenames = unlist(c(tipnames, rep("", nrow(internal))))
	xy = rbind(tips, internal)
	xy2 = cbind(xy, nodetypes, nodenames)
	xy2 = as.data.frame(xy2)
	names(xy2) = c("x", "y", "nodetypes", "nodenames")
	
	xy2 = df_nonum_factors_to_char(xy2, max_NAs=0.5)
	
	xy2$nodetypes[xy2$x == min(xy2$x)] = "root"
	
	if (nrow(corners) != (nrow(xy2)-1))
		{
		print("ERROR: the number of corners must equal the number of nodes-1  !")
		print(paste("Instead, length(nodes) =", nrow(xy2), "and nrow(corners) =", nrow(corners), sep=" "))
		}

	# sort file so that the tips are first
	# tip nodes in order from top to bottom:
	xytips = xy[xy$tipname != "", ]
	
	tips_in_order = xy2$tipname[xy2$tipname != ""]

	return(xy2)

	}


treerogue_associate_branch_bottoms_with_nodes <- function(xy2)
	{
	# Load the coordinates of the corners
	corners = read_table_good("corners.txt")
	bots = corners
	

	# Get xy data
	xx = xy2$x
	yy = xy2$y
	
	df = associate_branch_bottoms_with_nodes(xx, yy, bots)

	xy3 = cbind(xy2, df$chosen_node_bottoms_xx, df$chosen_node_bottoms_yy)
	names(xy3) = c(names(xy2), "cx", "cy")
	
	write.table(xy3, "linked_corners.txt", quote=FALSE, sep="	", row.names = FALSE,
            col.names = TRUE)

	return(xy3)
	}

treerogue_associate_branch_bottoms_with_nodes2 <- function(xy1)
	{
	# Load the coordinates of the corners
	corners = read_table_good("corners.txt")
	bots = corners
	

	# Get xy data
	xx = xy1$x
	yy = xy1$y
	
	xy22 = associate_branch_bottoms_with_nodes2(xx, yy, bots)
	
	b_lab = xy22$chosen_node_bottoms_labels
	t_lab = paste("t", xy22$topnodes, sep="")
	xy23 = cbind(xy22[,1:4], xy22[,6:7], b_lab, t_lab)
	names(xy23) = c("num", "tx", "ty", "nobot", "bx", "by", "b_lab", "t_lab")
	
	write.table(xy23, "linked_corners.txt", quote=FALSE, sep="	", row.names = FALSE,
            col.names = TRUE)

	return(xy23)
	}


# Associate branch bottom coordinates with nodes, and plot the results;
# The user may then edit the output associations if they so desire.
associate_branch_bottoms_with_nodes2 <- function(xx, yy, bots)
	{
	# There should be one less branch bottom than there are internal nodes
	# (because the root node (should) have no digitized "corner" beneath it)
	nodes = 1:length(xx)
	if (length(nodes) != nrow(bots) +1)
		{
		print("ERROR: the number of corners must equal the number of nodes-1  !")
		print(paste("Instead, length(nodes) =", length(nodes), "and nrow(bots) =", nrow(bots), sep=" "))
		} else {
		# OK, find bottom of branch to go with each top of branch
		# an array saying which branches have a bottom
		node_with_no_bottom = rep(TRUE, length(nodes))

		# these are the remaining branch bottoms that have not been associated yet
		bots_with_no_top = rep(TRUE, nrow(bots))
		bxx = bots$x
		byy = bots$y

		# Empty values to hold nodes
		chosen_node_bottoms_xx = rep(NA, length(xx))
		chosen_node_bottoms_yy = rep(NA, length(yy))
		chosen_node_bottoms_labels = rep(NA, length(yy))
		chosen_node_tops_labels = rep(NA, length(byy))

		
		botsnodes = 1:nrow(bots)
		botvals = cbind(botsnodes, bxx, byy, bots_with_no_top, chosen_node_tops_labels)
		botvals = adf(botvals)
		botvals$bots_with_no_top = TRUE
		
		topnodes = 1:length(xx)
		topvals = cbind(topnodes, xx, yy, node_with_no_bottom, chosen_node_bottoms_labels)
		topvals = adf(topvals)
		topvals$node_with_no_bottom = TRUE
		
		# remove minx row
		#minrow = topvals[topvals$xx == min(topvals$xx, na.rm=TRUE), ]
		#minrownum = as.numeric(minrow[1])
		#print(nrow(topvals))
		#topvals = topvals[-minrownum, ]
		#print(nrow(topvals))		
		
		
		i=1
		#while (sum(node_with_no_bottom) > 1)
		num_iterations = length(xx)
		for (i in 1:num_iterations)
			{
			#i=i+1
			#print(i)
			# look for branch bottom coordinates in a narrow window to the left of the node
			# basically (a) smallest slope and (b) closest distance in x
			
			## find next node to include (the rightmost internal node)
			available_rows = topvals[topvals$node_with_no_bottom==TRUE, ]
			maxrow = which(available_rows$xx == max(available_rows$xx, na.rm=TRUE))[1]
			
			#nextnode <- which( maxvals == max(maxvals, na.rm=TRUE) )[1]
			#nextnode = i
			nextnode = available_rows$topnodes[maxrow]
			
			####################################
			# Find the best matching branch bottom/corner for each node
			####################################
			# This is trial-and-error, you may have to plink to find a function 
			# that works.
			# That, or do manual edits to the tree later...
			####################################
			ydist <- (botvals$byy[bots_with_no_top==TRUE]) - topvals$yy[nextnode]
			xdist <- (botvals$bxx[bots_with_no_top==TRUE]) - topvals$xx[nextnode]
			#xdist <- (bxx*bots_with_no_top) - xx[nextnode]
			#tmp_botsnodes <- botsnodes*bots_with_no_top


			
			# Rank of the y distances
			rank_ydist = rank(abs(ydist))

			# calculate the slops
			xyslopes <- abs(ydist/xdist)

			
			# the best ancestor will have a low slope to the branch bottom, and a short negative distance in x
			xdist_neg = xdist

			xdist_neg[xdist > 0] = 0

			#print("hi0")
			#print(xdist)
			xdist_neg[xdist < 0] = -1 * xdist_neg[xdist < 0]

			# normalize to units of minimum absolute distance
			#min_dist = (min(abs(xdist[xdist!=0]), na.rm=TRUE))
			min_dist = 0.001
			xdist_neg_norm = (xdist_neg / min_dist)


			# short positive distances are less good (half as good) than short negative distances
			xdist_pos = xdist
			xdist_pos[xdist < 0] = 0
			xdist_pos[xdist > 0] = xdist_pos[xdist > 0]
			#xdist_pos_norm = (xdist_pos / min_dist) * 100
			xdist_pos_norm = (xdist_pos / min_dist) * 100

			rank_xdist = rank_ydist
			rank_xdist[xdist <= 0] = 1
			rank_xdist[xdist > 0] = 2


			
			###########################
			# Plink here especially...
			###########################
			rank_slope = (xyslopes^2)
			#final_rank = rank_ydist * abs(ydist) + 1*xyslopes^0.5 * xdist_neg_norm + xdist_pos_norm
			final_rank = (rank_ydist * abs(ydist)) *  xdist_neg_norm * xdist_pos_norm
			
			#xyslopes *
			###########################

			branch_bot_fits = final_rank			
			best_fit = which(branch_bot_fits == min(branch_bot_fits, na.rm=TRUE))[1]
			
			best_fit_row = botvals[bots_with_no_top==TRUE, ][best_fit, ]
			best_fit_rownum = as.numeric(best_fit_row[1])
			
			#print(best_fit_row)
			#print(best_fit_rownum)
			
			#bottom_to_add = botsnodes[bots_with_no_top][best_fit]
			#bottom_to_add = tmp_botsnodes[best_fit]
			#bottom_to_add_rownum = botvals$botsnodes[best_fit_rownum]
			topvals$chosen_node_bottoms_labels[nextnode] = paste("b", best_fit_rownum, sep="")
			
			#print("hi1")
			#print(paste("n", nextnode, sep=""))
			#botvals$chosen_node_tops_labels[best_fit_rownum] = paste("n", nextnode, sep="")
			#print("hi2")
			
			chosen_node_bottoms_xx[nextnode] = botvals$bxx[best_fit_rownum]
			chosen_node_bottoms_yy[nextnode] = botvals$byy[best_fit_rownum]
			
			
			#xx = xx[-bottom_to_add]
			#yy = yy[-bottom_to_add]
			#bxx = bxx[-bottom_to_add]
			#byy = byy[-bottom_to_add]
			#tmp_botsnodes = tmp_botsnodes[tmp_botsnodes != bottom_to_add]
			
			# remove the node from the list needing branch bottoms
			topvals$node_with_no_bottom[nextnode] = FALSE
			botvals$bots_with_no_top[best_fit_rownum] = FALSE

			
			}
		}
	
	#tmb_bot_labels = paste("b", 1:nrow(bots), sep="")
	
	#NA_index = which(is.na(chosen_node_bottoms_yy))
	#tmb_bot_labels2 = c(tmb_bot_labels[1:NA_index], NA, tmb_bot_labels[(NA_index+1) : length(tmb_bot_labels)] )
	
	
	#b_lab = t(t(tmb_bot_labels2))
	#b_lab = chosen_node_bottoms_labels
	#t_lab = t(t(paste("n", 1:length(xx), sep="")))
	#t_lab = chosen_node_tops_labels

	#tmpdata = cbind(nodes, zxx, zyy, chosen_node_bottoms_xx, chosen_node_bottoms_yy)
	tmpdata = cbind(topvals, chosen_node_bottoms_xx, chosen_node_bottoms_yy)
	df2 = adf(tmpdata)
	names(df2) = c("topnodes", "xx", "yy", "node_with_no_bottom", "chosen_node_bottoms_labels", "chosen_node_bottoms_xx", "chosen_node_bottoms_yy")
	#print(tmpdata)
	
	#df = as.data.frame(tmpdata)
	#startnames = length(names(df))-1
	#endnames = length(names(df))-0
	#names(df)[startnames:endnames] = c("b_lab", "t_lab")
	plot(df2$yy, df2$chosen_node_bottoms_yy)

	
	plot(c(df2$xx, bots$x), c(df2$yy, bots$y), pch="")
	points(df2$xx, df2$yy, pch="n")
	points(df2$chosen_node_bottoms_xx, df2$chosen_node_bottoms_yy, pch="b")
	title("Use this plot to check if branch bottoms match nodes")
	segments(df2$xx, df2$yy, df2$chosen_node_bottoms_xx, df2$chosen_node_bottoms_yy)
	
	plot(c(df2$xx, bots$x), c(df2$yy, bots$y), pch="")
	text(df2$xx, df2$yy, label=paste("t", df2$topnodes, sep=""), cex=0.6, col="blue", pos=4)
	text(df2$chosen_node_bottoms_xx, df2$chosen_node_bottoms_yy, label=df2$chosen_node_bottoms_labels, cex=0.6, col="red", pos=2)
	title("This plot has node numbers.")
	segments(df2$xx, df2$yy, df2$chosen_node_bottoms_xx, df2$chosen_node_bottoms_yy, cex=0.2)
	

	return(df2)
	}



build_tree_using_corners2 <- function(xy23)
	{
	# define the tip.labels
	tip.labels = xy23$nodenames
	tip.labels = tip.labels[tip.labels != ""]
	if (!missing(tip.labels))
		{
		ntips <- length(tip.labels)
		}
	
	
	xx = xy23$tx
	yy = xy23$ty
	cx = xy23$bx
	cy = xy23$by
	
	nodes <- 1:length(xx)
	is.tip <- nodes <= ntips

	# keep track of the nodes which are unlinked
	unlinked_nodes = rep(TRUE, length(nodes))

	
	# Checks (kinda) if internal nodes are ordered from left-to-right
	if (which.min(xx) != (ntips+1))
		{
		## 
		print("ERROR: reorder nodes the way ape/phylo expects! (tips first, then internals in order from left-to-right.")
		#yy[internal] <- rev(yy[!is.tip])[order(xx[!is.tip])]
		#xx[internal] <- rev(yy[!is.tip])[order(xx[!is.tip])]
		}

	edges <- matrix(nrow=0,ncol=2)
	edge.length <- numeric(0)
	nnode <- length(xx)-ntips

	while (sum(unlinked_nodes) > 1)
		{
		## find next node to include (the rightmost internal node)
		nextnode <- which(!is.tip & xx==max(xx[!is.tip]))[1]

		## find daughters
		
		# get the distance (in y) to all of the other corners
		ydist <- yy-yy[nextnode]
		xdist <- xx-xx[nextnode]

		# Check if it's the root
		if ( is.na(cy[nextnode]) )
			{
			cy[nextnode] = yy[nextnode]
			cx[nextnode] = 0			# leftmost coordinate must be 0!
			}
		
		cydist <- yy-cy[nextnode]
		cxdist <- xx-cx[nextnode]

		
		# find the TWO tips closest to this internal node, which are RIGHT of this node
		# this only finds the CLOSEST tip in Y, we want the TWO closest tips!!
		#daughters <- which(is.tip & dist==min(dist[is.tip]))
		
		# rank the ydistances in the y direction
		rank_of_ydists = order(cydist)
		
		# rank the xdistances in the x direction
		rank_of_xdists = order(cxdist)
		
		# get the node numbers in order; delete from this list as they are eliminated
		nodes_to_keep = nodes
		
		# daughter nodes must be to the right (in x) of the nextnode
		# (and they must be unlinked)
		nodes_to_keep = nodes_to_keep[unlinked_nodes][xdist[unlinked_nodes] > 0]
		
		# daughter nodes should be the two closest corners to nextnode (in y, mostly AND x)
		absolute_dist_from_node = 100*abs(cydist[nodes_to_keep]) + 1*abs(cxdist[nodes_to_keep])
		
		# sort the distances
		absolute_dist_from_node_sorted = sort(absolute_dist_from_node)
		
		# take the 2nd smallest absolute distance
		y_abs_dist_tokeep = absolute_dist_from_node_sorted[2]
		
		nodes_to_keep_final = nodes_to_keep[absolute_dist_from_node <= y_abs_dist_tokeep]
		print(paste("Linking: #", nodes_to_keep_final[1], " ", tip.labels[nodes_to_keep_final[1]], ", #", nodes_to_keep_final[2], " ", tip.labels[nodes_to_keep_final[2]], sep=""))
		
		#daughters <- which(is.tip & dist==min(dist[is.tip]))
		daughters = nodes_to_keep_final

		## be careful with numeric fuzz?
		edges <- rbind(edges,
					   nodes[c(nextnode,daughters[1])],
					   nodes[c(nextnode,daughters[2])])
		edge.length <- c(edge.length,xx[daughters]-xx[nextnode])

		# add nextnode to the list of tips (which are not available nodes for the nextnode)
		is.tip[nextnode] <- TRUE

		# remove the daughters & coordinates from the list of available nodes
		unlinked_nodes[daughters] = FALSE
		print(sum(unlinked_nodes))
		
		#xx <- xx[-daughters]
		#yy <- yy[-daughters]

		
		# remove the daughters from the list of possible nodes to link
		#unlinked_nodes
		
		#is.tip <- is.tip[-daughters]
		#nodes <- nodes[-daughters]
		}
	tr <- list(tip.labels=tip.labels,
			 edge=edges,
			 edge.length=edge.length,
			 Nnode=nnode)

	class(tr) <- "phylo"
	tr <- reorder(tr)
	tr$tip.labels = tip.labels
	return(tr)
	}



read_table_good <- function(fn)
	{
	# Read in the data, store in variable d
	# This has all of Nick's preferred options
	dtf = read.table(fn, header=TRUE, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
	return(dtf)
	}


# Remove factors from the character-like columns of a data.frame
# (leaves the numbers as numbers)
df_nonum_factors_to_char <- function(dtf, max_NAs=0.5)
	{
	dtf_classes = cls.df(dtf, printout=TRUE)
	
	dtf_names = names(dtf)
	numcols = ncol(dtf)
	
	cls_col_list = c()
	for (i in 1:numcols)
		{
		# Get one column:
		cmdstr = paste("cls_col = class(dtf$", dtf_names[i], ")", sep="")
		eval(parse(text = cmdstr))
		
		#cat(i, ": ", dtf_names[i], "	=	", cls_col, "\n", sep="")
		cls_col_list[i] = cls_col
		}
	
	for (i in 1:numcols)
		{
		if (cls_col_list[i] == "factor")
			{
			# Get one column, convert to character:
			cmdstr = paste("newcol = as.character(dtf$", dtf_names[i], ")", sep="")
			eval(parse(text = cmdstr))			
			
			cmdstr = paste("dtf$", dtf_names[i], " = newcol", sep="")
			eval(parse(text = cmdstr))				
			}
		}
	tmp_classes = cls.df(dtf)
	dtf_classes$newclasses = tmp_classes[,ncol(tmp_classes)]
	cat("\n")
	cat("dfnums_to_numeric(dtf, max_NAs=", max_NAs, ") reports: dataframe 'dtf_classes' has ", nrow(dtf_classes), " rows, ", ncol(dtf_classes), " columns.\n", sep="")
	cat("...names() and classes() of each column below...\n", sep="")
	cat("\n")
	print(dtf_classes)
	
	return(dtf)
	}

data = '

==================================================================
==================================================================


Example text files:

tipnames.txt:
=============
tipnames
ammiralis
dalli
textile
gloriamaris
aulicus
episcopatus
crocatus
omaria
furvus
bandanus
marmoreus
tessulatus
arenatus
pulicarius
radiatus
parius
laterculatus
aurisiacus
stercusmuscarum
consors
jaspideus
orbignyi
vimineus
OUTGROUPS
=============


tips.txt:
================
x	y
5.043	9.051
5.354	8.724
5.58	8.374
5.269	8.016
5.382	7.697
5.411	7.348
5.524	7.064
5.597	6.802
5.524	6.562
5.495	6.234
5.495	5.885
5.326	5.557
5.948	5.251
5.552	4.88
5.665	4.53
5.411	4.181
5.298	3.875
5.609	3.548
5.58	3.133
5.524	2.892
6.202	2.609
5.524	2.325
5.411	1.953
4.987	1.713
================


internal.txt:
============
x	y
4.732	8.658
4.874	8.549
4.591	8.309
4.421	7.632
4.93	7.523
4.761	7.195
4.789	6.933
4.224	6.889
4.591	6.824
4.082	6.212
5.298	6.059
5.298	5.076
3.913	5.033
4.28	4.421
4.987	4.356
4.676	4.115
2.697	4.006
4.082	3.744
4.902	3.329
4.619	3.133
3.969	2.39
4.506	2.128
1.199	2.63
============

corners.txt:
=================
x	y
4.732	9.061
4.874	8.724
4.874	8.374
4.732	8.549
4.45	8.309
4.591	8.025
4.93	7.697
4.93	7.348
4.761	7.523
4.619	7.195
4.817	7.086
4.817	6.824
4.619	6.562
4.252	6.059
4.082	5.557
4.308	5.076
4.308	4.115
4.704	4.356
4.987	4.53
4.987	4.181
4.676	3.875
4.93	3.548
4.93	3.133
4.619	2.892
4.619	3.351
1.199	4.006
1.228	1.713
2.714	2.387
2.697	5.055
3.981	2.605
3.969	2.137
4.506	2.325
4.506	1.953
4.082	3.133
4.082	4.421
4.252	7.632
3.884	6.212
4.082	6.911
5.298	6.212
5.298	5.863
5.354	5.251
5.354	4.902
3.884	3.744
4.591	8.658
4.761	6.955
4.45	6.824
=================
'