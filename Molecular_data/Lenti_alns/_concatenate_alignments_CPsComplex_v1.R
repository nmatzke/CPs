#######################################################
# Purpose: example of horizontally concatenating
# sequence alignments (hcat_align)
#
# Location: 
# ~/GitHub/pr/tutorials/concat_align/
# _concatenate_alignments_example_v1.R
#######################################################

installs='
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
install.packages("seqinr")
'

library(ape)
library(phytools)
library(BioGeoBEARS)
library(Biostrings)
library(seqinr) # for amb

sourceall("~/GitHub/str2phy/Rsrc/")

wd = "~/GitHub/pr/tutorials/concat_align/"
#wd = "~/GitHub/CPs/Molecular_data/Lenti_alns/"
setwd(wd)
list.files()

setwd(wd)
list.files()

# Load the 3 alignments (for viewing)
fasta_fn1 = "Lentibulariacea_matK_aln_v1.fasta"
fasta_fn2 = "Lentibulariacea_rps16_aln_v1.fasta"
fasta_fn3 = "Lentibulariacea_rbcL_aln_v1.fasta"

aln1 = read_FASTA_safe(fasta_fn1, type="DNA")
aln2 = read_FASTA_safe(fasta_fn2, type="DNA")
aln3 = read_FASTA_safe(fasta_fn3, type="DNA")
aln1
length(aln1)
length(aln1[[1]])

# Get sequences as strings
chars1 = as.character(aln1)
seqnames1 = names(chars1)
strings1 = sapply(FUN=paste0, X=chars1, collapse="")
chars2 = as.character(aln2)
seqnames2 = names(chars2)
strings2 = sapply(FUN=paste0, X=chars2, collapse="")
chars3 = as.character(aln3)
seqnames3 = names(chars3)
strings3 = sapply(FUN=paste0, X=chars3, collapse="")
ncols1 = length(aln1[[1]])
ncols2 = length(aln2[[1]])
ncols3 = length(aln3[[1]])

#######################################################
# Challenge: merge sequences based on some "key"
# in each sequence
#######################################################
seqnames = c(seqnames1, seqnames2, seqnames3)
sort(seqnames)

# In Masafumi's sequences, it appears that it is usually:
# Genus_species_seqID_gene
# ...but sometimes...
# _R_Genus_species_seqID_gene
# Genus_species_SB-2009_seqID_gene
# 
# Best guess for keys:
# 1. split off last fields, take the rest
# 2. split off the seqID (special function required due to NC_320982 seqIDs
# 3. remove these last 2 fields to get name
# 4. remove _R_ prefix if desired
seqnames1_minus_last = get_all_but_suffixes(tmptxts=seqnames1, split="_")
seqnames1_minus_last
lastSeqIDs1 = get_suffix_seqids_from_name(strings=seqnames1_minus_last, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
lastSeqIDs1
seqnames1_minus_last_minus_lastSeqID  = get_strings_without_seqID_suffix(strings=seqnames1_minus_last,  split="_", recodes=genbank_prefixes(), removetxt=c(">"))
seqnames1_minus_last_minus_lastSeqID
sort(seqnames1_minus_last_minus_lastSeqID)
rev(sort(table(seqnames1_minus_last_minus_lastSeqID)))
#   [1] "_R_Utricularia_albocaerulea"          "_R_Utricularia_purpurascens"          "_R_Utricularia_reticulata"           
#   [4] "Genlisea_aff_pygmaea_SB-2009"         "Genlisea_aff_repens_SB-2009"          "Genlisea_aff_violacea_SB-2009"       

seqnames2_minus_last = get_all_but_suffixes(tmptxts=seqnames2, split="_")
lastSeqIDs2 = get_suffix_seqids_from_name(strings=seqnames2_minus_last, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
seqnames2_minus_last_minus_lastSeqID  = get_strings_without_seqID_suffix(strings=seqnames2_minus_last,  split="_", recodes=genbank_prefixes(), removetxt=c(">"))
sort(seqnames2_minus_last_minus_lastSeqID)
rev(sort(table(seqnames2_minus_last_minus_lastSeqID)))

seqnames3_minus_last = get_all_but_suffixes(tmptxts=seqnames3, split="_")
lastSeqIDs3 = get_suffix_seqids_from_name(strings=seqnames3_minus_last, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
seqnames3_minus_last_minus_lastSeqID  = get_strings_without_seqID_suffix(strings=seqnames3_minus_last,  split="_", recodes=genbank_prefixes(), removetxt=c(">"))
sort(seqnames3_minus_last_minus_lastSeqID)
rev(sort(table(seqnames3_minus_last_minus_lastSeqID)))

# Combine to master list of OTUs (Operational Taxonomic Unit) (species but also guessed species, subspecies, molecular samples that no one has named yet, etc.)
keys = c(seqnames1_minus_last_minus_lastSeqID, seqnames2_minus_last_minus_lastSeqID, seqnames3_minus_last_minus_lastSeqID)
head(rev(sort(table(keys))))

keys = gsub(pattern="_R_", replacement="", x=keys, ignore.case=FALSE)
rev(sort(table(keys)))
length(keys) # 475 sequences

# Let's remove the _R_
keys = sort(unique(keys)) # alphabetical
keys = keys[rev(order(nchar(keys)))]
nchar(keys)
length(keys) # 301 unique keys


#######################################################
# Now that we have the keys, let's assemble sequences for each key
#######################################################
list_of_seqstrings = list(strings1, strings2, strings3)
list_of_seqnames = list(seqnames1, seqnames2, seqnames3)
list_of_lastSeqIDs = list(lastSeqIDs1, lastSeqIDs2, lastSeqIDs3)
list_of_ncols = list(ncols1, ncols2, ncols3)
merged_list_of_seqstrings = list()
merged_list_of_seqnames = list()

# Loop through gene alignments (j)
for (j in 1:length(list_of_seqstrings))
	{
	strings = list_of_seqstrings[[j]]
	seqnames = list_of_seqnames[[j]]
	lastSeqIDs = list_of_lastSeqIDs[[j]]
	ncols = list_of_ncols[[j]]
	tmp_strings = NULL
	tmp_names = NULL
	# Loop through keys (i)
	for (i in 1:length(keys))
		{
		TF = grepl(pattern=keys[i], x=seqnames)
		if (sum(TF) == 0)
			{
			tmp = paste0(rep("-", times=ncols), collapse="")
			tmpname = "NA"
			}
		if (sum(TF) == 1)
			{
			tmp = strings[TF]
			tmpname = lastSeqIDs[TF]
			
			# If you get 1 and only 1 hit, remove from the list of thing to search next time
			# (this works if we work from the longest key to the shortest key, meaning we find
			#  Pinguicula_vulgaris_whatever before we find Pinguicula_vulgaris)
			# reduce seqnames
			hitnum_to_remove = (1:length(TF))[TF]
			seqnames = seqnames[-1*hitnum_to_remove]
			lastSeqIDs = lastSeqIDs[-1*hitnum_to_remove]
			strings = strings[-1*hitnum_to_remove]
			}
		if (sum(TF) > 1)
			{
			tmp = strings[TF][1]
			tmpname = lastSeqIDs[TF][1]
			txt = paste0("WARNING: for gene #j=", j, ", key #i=", i, ", key='", keys[i], "', ", sum(TF), " seqnames matched. For now, taking 1st hit.")
			cat("\n")
			cat(txt)
			cat("\n")
			cat("lastSeqIDs[TF]:\n")
			print(lastSeqIDs[TF])
			warning(txt)
			}
		tmp_strings = c(tmp_strings, tmp)
		tmp_names = c(tmp_names, tmpname)
		} # END for (i in 1:length(keys))
	merged_list_of_seqstrings[[j]] = tmp_strings
	merged_list_of_seqnames[[j]] = tmp_names
	} # END for (j in 1:length(list_of_seqstrings))

merged_list_of_seqnames
IDs = unlist(merged_list_of_seqnames)
IDs_mat = matrix(data=IDs, ncol=length(merged_list_of_seqnames), byrow=FALSE)
final_seqnames_df = cbind(keys, IDs_mat)
final_seqnames = apply(X=final_seqnames_df, MARGIN=1, FUN=paste0, collapse="|")
final_seqnames
final_seqnames = paste0(">", final_seqnames, collapse=NULL)
final_seqnames
final_seqnames = paste0(final_seqnames, "\n", collapse=NULL)
final_seqnames

seqs = unlist(merged_list_of_seqstrings)
seqs_mat = matrix(data=seqs, ncol=length(merged_list_of_seqnames), byrow=FALSE)
final_seqs_df = cbind(final_seqnames, seqs_mat)
final_seqs = apply(X=final_seqs_df, MARGIN=1, FUN=paste0, collapse="")
final_seqs

zz = write.table(x=final_seqs, file="merged_3genes.fasta", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
