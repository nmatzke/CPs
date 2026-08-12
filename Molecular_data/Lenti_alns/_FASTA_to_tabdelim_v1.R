

#######################################################
# Purpose: read in FASTA sequences,
# write out: tab-delimited text file for use
# in Excel
#
# Location: 
# ~/GitHub/CPs/Molecular_data/Lenti_alns/_FASTA_to_tabdelim_v1.R
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

wd = "~/GitHub/CPs/Molecular_data/Lenti_alns/"
setwd(wd)
list.files()

#######################################################
# Concatenate 2 alignments
# (here, amino-acid and 3Di alignments, but it shouldn't matter)
#
# Assumptions:
# 1. In an alignment, each row/sequence has the same length
#######################################################

# Load the two alignments (for viewing)
fasta_fn1 = "Lentibulariacea_rps16_aln_v1.fasta"
fasta_fn2 = "Lentibulariacea_rbcL_aln_v1.fasta"
fasta_fn3 = "Lentibulariacea_matK_aln_v1.fasta"

aln1 = read_FASTA_safe(fasta_fn1, type="DNA")
aln2 = read_FASTA_safe(fasta_fn2, type="DNA")
aln3 = read_FASTA_safe(fasta_fn3, type="DNA")
aln1
length(aln1)
length(aln1[[1]])

# Goal: write these out to tab-delimited text files

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

# Write to text file via temporary data.frame (df)
tmpdf = cbind(seqnames1, strings1)
tmp_outfn = gsub(pattern=".fasta", replacement="_tab.txt", x=fasta_fn1, ignore.case=TRUE)
tmp = write.table(x=tmpdf, file=tmp_outfn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

tmpdf = cbind(seqnames2, strings2)
tmp_outfn = gsub(pattern=".fasta", replacement="_tab.txt", x=fasta_fn2, ignore.case=TRUE)
tmp = write.table(x=tmpdf, file=tmp_outfn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

tmpdf = cbind(seqnames3, strings3)
tmp_outfn = gsub(pattern=".fasta", replacement="_tab.txt", x=fasta_fn3, ignore.case=TRUE)
tmp = write.table(x=tmpdf, file=tmp_outfn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
