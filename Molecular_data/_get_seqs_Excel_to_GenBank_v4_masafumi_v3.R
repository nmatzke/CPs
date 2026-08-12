# Installation
# (hit 'n' for no when asked to update a zillion other packages)
install_txt='
library(BiocManager)
BiocManager::install(c("annotate", "Biostrings", "genomes"))
install.packages("rentrez") # for entrez_fetch
'



# Setup
library(seqinr)
library(BioGeoBEARS)	# for cls.df
library(gdata)			# for trim
library(XML)
library(parallel)
library(S4Vectors)
library(stats4)
library(IRanges)
library(XVector)
library(Biostrings)
library(BiocGenerics)
library(rentrez) # Fetching full records: entrez_fetch()
library(openxlsx)
library(ape)
library(annotate)		# for Bioconductor's 
						#"genbank" (efetch by 
						# gene Accession #s)
# BLAST and download related sequences


# HELP:
# http://www.ncbi.nlm.nih.gov/books/NBK3837/figure/EntrezHelp.F4/?report=objectonly
# http://news.open-bio.org/news/2009/06/ncbi-einfo-biopython/
# http://www.ncbi.nlm.nih.gov/books/NBK179288/
# http://biopython-documentation-chinese-translate.googlecode.com/git/_build/epub/ch8.html
# http://www.ncbi.nlm.nih.gov/Class/NAWBIS/Modules/InfoHubs/Exercises/infohubs_qa_taxon.html
# http://www.ncbi.nlm.nih.gov/refseq/rsg/about/

# Source this function:
# http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html
#source("/drives/GDrive/__GDrive_projects/2017-08-08_Matt_Baker_flagellum_speaker/01_BLAST/blastsequences_v4.R")

# Code to use:
source("~/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")


#Download Utricularia sequences
#wd = "/drives/GDrive/__GDrive_projects/2017-08-08_Matt_Baker_flagellum_speaker/01_BLAST/startdata/"
wd = "~/GitHub/masafumiPhD/Lentibulariaceae/Molecular_data/"
setwd(wd)


##############################################################
# Read in Masafumi's CP sequence IDs from Excel "Utricularia"
##############################################################
library(openxlsx)

xlsfn = "Lenti_seqs.xlsx"
xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet="Utricularia", startRow=1, colNames=TRUE)

######################################
# Things we want to download Utricularia$matK
######################################
IDs = xls$matK

# Subset xls table to just the rows for which we have matKs

# Convert NA to "/"
xls2 = xls
naTA = is.na(xls2$matK)
xls2$matK[naTA == TRUE] = "/"
TF = xls2$matK != "/"
xls2 = xls2[TF,]
xls2$matK
xls2$Species

#######################################################
# Download a list of sequence IDs, via a loop to label things
#######################################################

list_of_download_fns = NULL

cat(paste0("\nDownloading ", length(xls2$matK), " sequence IDs using rentrez::entrez_fetch()\n"))
for (i in 1:length(xls2$matK))
#for (i in 1:5) # test on just 5
	{
	cat(i, ",", sep="")
	seq_ids = c(xls2$matK[i])
	downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
	outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
	list_of_download_fns = c(list_of_download_fns, outfn)
	write(x=downloaded_seqs, file=outfn, sep="")
	} # END for (i in 1:length(xls2$matK))


# Check Utricularia gibba OR619942.1
downloaded_seqs = entrez_fetch(db="nuccore", id=c("OR619942.1"), rettype="fasta_cds_na")
# Write out as a once off
write(x=downloaded_seqs, file="Utricularia_gibba_OR619942.1.fasta", sep="")


# Look at example download
moref("Utricularia_tenuicaulis_NC_058517.1.fasta")
list.files()
list_of_download_fns

list_of_fastas = list.files()
TF1 = grepl(pattern="Utricularia_", x=list_of_fastas)
TF2 = grepl(pattern=".fasta", x=list_of_fastas)
TF = (TF1 + TF2) == 2
list_of_fastas = list_of_fastas[TF]
list_of_fastas

# Pull out the matKs
matK_seqs_txt = list()
list_i = 0
matK_seqs_names = NULL
for (i in 1:length(list_of_fastas))
	{
	fn = list_of_fastas[i]
	#seqs = ape::read.FASTA(file=fn, type="DNA") # THE 'type' IS VERY IMPORTANT HERE
	seqs = seqinr::read.fasta(file=fn, seqtype="DNA")
	#seqnames = names(seqs)
	
	# Just names doesnt work for seqinr seqs, which uses "attributes", specifically "Annot", to hold the full name
	# to get 1 attribute for the 1st sequence in the list:
	attr(x=seqs[[1]], which="Annot")
	
	# sapply applies/runs a function over each item in a list:
	seqnames = sapply(X=seqs, FUN=attr, which="Annot")

	
	# Search for gene=matK in the sequence headers
	TF1 = grepl(pattern="gene=matK", x=seqnames)
	TF2 = grepl(pattern="protein=maturase K", x=seqnames) # alternative way to identify matKs
	TF = (TF1 + TF2) >= 1
	
	# Error trap
	if (sum(TF) < 1)
		{
		stop_txt = paste0("Stop: no matK found for i=", i, ", fn=", fn)
		stop(stop_txt)
		}
	matching_names = seqnames[TF]
	print(matching_names)
	
	# Extract the matK
	nums = 1:length(TF)
	num_for_matK = nums[TF]
	#seqs # from APE, sequences in a "binary" format to save space
	#seqstxt = as.character(seqs) # sequences as letters
	#seqstxt = sapply(X=as.character(seqs), FUN=paste0, collapse="") # sequences as strings (desired)
	#matKtxt = seqstxt[[num_for_matK]]
	matKtxt = seqs[[num_for_matK]]
	matK_seqs_txt[[(list_i=list_i+1)]] = matKtxt
	matK_seqs_names = c(matK_seqs_names, matching_names)
	}
matK_seqs_names = unname(matK_seqs_names) # take off the names on the names
matK_seqs_names

# Remove first character (which is a >)
substr(x=matK_seqs_names[1], start=2, stop=nchar(matK_seqs_names[1]))
matK_seqs_names = sapply(X=matK_seqs_names, FUN=substr, start=2, stop=nchar(matK_seqs_names))

names(matK_seqs_txt) = matK_seqs_names
#matK_seqs_txt

#matK_seqs = as.DNAbin(matK_seqs_txt)
matK_seqs = matK_seqs_txt
matK_seqs

seqinr::write.fasta(sequences=matK_seqs, names=matK_seqs_names, file.out="all_Utric_matK_origLabels.fasta")

system("open all_Utric_matK_origLabels.fasta")



# Put in species names as sequence labels
spnames = gsub(pattern=".fasta", replacement="_matK", x=list_of_fastas)
spnames

seqinr::write.fasta(sequences=matK_seqs, names=spnames, file.out="all_Utric_matK.fasta")
system("open all_Utric_matK.fasta")





######################################
# Things we want to download Utricularia$rbcL
######################################
IDsrbcL = xls$rbcL
IDsrbcL
# Subset xls table to just the rows for which we have matKs

# Convert NA to "/"
xls2 = xls
naTA = is.na(xls2$rbcL)
xls2$rbcL[naTA == TRUE] = "/"
TF = xls2$rbcL != "/"
xls2 = xls2[TF,]
xls2$rbcL
xls2$Species

#######################################################
# Download a list of sequence IDs, via a loop to label things
#######################################################

list_of_download_fns = NULL

cat(paste0("\nDownloading ", length(xls2$rbcL), " sequence IDs using rentrez::entrez_fetch()\n"))
for (i in 1:length(xls2$rbcL))
  #for (i in 1:5) # test on just 5
{
  cat(i, ",", sep="")
  seq_ids = c(xls2$rbcL[i])
  downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
  outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
  list_of_download_fns = c(list_of_download_fns, outfn)
  write(x=downloaded_seqs, file=outfn, sep="")
} # END for (i in 1:length(xls2$rbcL))

# Look at example download
moref("Utricularia_tenuicaulis_NC_058517.1.fasta")
list.files()
list_of_download_fns

list_of_fastas = list.files()
TF1 = grepl(pattern="Utricularia_", x=list_of_fastas)
TF2 = grepl(pattern=".fasta", x=list_of_fastas)
TF = (TF1 + TF2) == 2
list_of_fastas = list_of_fastas[TF]
list_of_fastas

# Pull out the rbcL
rbcL_seqs_txt = list()
list_i = 0
rbcL_seqs_names = NULL
for (i in 1:length(list_of_fastas))
{
  fn = list_of_fastas[i]
  #seqs = ape::read.FASTA(file=fn, type="DNA") # THE 'type' IS VERY IMPORTANT HERE
  seqs = seqinr::read.fasta(file=fn, seqtype="DNA")
  #seqnames = names(seqs)
  
  # Just names doesnt work for seqinr seqs, which uses "attributes", specifically "Annot", to hold the full name
  # to get 1 attribute for the 1st sequence in the list:
  attr(x=seqs[[1]], which="Annot")
  
  # sapply applies/runs a function over each item in a list:
  seqnames = sapply(X=seqs, FUN=attr, which="Annot")
  
  
  # Search for gene=rbcL in the sequence headers
  TF1 = grepl(pattern="gene=rbcL", x=seqnames)
  TF2 = grepl(pattern="protein=ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit", x=seqnames) # alternative way to identify rbcLs
  TF3 = grepl(pattern="product=ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit", x=seqnames)
  TF = (TF1 + TF2 + TF3) >= 1
  
  # Error trap
  if (sum(TF) < 1)
  {
    stop_txt = paste0("Stop: no rbcL found for i=", i, ", fn=", fn)
    stop(stop_txt)
  }
  matching_names = seqnames[TF]
  print(matching_names)
  
  # Extract the rbcL
  nums = 1:length(TF)
  num_for_rbcL = nums[TF]
  #seqs # from APE, sequences in a "binary" format to save space
  #seqstxt = as.character(seqs) # sequences as letters
  #seqstxt = sapply(X=as.character(seqs), FUN=paste0, collapse="") # sequences as strings (desired)
  #matKtxt = seqstxt[[num_for_matK]]
  rbcLtxt = seqs[[num_for_rbcL]]
  rbcL_seqs_txt[[(list_i=list_i+1)]] = rbcLtxt
  rbcL_seqs_names = c(rbcL_seqs_names, matching_names)
}
rbcL_seqs_names = unname(rbcL_seqs_names) # take off the names on the names
rbcL_seqs_names

# Remove first character (which is a >)
substr(x=rbcL_seqs_names[1], start=2, stop=nchar(rbcL_seqs_names[1]))
rbcL_seqs_names = sapply(X=rbcL_seqs_names, FUN=substr, start=2, stop=nchar(rbcL_seqs_names))

names(rbcL_seqs_txt) = rbcL_seqs_names
#rbcL_seqs_txt

#matK_seqs = as.DNAbin(matK_seqs_txt)
rbcL_seqs = rbcL_seqs_txt
rbcL_seqs

seqinr::write.fasta(sequences=rbcL_seqs, names=rbcL_seqs_names, file.out="all_Utric_rbcL_origLabels.fasta")

system("open all_Utric_rbcL_origLabels.fasta")



# Put in species names as sequence labels
spnames = gsub(pattern=".fasta", replacement="_rbcL", x=list_of_fastas)
spnames

seqinr::write.fasta(sequences=rbcL_seqs, names=spnames, file.out="all_Utric_rbcL.fasta")
system("open all_Utric_rbcL.fasta")






######################################
# Things we want to download Utricularia$rps16
######################################

#not working, cannot downaload fast file from genbank, review codings
IDsrps16 = xls$rps16

# Subset xls table to just the rows for which we have matKs

# Convert NA to "/"
xls2 = xls
naTA = is.na(xls2$rps16)
xls2$rps16[naTA == TRUE] = "/"
TF = xls2$rps16 != "/"
xls2 = xls2[TF,]
xls2$rps16
xls2$Species

#######################################################
# Download a list of sequence IDs, via a loop to label things
#######################################################

list_of_download_fns = NULL

cat(paste0("\nDownloading ", length(xls2$rps16), " sequence IDs using rentrez::entrez_fetch()\n"))

for (i in 1:length(xls2$rps16)) {
  cat(i, ",", sep="")
  seq_ids = c(xls2$rps16[i])
  
  # 1. Try downloading CDS first
  downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
  
  # 2. Check if it's empty (trimming whitespace/newlines to be safe)
  if (trimws(downloaded_seqs) == "") {
    # If empty, fall back to standard fasta
    downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta")
  }
  
  outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
  list_of_download_fns = c(list_of_download_fns, outfn)
  write(x=downloaded_seqs, file=outfn, sep="")
} # END for


cat(paste0("\nDownloading ", length(xls2$rps16), " sequence IDs using rentrez::entrez_fetch()\n"))
for (i in 1:length(xls2$rps16))
  #for (i in 1:5) # test on just 5
{
  cat(i, ",", sep="")
  seq_ids = c(xls2$rps16[i])
  downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
  outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
  list_of_download_fns = c(list_of_download_fns, outfn)
  write(x=downloaded_seqs, file=outfn, sep="")
} # END for (i in 1:length(xls2$rps16))


# Look at example download
moref("Utricularia_adpressa_MZ981773.1.fasta")
list.files()
list_of_download_fns

list_of_fastas = list.files()
TF1 = grepl(pattern="Utricularia_", x=list_of_fastas)
TF2 = grepl(pattern=".fasta", x=list_of_fastas)
TF = (TF1 + TF2) == 2
list_of_fastas = list_of_fastas[TF]
list_of_fastas

# Pull out the rps16
rps16_seqs_txt = list()
list_i = 0
rps16_seqs_names = NULL
for (i in 1:length(list_of_fastas))
{
  fn = list_of_fastas[i]
  #seqs = ape::read.FASTA(file=fn, type="DNA") # THE 'type' IS VERY IMPORTANT HERE
  seqs = seqinr::read.fasta(file=fn, seqtype="DNA")
  #seqnames = names(seqs)
  
  # Just names doesnt work for seqinr seqs, which uses "attributes", specifically "Annot", to hold the full name
  # to get 1 attribute for the 1st sequence in the list:
  attr(x=seqs[[1]], which="Annot")
  
  # sapply applies/runs a function over each item in a list:
  seqnames = sapply(X=seqs, FUN=attr, which="Annot")
  seqnames
  
  # Search for gene=matK in the sequence headers
  TF1 = grepl(pattern="gene=rps16", x=seqnames) # alternative way to identify rps16
  TF2 = grepl(pattern="product=ribosomal protein S16", x=seqnames)
  TF3 = grepl(pattern="protein=ribosomal protein S16", x=seqnames)
  TF4 = grepl(pattern="rps16", x=seqnames)
  TF = (TF1 + TF2 + TF3 + TF4) >= 1
  
  # Error trap
  if (sum(TF) < 1)
  {
    stop_txt = paste0("Stop: no rps16 found for i=", i, ", fn=", fn)
    stop(stop_txt)
  }
  matching_names = seqnames[TF]
  print(matching_names)
  
  # Extract the rbcL
  nums = 1:length(TF)
  num_for_rps16 = nums[TF]
  #seqs # from APE, sequences in a "binary" format to save space
  #seqstxt = as.character(seqs) # sequences as letters
  #seqstxt = sapply(X=as.character(seqs), FUN=paste0, collapse="") # sequences as strings (desired)
  #matKtxt = seqstxt[[num_for_matK]]
  rps16txt = seqs[[num_for_rps16]]
  rps16_seqs_txt[[(list_i=list_i+1)]] = rps16txt
  rps16_seqs_names = c(rps16_seqs_names, matching_names)
}
rps16_seqs_names = unname(rps16_seqs_names) # take off the names on the names
rps16_seqs_names

# Remove first character (which is a >)
substr(x=rps16_seqs_names[1], start=2, stop=nchar(rps16_seqs_names[1]))
rps16_seqs_names = sapply(X=rps16_seqs_names, FUN=substr, start=2, stop=nchar(rps16_seqs_names))

names(rps16_seqs_txt) = rps16_seqs_names
#rps16_seqs_txt

#rps16_seqs = as.DNAbin(rps16_seqs_txt)
rps16_seqs = rps16_seqs_txt
rps16_seqs

seqinr::write.fasta(sequences=rps16_seqs, names=rps16_seqs_names, file.out="all_Utric_rps16_origLabels.fasta")

system("open all_Utric_rps16_origLabels.fasta")



# Put in species names as sequence labels
spnames = gsub(pattern=".fasta", replacement="_rps16", x=list_of_fastas)
spnames

seqinr::write.fasta(sequences=rps16_seqs, names=spnames, file.out="all_Utric_rps16.fasta")
system("open all_Utric_rps16.fasta")

























##############################################
# Things we want to download  Pinguicula$matK
##############################################

xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet="Pinguicula", startRow=1, colNames=TRUE)

IDs = xls$matK

# Subset xls table to just the rows for which we have matKs

# Convert NA to "/"
xls2 = xls
naTA = is.na(xls2$matK)
xls2$matK[naTA == TRUE] = "/"
TF = xls2$matK != "/"
xls2 = xls2[TF,]
xls2$matK
xls2$Species

#######################################################
# Download a list of sequence IDs, via a loop to label things
#######################################################

list_of_download_fns = NULL

cat(paste0("\nDownloading ", length(xls2$matK), " sequence IDs using rentrez::entrez_fetch()\n"))

for (i in 1:length(xls2$matK)) {
  cat(i, ",", sep="")
  seq_ids = c(xls2$matK[i])
  
  # 1. Try downloading CDS first
  downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
  
  # 2. Check if it's empty (trimming whitespace/newlines to be safe)
  if (trimws(downloaded_seqs) == "") {
    # If empty, fall back to standard fasta
    downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta")
  }
  
  outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
  list_of_download_fns = c(list_of_download_fns, outfn)
  write(x=downloaded_seqs, file=outfn, sep="")
} # END for


# Look at example download
list.files()
list_of_download_fns

list_of_fastas = list.files()
TF1 = grepl(pattern="Pinguicula_", x=list_of_fastas)
TF2 = grepl(pattern=".fasta", x=list_of_fastas)
TF = (TF1 + TF2) == 2
list_of_fastas = list_of_fastas[TF]
list_of_fastas

# Pull out the matK
matK_seqs_txt = list()
list_i = 0
matK_seqs_names = NULL
for (i in 1:length(list_of_fastas))
{
  fn = list_of_fastas[i]
  #seqs = ape::read.FASTA(file=fn, type="DNA") # THE 'type' IS VERY IMPORTANT HERE
  seqs = seqinr::read.fasta(file=fn, seqtype="DNA")
  #seqnames = names(seqs)
  
  # Just names doesnt work for seqinr seqs, which uses "attributes", specifically "Annot", to hold the full name
  # to get 1 attribute for the 1st sequence in the list:
  attr(x=seqs[[1]], which="Annot")
  
  # sapply applies/runs a function over each item in a list:
  seqnames = sapply(X=seqs, FUN=attr, which="Annot")
  seqnames
  
  # Search for gene=matK in the sequence headers
  TF1 = grepl(pattern="gene=matK", x=seqnames) # alternative way to identify matK
  TF2 = grepl(pattern="product=maturase K", x=seqnames)
  TF3 = grepl(pattern="protein=maturase K", x=seqnames)
  TF4 = grepl(pattern="matK", x=seqnames)
  TF = (TF1 + TF2 + TF3 + TF4) >= 1
  
  # Error trap
  if (sum(TF) < 1)
  {
    stop_txt = paste0("Stop: no matK found for i=", i, ", fn=", fn)
    stop(stop_txt)
  }
  matching_names = seqnames[TF]
  print(matching_names)
  
  # Extract the rbcL
  nums = 1:length(TF)
  num_for_matK = nums[TF]
  #seqs # from APE, sequences in a "binary" format to save space
  #seqstxt = as.character(seqs) # sequences as letters
  #seqstxt = sapply(X=as.character(seqs), FUN=paste0, collapse="") # sequences as strings (desired)
  #matKtxt = seqstxt[[num_for_matK]]
  matKtxt = seqs[[num_for_matK]]
  matK_seqs_txt[[(list_i=list_i+1)]] = matKtxt
  matK_seqs_names = c(matK_seqs_names, matching_names)
}
matK_seqs_names = unname(matK_seqs_names) # take off the names on the names
matK_seqs_names

# Remove first character (which is a >)
substr(x=matK_seqs_names[1], start=2, stop=nchar(matK_seqs_names[1]))
matK_seqs_names = sapply(X=matK_seqs_names, FUN=substr, start=2, stop=nchar(matK_seqs_names))

names(matK_seqs_txt) = matK_seqs_names
#matK_seqs_txt

#matK_seqs = as.DNAbin(matK_seqs_txt)
matK_seqs = matK_seqs_txt
matK_seqs

seqinr::write.fasta(sequences=matK_seqs, names=matK_seqs_names, file.out="all_Ping_matK_origLabels.fasta")

system("open all_Ping_matK_origLabels.fasta")



# Put in species names as sequence labels
spnames = gsub(pattern=".fasta", replacement="_matK", x=list_of_fastas)
spnames

seqinr::write.fasta(sequences=matK_seqs, names=spnames, file.out="all_Ping_matK.fasta")
system("open all_Ping_matK.fasta")




##############################################
# Things we want to download  Pinguicula$rbcL
##############################################
xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet="Pinguicula", startRow=1, colNames=TRUE)

IDs = xls$rbcL

# Subset xls table to just the rows for which we have rbcLs

# Convert NA to "/"
xls2 = xls
naTA = is.na(xls2$rbcL)
xls2$rbcL[naTA == TRUE] = "/"
TF = xls2$rbcL != "/"
xls2 = xls2[TF,]
xls2$rbcL
xls2$Species

#######################################################
# Download a list of sequence IDs, via a loop to label things
#######################################################

list_of_download_fns = NULL

cat(paste0("\nDownloading ", length(xls2$rbcL), " sequence IDs using rentrez::entrez_fetch()\n"))

for (i in 1:length(xls2$rbcL)) {
  cat(i, ",", sep="")
  seq_ids = c(xls2$rbcL[i])
  
  # 1. Try downloading CDS first
  downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
  
  # 2. Check if it's empty (trimming whitespace/newlines to be safe)
  if (trimws(downloaded_seqs) == "") {
    # If empty, fall back to standard fasta
    downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta")
  }
  
  outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
  list_of_download_fns = c(list_of_download_fns, outfn)
  write(x=downloaded_seqs, file=outfn, sep="")
} # END for


# Look at example download
list.files()
list_of_download_fns

list_of_fastas = list.files()
TF1 = grepl(pattern="Pinguicula_", x=list_of_fastas)
TF2 = grepl(pattern=".fasta", x=list_of_fastas)
TF = (TF1 + TF2) == 2
list_of_fastas = list_of_fastas[TF]
list_of_fastas

# Pull out the rbcL
rbcL_seqs_txt = list()
list_i = 0
rbcL_seqs_names = NULL
for (i in 1:length(list_of_fastas))
{
  fn = list_of_fastas[i]
  #seqs = ape::read.FASTA(file=fn, type="DNA") # THE 'type' IS VERY IMPORTANT HERE
  seqs = seqinr::read.fasta(file=fn, seqtype="DNA")
  #seqnames = names(seqs)
  
  # Just names doesnt work for seqinr seqs, which uses "attributes", specifically "Annot", to hold the full name
  # to get 1 attribute for the 1st sequence in the list:
  attr(x=seqs[[1]], which="Annot")
  
  # sapply applies/runs a function over each item in a list:
  seqnames = sapply(X=seqs, FUN=attr, which="Annot")
  seqnames
  
  # Search for gene=rbcL in the sequence headers
  TF1 = grepl(pattern="gene=rbcL", x=seqnames) # alternative way to identify rbcL
  TF2 = grepl(pattern="protein=ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit", x=seqnames)
  TF3 = grepl(pattern="product=ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit", x=seqnames)
  TF4 = grepl(pattern="rbcL", x=seqnames)
  TF = (TF1 + TF2 + TF3 + TF4) >= 1
  
  # Error trap
  if (sum(TF) < 1)
  {
    stop_txt = paste0("Stop: no rbcL found for i=", i, ", fn=", fn)
    stop(stop_txt)
  }
  matching_names = seqnames[TF]
  print(matching_names)
  
  # Extract the rbcL
  nums = 1:length(TF)
  num_for_rbcL = nums[TF]
  #seqs # from APE, sequences in a "binary" format to save space
  #seqstxt = as.character(seqs) # sequences as letters
  #seqstxt = sapply(X=as.character(seqs), FUN=paste0, collapse="") # sequences as strings (desired)
  #rbcLtxt = seqstxt[[num_for_rbcL]]
  rbcLtxt = seqs[[num_for_rbcL]]
  rbcL_seqs_txt[[(list_i=list_i+1)]] = rbcLtxt
  rbcL_seqs_names = c(rbcL_seqs_names, matching_names)
}
rbcL_seqs_names = unname(rbcL_seqs_names) # take off the names on the names
rbcL_seqs_names

# Remove first character (which is a >)
substr(x=rbcL_seqs_names[1], start=2, stop=nchar(rbcL_seqs_names[1]))
rbcL_seqs_names = sapply(X=rbcL_seqs_names, FUN=substr, start=2, stop=nchar(rbcL_seqs_names))

names(rbcL_seqs_txt) = rbcL_seqs_names
#rbcL_seqs_txt

#rbcL_seqs = as.DNAbin(rbcL_seqs_txt)
rbcL_seqs = rbcL_seqs_txt
rbcL_seqs

seqinr::write.fasta(sequences=rbcL_seqs, names=rbcL_seqs_names, file.out="all_Ping_rbcL_origLabels.fasta")

system("open all_Ping_rbcL_origLabels.fasta")



# Put in species names as sequence labels
spnames = gsub(pattern=".fasta", replacement="_rbcL", x=list_of_fastas)
spnames

seqinr::write.fasta(sequences=rbcL_seqs, names=spnames, file.out="all_Ping_rbcL.fasta")
system("open all_Ping_rbcL.fasta")





######################################
# Things we want to download Pinguicula$rps16
######################################
xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet="Pinguicula", startRow=1, colNames=TRUE)

IDsrps16 = xls$rps16

# Subset xls table to just the rows for which we have matKs

# Convert NA to "/"
xls2 = xls
naTA = is.na(xls2$rps16)
xls2$rps16[naTA == TRUE] = "/"
TF = xls2$rps16 != "/"
xls2 = xls2[TF,]
xls2$rps16
xls2$Species

#######################################################
# Download a list of sequence IDs, via a loop to label things
#######################################################

list_of_download_fns = NULL

cat(paste0("\nDownloading ", length(xls2$rps16), " sequence IDs using rentrez::entrez_fetch()\n"))

for (i in 1:length(xls2$rps16)) {
  cat(i, ",", sep="")
  seq_ids = c(xls2$rps16[i])
  
  # 1. Try downloading CDS first
  downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
  
  # 2. Check if it's empty (trimming whitespace/newlines to be safe)
  if (trimws(downloaded_seqs) == "") {
    # If empty, fall back to standard fasta
    downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta")
  }
  
  outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
  list_of_download_fns = c(list_of_download_fns, outfn)
  write(x=downloaded_seqs, file=outfn, sep="")
} # END for


list.files()
list_of_download_fns

list_of_fastas = list.files()
TF1 = grepl(pattern="Pinguicula_", x=list_of_fastas)
TF2 = grepl(pattern=".fasta", x=list_of_fastas)
TF = (TF1 + TF2) == 2
list_of_fastas = list_of_fastas[TF]
list_of_fastas

# Pull out the rps16
rps16_seqs_txt = list()
list_i = 0
rps16_seqs_names = NULL
for (i in 1:length(list_of_fastas))
{
  fn = list_of_fastas[i]
  #seqs = ape::read.FASTA(file=fn, type="DNA") # THE 'type' IS VERY IMPORTANT HERE
  seqs = seqinr::read.fasta(file=fn, seqtype="DNA")
  #seqnames = names(seqs)
  
  # Just names doesnt work for seqinr seqs, which uses "attributes", specifically "Annot", to hold the full name
  # to get 1 attribute for the 1st sequence in the list:
  attr(x=seqs[[1]], which="Annot")
  
  # sapply applies/runs a function over each item in a list:
  seqnames = sapply(X=seqs, FUN=attr, which="Annot")
  seqnames
  
  # Search for gene=matK in the sequence headers
  TF1 = grepl(pattern="gene=rps16", x=seqnames) # alternative way to identify rps16
  TF2 = grepl(pattern="product=ribosomal protein S16", x=seqnames)
  TF3 = grepl(pattern="protein=ribosomal protein S16", x=seqnames)
  TF4 = grepl(pattern="rps16", x=seqnames)
  TF = (TF1 + TF2 + TF3 + TF4) >= 1
  
  # Error trap
  if (sum(TF) < 1)
  {
    stop_txt = paste0("Stop: no rps16 found for i=", i, ", fn=", fn)
    stop(stop_txt)
  }
  matching_names = seqnames[TF]
  print(matching_names)
  
  # Extract the rbcL
  nums = 1:length(TF)
  num_for_rps16 = nums[TF]
  #seqs # from APE, sequences in a "binary" format to save space
  #seqstxt = as.character(seqs) # sequences as letters
  #seqstxt = sapply(X=as.character(seqs), FUN=paste0, collapse="") # sequences as strings (desired)
  #matKtxt = seqstxt[[num_for_matK]]
  rps16txt = seqs[[num_for_rps16]]
  rps16_seqs_txt[[(list_i=list_i+1)]] = rps16txt
  rps16_seqs_names = c(rps16_seqs_names, matching_names)
}
rps16_seqs_names = unname(rps16_seqs_names) # take off the names on the names
rps16_seqs_names

# Remove first character (which is a >)
substr(x=rps16_seqs_names[1], start=2, stop=nchar(rps16_seqs_names[1]))
rps16_seqs_names = sapply(X=rps16_seqs_names, FUN=substr, start=2, stop=nchar(rps16_seqs_names))

names(rps16_seqs_txt) = rps16_seqs_names
#rps16_seqs_txt

#rps16_seqs = as.DNAbin(rps16_seqs_txt)
rps16_seqs = rps16_seqs_txt
rps16_seqs

seqinr::write.fasta(sequences=rps16_seqs, names=rps16_seqs_names, file.out="all_Ping_rps16_origLabels.fasta")

system("open all_Ping_rps16_origLabels.fasta")



# Put in species names as sequence labels
spnames = gsub(pattern=".fasta", replacement="_rps16", x=list_of_fastas)
spnames

seqinr::write.fasta(sequences=rps16_seqs, names=spnames, file.out="all_Ping_rps16.fasta")
system("open all_Ping_rps16.fasta")






















##############################################################
# Read in Masafumi's CP sequence IDs from Excel "Genlisea"
##############################################################
library(openxlsx)

xlsfn = "Lenti_seqs.xlsx"

##############################################
# Things we want to download  Genlisea$matK
##############################################

xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet="Genlisea", startRow=1, colNames=TRUE)

IDs = xls$matK

# Subset xls table to just the rows for which we have matKs

# Convert NA to "/"
xls2 = xls
naTA = is.na(xls2$matK)
xls2$matK[naTA == TRUE] = "/"
TF = xls2$matK != "/"
xls2 = xls2[TF,]
xls2$matK
xls2$Species

#######################################################
# Download a list of sequence IDs, via a loop to label things
#######################################################

list_of_download_fns = NULL

cat(paste0("\nDownloading ", length(xls2$matK), " sequence IDs using rentrez::entrez_fetch()\n"))

for (i in 1:length(xls2$matK)) {
  cat(i, ",", sep="")
  seq_ids = c(xls2$matK[i])
  
  # 1. Try downloading CDS first
  downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
  
  # 2. Check if it's empty (trimming whitespace/newlines to be safe)
  if (trimws(downloaded_seqs) == "") {
    # If empty, fall back to standard fasta
    downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta")
  }
  
  outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
  list_of_download_fns = c(list_of_download_fns, outfn)
  write(x=downloaded_seqs, file=outfn, sep="")
} # END for


# Look at example download
list.files()
list_of_download_fns

list_of_fastas = list.files()
TF1 = grepl(pattern="Genlisea_", x=list_of_fastas)
TF2 = grepl(pattern=".fasta", x=list_of_fastas)
TF = (TF1 + TF2) == 2
list_of_fastas = list_of_fastas[TF]
list_of_fastas

# Pull out the matK
matK_seqs_txt = list()
list_i = 0
matK_seqs_names = NULL
for (i in 1:length(list_of_fastas))
{
  fn = list_of_fastas[i]
  #seqs = ape::read.FASTA(file=fn, type="DNA") # THE 'type' IS VERY IMPORTANT HERE
  seqs = seqinr::read.fasta(file=fn, seqtype="DNA")
  #seqnames = names(seqs)
  
  # Just names doesnt work for seqinr seqs, which uses "attributes", specifically "Annot", to hold the full name
  # to get 1 attribute for the 1st sequence in the list:
  attr(x=seqs[[1]], which="Annot")
  
  # sapply applies/runs a function over each item in a list:
  seqnames = sapply(X=seqs, FUN=attr, which="Annot")
  seqnames
  
  # Search for gene=matK in the sequence headers
  TF1 = grepl(pattern="gene=matK", x=seqnames) # alternative way to identify matK
  TF2 = grepl(pattern="product=maturase K", x=seqnames)
  TF3 = grepl(pattern="protein=maturase K", x=seqnames)
  TF4 = grepl(pattern="matK", x=seqnames)
  TF = (TF1 + TF2 + TF3 + TF4) >= 1
  
  # Error trap
  if (sum(TF) < 1)
  {
    stop_txt = paste0("Stop: no matK found for i=", i, ", fn=", fn)
    stop(stop_txt)
  }
  matching_names = seqnames[TF]
  print(matching_names)
  
  # Extract the matK
  nums = 1:length(TF)
  num_for_matK = nums[TF]
  #seqs # from APE, sequences in a "binary" format to save space
  #seqstxt = as.character(seqs) # sequences as letters
  #seqstxt = sapply(X=as.character(seqs), FUN=paste0, collapse="") # sequences as strings (desired)
  #matKtxt = seqstxt[[num_for_matK]]
  matKtxt = seqs[[num_for_matK]]
  matK_seqs_txt[[(list_i=list_i+1)]] = matKtxt
  matK_seqs_names = c(matK_seqs_names, matching_names)
}
matK_seqs_names = unname(matK_seqs_names) # take off the names on the names
matK_seqs_names

# Remove first character (which is a >)
substr(x=matK_seqs_names[1], start=2, stop=nchar(matK_seqs_names[1]))
matK_seqs_names = sapply(X=matK_seqs_names, FUN=substr, start=2, stop=nchar(matK_seqs_names))

names(matK_seqs_txt) = matK_seqs_names
#matK_seqs_txt

#matK_seqs = as.DNAbin(matK_seqs_txt)
matK_seqs = matK_seqs_txt
matK_seqs

seqinr::write.fasta(sequences=matK_seqs, names=matK_seqs_names, file.out="all_Genli_matK_origLabels.fasta")

system("open all_Genli_matK_origLabels.fasta")



# Put in species names as sequence labels
spnames = gsub(pattern=".fasta", replacement="_matK", x=list_of_fastas)
spnames

seqinr::write.fasta(sequences=matK_seqs, names=spnames, file.out="all_Genli_matK.fasta")
system("open all_Genli_matK.fasta")





##############################################
# Things we want to download  Genlisea$rbcL
##############################################

xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet="Genlisea", startRow=1, colNames=TRUE)

IDs = xls$rbcL

# Subset xls table to just the rows for which we have rbcLs

# Convert NA to "/"
xls2 = xls
naTA = is.na(xls2$rbcL)
xls2$rbcL[naTA == TRUE] = "/"
TF = xls2$rbcL != "/"
xls2 = xls2[TF,]
xls2$rbcL
xls2$Species

#######################################################
# Download a list of sequence IDs, via a loop to label things
#######################################################

list_of_download_fns = NULL

cat(paste0("\nDownloading ", length(xls2$rbcL), " sequence IDs using rentrez::entrez_fetch()\n"))

for (i in 1:length(xls2$rbcL)) {
  cat(i, ",", sep="")
  seq_ids = c(xls2$rbcL[i])
  
  # 1. Try downloading CDS first
  downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
  
  # 2. Check if it's empty (trimming whitespace/newlines to be safe)
  if (trimws(downloaded_seqs) == "") {
    # If empty, fall back to standard fasta
    downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta")
  }
  
  outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
  list_of_download_fns = c(list_of_download_fns, outfn)
  write(x=downloaded_seqs, file=outfn, sep="")
} # END for


# Look at example download
list.files()
list_of_download_fns

list_of_fastas = list.files()
TF1 = grepl(pattern="Genlisea_", x=list_of_fastas)
TF2 = grepl(pattern=".fasta", x=list_of_fastas)
TF = (TF1 + TF2) == 2
list_of_fastas = list_of_fastas[TF]
list_of_fastas

# Pull out the rbcL
rbcL_seqs_txt = list()
list_i = 0
rbcL_seqs_names = NULL
for (i in 1:length(list_of_fastas))
{
  fn = list_of_fastas[i]
  #seqs = ape::read.FASTA(file=fn, type="DNA") # THE 'type' IS VERY IMPORTANT HERE
  seqs = seqinr::read.fasta(file=fn, seqtype="DNA")
  #seqnames = names(seqs)
  
  # Just names doesnt work for seqinr seqs, which uses "attributes", specifically "Annot", to hold the full name
  # to get 1 attribute for the 1st sequence in the list:
  attr(x=seqs[[1]], which="Annot")
  
  # sapply applies/runs a function over each item in a list:
  seqnames = sapply(X=seqs, FUN=attr, which="Annot")
  seqnames
  
  # Search for gene=rbcL in the sequence headers
  TF1 = grepl(pattern="gene=rbcL", x=seqnames) # alternative way to identify rbcL
  TF2 = grepl(pattern="protein=ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit", x=seqnames)
  TF3 = grepl(pattern="product=ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit", x=seqnames)
  TF4 = grepl(pattern="rbcL", x=seqnames)
  TF = (TF1 + TF2 + TF3 + TF4) >= 1
  
  # Error trap
  if (sum(TF) < 1)
  {
    stop_txt = paste0("Stop: no rbcL found for i=", i, ", fn=", fn)
    stop(stop_txt)
  }
  matching_names = seqnames[TF]
  print(matching_names)
  
  # Extract the rbcL
  nums = 1:length(TF)
  num_for_rbcL = nums[TF]
  #seqs # from APE, sequences in a "binary" format to save space
  #seqstxt = as.character(seqs) # sequences as letters
  #seqstxt = sapply(X=as.character(seqs), FUN=paste0, collapse="") # sequences as strings (desired)
  #rbcLtxt = seqstxt[[num_for_rbcL]]
  rbcLtxt = seqs[[num_for_rbcL]]
  rbcL_seqs_txt[[(list_i=list_i+1)]] = rbcLtxt
  rbcL_seqs_names = c(rbcL_seqs_names, matching_names)
}
rbcL_seqs_names = unname(rbcL_seqs_names) # take off the names on the names
rbcL_seqs_names

# Remove first character (which is a >)
substr(x=rbcL_seqs_names[1], start=2, stop=nchar(rbcL_seqs_names[1]))
rbcL_seqs_names = sapply(X=rbcL_seqs_names, FUN=substr, start=2, stop=nchar(rbcL_seqs_names))

names(rbcL_seqs_txt) = rbcL_seqs_names
#rbcL_seqs_txt

#rbcL_seqs = as.DNAbin(rbcL_seqs_txt)
rbcL_seqs = rbcL_seqs_txt
rbcL_seqs

seqinr::write.fasta(sequences=rbcL_seqs, names=rbcL_seqs_names, file.out="all_Genli_rbcL_origLabels.fasta")

system("open all_Genli_rbcL_origLabels.fasta")



# Put in species names as sequence labels
spnames = gsub(pattern=".fasta", replacement="_rbcL", x=list_of_fastas)
spnames

seqinr::write.fasta(sequences=rbcL_seqs, names=spnames, file.out="all_Genli_rbcL.fasta")
system("open all_Genli_rbcL.fasta")






##############################################
# Things we want to download  Genlisea$rps16
##############################################

xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet="Genlisea", startRow=1, colNames=TRUE)

IDs = xls$rps16

# Subset xls table to just the rows for which we have rps16s

# Convert NA to "/"
xls2 = xls
naTA = is.na(xls2$rps16)
xls2$rps16[naTA == TRUE] = "/"
TF = xls2$rps16 != "/"
xls2 = xls2[TF,]
xls2$rps16
xls2$Species

#######################################################
# Download a list of sequence IDs, via a loop to label things
#######################################################

list_of_download_fns = NULL

cat(paste0("\nDownloading ", length(xls2$rps16), " sequence IDs using rentrez::entrez_fetch()\n"))

for (i in 1:length(xls2$rps16)) {
  cat(i, ",", sep="")
  seq_ids = c(xls2$rps16[i])
  
  # 1. Try downloading CDS first
  downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
  
  # 2. Check if it's empty (trimming whitespace/newlines to be safe)
  if (trimws(downloaded_seqs) == "") {
    # If empty, fall back to standard fasta
    downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta")
  }
  
  outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
  list_of_download_fns = c(list_of_download_fns, outfn)
  write(x=downloaded_seqs, file=outfn, sep="")
} # END for


# Look at example download
list.files()
list_of_download_fns

list_of_fastas = list.files()
TF1 = grepl(pattern="Genlisea_", x=list_of_fastas)
TF2 = grepl(pattern=".fasta", x=list_of_fastas)
TF = (TF1 + TF2) == 2
list_of_fastas = list_of_fastas[TF]
list_of_fastas

# Pull out the rps16
rps16_seqs_txt = list()
list_i = 0
rps16_seqs_names = NULL
for (i in 1:length(list_of_fastas))
{
  fn = list_of_fastas[i]
  #seqs = ape::read.FASTA(file=fn, type="DNA") # THE 'type' IS VERY IMPORTANT HERE
  seqs = seqinr::read.fasta(file=fn, seqtype="DNA")
  #seqnames = names(seqs)
  
  # Just names doesnt work for seqinr seqs, which uses "attributes", specifically "Annot", to hold the full name
  # to get 1 attribute for the 1st sequence in the list:
  attr(x=seqs[[1]], which="Annot")
  
  # sapply applies/runs a function over each item in a list:
  seqnames = sapply(X=seqs, FUN=attr, which="Annot")
  seqnames
  
  # Search for gene=rps16 in the sequence headers
  TF1 = grepl(pattern="gene=rps16", x=seqnames) # alternative way to identify rps16
  TF2 = grepl(pattern="product=ribosomal protein S16", x=seqnames)
  TF3 = grepl(pattern="protein=ribosomal protein S16", x=seqnames)
  TF4 = grepl(pattern="rps16", x=seqnames)
  TF = (TF1 + TF2 + TF3 + TF4) >= 1
  
  # Error trap
  if (sum(TF) < 1)
  {
    stop_txt = paste0("Stop: no rps16 found for i=", i, ", fn=", fn)
    stop(stop_txt)
  }
  matching_names = seqnames[TF]
  print(matching_names)
  
  # Extract the rps16
  nums = 1:length(TF)
  num_for_rps16 = nums[TF]
  #seqs # from APE, sequences in a "binary" format to save space
  #seqstxt = as.character(seqs) # sequences as letters
  #seqstxt = sapply(X=as.character(seqs), FUN=paste0, collapse="") # sequences as strings (desired)
  #rps16txt = seqstxt[[num_for_rps16]]
  rps16txt = seqs[[num_for_rps16]]
  rps16_seqs_txt[[(list_i=list_i+1)]] = rps16txt
  rps16_seqs_names = c(rps16_seqs_names, matching_names)
}
rps16_seqs_names = unname(rps16_seqs_names) # take off the names on the names
rps16_seqs_names

# Remove first character (which is a >)
substr(x=rps16_seqs_names[1], start=2, stop=nchar(rps16_seqs_names[1]))
rps16_seqs_names = sapply(X=rps16_seqs_names, FUN=substr, start=2, stop=nchar(rps16_seqs_names))

names(rps16_seqs_txt) = rps16_seqs_names
#rps16_seqs_txt

#rps16_seqs = as.DNAbin(rps16_seqs_txt)
rps16_seqs = rps16_seqs_txt
rps16_seqs

seqinr::write.fasta(sequences=rps16_seqs, names=rps16_seqs_names, file.out="all_Genli_rps16_origLabels.fasta")

system("open all_Genli_rps16_origLabels.fasta")



# Put in species names as sequence labels
spnames = gsub(pattern=".fasta", replacement="_rps16", x=list_of_fastas)
spnames

seqinr::write.fasta(sequences=rps16_seqs, names=spnames, file.out="all_Genli_rps16.fasta")
system("open all_Genli_rps16.fasta")




















##############################################
# Combine 3 genera into one fasta file "matK"
##############################################

# 1. Read the three separate genus FASTA files
genus1 <- readDNAStringSet("all_Utric_matK.fasta")
genus2 <- readDNAStringSet("all_Ping_matK.fasta")
genus3 <- readDNAStringSet("all_Genli_matK.fasta")

# 2. Combine them into one single object
all_genera <- c(genus1, genus2, genus3)

# 3. Take a quick look at your sequence names
names(all_genera)

# 4. Write the combined data to a new FASTA file
writeXStringSet(all_genera, filepath = "Lentibulariacea_matK.fasta")

system("open Lentibulariacea_matK.fasta")






##############################################
# Combine 3 genera into one fasta file "rbcL"
##############################################

# 1. Read the three separate genus FASTA files
genus1 <- readDNAStringSet("all_Utric_rbcL.fasta")
genus2 <- readDNAStringSet("all_Ping_rbcL.fasta")
genus3 <- readDNAStringSet("all_Genli_rbcL.fasta")

# 2. Combine them into one single object
all_genera <- c(genus1, genus2, genus3)

# 3. Take a quick look at your sequence names
names(all_genera)

# 4. Write the combined data to a new FASTA file
writeXStringSet(all_genera, filepath = "Lentibulariacea_rbcL.fasta")

system("open Lentibulariacea_rbcL.fasta")






##############################################
# Combine 3 genera into one fasta file "rps16"
##############################################

# 1. Read the three separate genus FASTA files
genus1 <- readDNAStringSet("all_Utric_rps16.fasta")
genus2 <- readDNAStringSet("all_Ping_rps16.fasta")
genus3 <- readDNAStringSet("all_Genli_rps16.fasta")

# 2. Combine them into one single object
all_genera <- c(genus1, genus2, genus3)

# 3. Take a quick look at your sequence names
names(all_genera)

# 4. Write the combined data to a new FASTA file
writeXStringSet(all_genera, filepath = "Lentibulariacea_rps16.fasta")

system("open Lentibulariacea_rps16.fasta")




#Running IQtree, below is the command I used
#bin/iqtree3 -s Lentibulariacea_matK_aln_v1.fasta -m MFP -B 1000 -T AUTO
#bin/iqtree3 -s Lentibulariacea_rbcL_aln_v1.fasta -m MFP -B 1000 -T AUTO
#bin/iqtree3 -s Lentibulariacea_rps16_aln_v1.fasta -m MFP -B 1000 -T AUTO






###################################
###### Concatnate alignments ######
###################################

# ---- 1. Read in your gene alignments ----
files <- c("Lentibulariacea_matK_aln_v1.fasta", "Lentibulariacea_rbcL_aln_v1.fasta", "Lentibulariacea_rps16_aln_v1.fasta")
alignments <- lapply(files, read.FASTA)  # returns list of DNAbin objects
names(alignments) <- c("matK", "rbcL", "rps16")

# ---- 2. Get the full, sorted species list ----
all_species <- sort(unique(unlist(lapply(alignments, names))))

# ---- 3. Build a padded matrix for each gene ----
pad_alignment <- function(dnabin, species_list) {
  aln_len <- ncol(as.matrix(dnabin))       # alignment width for this gene
  present <- names(dnabin)
  
  # Create output list, filling gaps for missing species
  out <- vector("list", length(species_list))
  names(out) <- species_list
  
  for (sp in species_list) {
    if (sp %in% present) {
      out[[sp]] <- as.character(dnabin[[sp]])
    } else {
      out[[sp]] <- rep("-", aln_len)       # missing = all gaps
    }
  }
  
  # Convert to matrix (species x sites)
  do.call(rbind, out)
}

padded_list <- lapply(alignments, pad_alignment, species_list = all_species)

# ---- 4. Concatenate across genes (column-wise) ----
supermatrix <- do.call(cbind, padded_list)
rownames(supermatrix) <- all_species

# ---- 5. Convert back to DNAbin and write out ----
supermatrix_dnabin <- as.DNAbin(supermatrix)
write.FASTA(supermatrix_dnabin, "concatenated_alignment.fasta")

# Optional: write a partition file for downstream phylogenetics (RAxML/IQ-TREE)
gene_lengths <- sapply(alignments, function(x) ncol(as.matrix(x)))
starts <- cumsum(c(1, head(gene_lengths, -1) + cumsum(head(gene_lengths, -1))*0))







































#

######################################################
# OLD
#######################################################
# Database fields (db)
# Table 1
# – Valid values of &retmode and &rettype for EFetch (null = empty string)
# https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/

# protein
#downloaded_seqs = entrez_fetch(db="protein", id=seq_ids, rettype="fasta")

# get FASTA file of CDS (coding sequences) previously identified by genome project
# rettype = return type = type of file to return
# CDS nucleotide FASTA
downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")

outfn = "CDS_nucleotide_FASTA.fasta"
write(x=downloaded_seqs, file=outfn, sep="")
opd()

#######################################################
# Get sequences by BLAST
#######################################################

install_cmds='
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("rBLAST", repos = "https://mhahsler.r-universe.dev")
' # END install_cmds

library(Biostrings)

library(rBLAST)