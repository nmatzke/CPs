
notes='
2026-09-28
~/GitHub/CPs/ex_cpDNA_genome/NC_023463.1_Pinguicula_ehlersiae_chloroplast/

Annotated
NC_023463.1	Pinguicula ehlersiae chloroplast, complete genome	Pinguicula ehlersiae	319	391	97%	5e-89	91.52%	147147	NC_023463.1

NC_023463.1_Pinguicula_ehlersiae_chloroplast.fasta

Annotated at:
CHLOROBOX
GeSeq
https://chlorobox.mpimp-golm.mpg.de/geseq.html

THis produced GFF (gene finder format) + FASTA file

GeSeqJob-20260901-120850_NC_023463.1-Pinguicula-ehlersiae-chloroplast\%2C-complete-genome_GFF3\ +\ FASTA.gff3 
->
NC_023463.1_Pinguicula_ehlersiae_chloroplast_gff3+FASTA.txt


Parse to individual FASTAs in R:

~/GitHub/CPs/ex_cpDNA_genome/NC_023463.1_Pinguicula_ehlersiae_chloroplast/
'



setup='
install.packages("readr")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
	{
	BiocManager::install(c("GenomicFeatures", "Biostrings"))
	BiocManager::install("txdbmaker")
	}
' # END setup

library(GenomicFeatures)
library(Biostrings)
library(BioGeoBEARS)
library(Rsamtools)
		library(readr)

wd = "~/GitHub/CPs/ex_cpDNA_genome/NC_023463.1_Pinguicula_ehlersiae_chloroplast/"
setwd(wd)

# 1. Define paths to your input files
gb_file <- "NC_023463.1_Pinguicula_ehlersiae_chloroplast_GeSeqAnnotated.gb"

gff_file   <- "NC_023463.1_Pinguicula_ehlersiae_chloroplast_gff3+FASTA.gff3"
fasta_file <- "NC_023463.1_Pinguicula_ehlersiae_chloroplast_gff3+FASTA.fasta"
output_file <- "extracted_sequences.fasta"

moref(gff_file)

# 1a. Create an index file (.fai)
fasta_file_indexed = indexFa(fasta_file)
# 1b. Open it as a FaFile object
fa_file <- FaFile(fasta_file)

# 2. Build a TxDb object from the GFF file
# This structures your genome annotations logically
txdb <- txdbmaker::makeTxDbFromGFF(gff_file, format = "gff")

# 3. Load the full genome FASTA file into memory
#genome_seqs <- readDNAStringSet(fa_file)

# 4. Extract your features of interest
# Change 'cdsBy' to 'genes', 'exonsBy', or 'transcripts' depending on what you need
features <- cdsBy(txdb, by = "tx", use.names = TRUE)

# 4. Now getSeq() will work perfectly!
# sub_records <- getSeq(fa_file, gr)

# 5. Extract the individual fasta sequences matching the GFF coordinates
extracted_seqs <- extractTranscriptSeqs(x=fa_file, transcripts=features)

# 6. Write the individual sequences to a new FASTA file
writeXStringSet(extracted_seqs, filepath = output_file)
moref(output_file)






#2. Load the 'devtools' package and install the development version of
#'AnnotationBustR' from GitHub:
library(devtools)
pak::pak("sborstein/AnnotationBustR")  # install the package from GitHub
library(AnnotationBustR)# load the package
IntronExon.example <- AnnotationBust(Accessions=c("KX687911.1", "KX687910.1"), Terms=IntronExonExampleTerms, Prefix="DemoIntronExon")


