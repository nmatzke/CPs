#######################################################
# Local version of AnnotationBustR::AnnotationBust
# Ie assume 1 accession, which is a local filename
#######################################################
AnnotationBust_local <- function (Accessions, Terms, Duplicates = NULL, DuplicateInstances = NULL, 
    TranslateSeqs = "None", DuplicateSpecies = FALSE, Prefix = NULL, 
    TidyAccessions = TRUE, Reference = TRUE, Verbose = TRUE, gbfn="insert_filename.gbb") 
		{
		setup='
		wd = "~/GitHub/CPs/ex_cpDNA_genome/NC_023463.1_Pinguicula_ehlersiae_chloroplast/"
setwd(wd)
gbfn <- "NC_023463.1_Pinguicula_ehlersiae_chloroplast_GeSeqAnnotated.gb"
Accessions = "NC_023463.1"
accession.index = 1
Terms = cpDNAterms

Duplicates = NULL
DuplicateInstances = NULL
TranslateSeqs = "Both"
Prefix=NULL
DuplicateSpecies = FALSE
Reference = TRUE
Verbose = TRUE
TidyAccessions = TRUE
		'

	unique.features <- unique(Terms$Feature)
	
	
	# Set up duplicates (not used for a single file)
	if (is.null(Duplicates) == FALSE)
		{
		if (!length(Duplicates) == length(DuplicateInstances))
			{
			stop("Length of Duplicates and DuplicateInstances is not equal. You must specify the number of duplicates for each feature with duplicates you would like to extract")
			}
		new.file.names <- character(length = 0)
		dup.frame <- data.frame(Duplicates, DuplicateInstances)
		singles <- unique.features[unique.features %in% dup.frame$Duplicates == 
				FALSE]
		doubles <- unique.features[unique.features %in% dup.frame$Duplicates == 
				TRUE]
		for (i in 1:length(doubles))
			{
			sub.dups <- subset(dup.frame, dup.frame$Duplicates %in% 
					doubles[i] == TRUE)
			number.names <- paste0(sub.dups$Duplicates, "_", 
					1:sub.dups$DuplicateInstances)
			new.file.names <- append(new.file.names, number.names)
			}
		file.names <- c(as.vector(singles), new.file.names)
		}	else {
		file.names <- as.vector(unique.features)
		} # END if (is.null(Duplicates) == FALSE)
		
	# Set up Accession Table
	Accession.Table <- data.frame(data.frame(matrix(NA, nrow = length(Accessions), ncol=1+length(file.names))))
	colnames(Accession.Table) <- c("Species", file.names)
	if (Reference == TRUE)
		{
		Accession.Table$Reference <- NA
		}

	if (is.null(Prefix))
		{
		File.Prefix <- NULL
		} else {
		File.Prefix <- paste0(Prefix, "_")
		} # END if (is.null(Prefix))


	# Big loop through accessions
	accession.index = 1
	for (accession.index in 1:length(Accessions))
		{
		Current.Accession <- Accessions[accession.index]

		# Install the package if you haven't already
		raw_gb_text <- readr::read_file(gbfn)

		#raw_gb_text <- try(rentrez::entrez_fetch(db = "nuccore", 
		#    id = Current.Accession, rettype = "gbwithparts", 
		#    retmode = "text"))
		gb_lines <- strsplit(raw_gb_text, "\n|\r\n")[[1]]
		org_line_index <- grep("^\\s*ORGANISM", gb_lines)
		org_line <- gb_lines[org_line_index]
		organism_name <- gsub(" ", "_", trimws(sub("^\\s*ORGANISM", "", org_line)))
		if (Reference == TRUE)
			{
			Ref2Store <- AnnotationBustR:::ParseReference(ReadGB = gb_lines)
			Accession.Table[accession.index, "Reference"] <- Ref2Store
			}
		if (Verbose == TRUE)
			{
			message(paste("Working On Accession ", accession.index, " of ", length(Accessions), ": ", Accessions[accession.index], ", ", organism_name, sep = ""))
			}

		Accession.Table[accession.index, "Species"] <- organism_name
		current_record <- AnnotationBustR:::parse_genbank(ReadGB = gb_lines, primary_accession = Current.Accession)
		ifelse(DuplicateSpecies == TRUE, seq.name <- paste(organism_name, Accessions[accession.index], sep = "_"), seq.name <- organism_name)
			
		for (term.index in seq_along(unique.features))
			{
			synonyms <- Terms[Terms$Feature == unique.features[term.index], ]
				for (synonym.index in 1:nrow(synonyms))
					{
					found.type <- base::grep(pattern = paste0("^", synonyms$Type[synonym.index], "$"), x = names(current_record))
					
					if (length(found.type) == 1)
						{
						term.search <- base::grep(pattern = paste0("^", synonyms$Name[synonym.index], "$"), x = current_record[[found.type]]$gene)
						ifelse(length(term.search) > 0, term.search <- term.search, 
							term.search <- base::grep(pattern = paste0("^", 
								synonyms$Name[synonym.index], "$"), x = current_record[[found.type]]$product))
						ifelse(length(term.search) > 0, term.search <- term.search, 
							term.search <- base::grep(pattern = paste0("^", 
								synonyms$Name[synonym.index], "$"), x = current_record[[found.type]]$note))
						if (length(term.search) == 0 && names(current_record)[[found.type]]=="D-loop")
							{
							term.search <- base::grep(pattern = paste0("^", synonyms$Name[synonym.index], "$"), x = current_record[[found.type]]$type)
							}
						if (synonyms[synonym.index, ]$Type %in% c("exon", "intron"))
							{
							temp.term <- current_record[[found.type]][term.search, ]
							IntronExonHits <- which(temp.term$number == synonyms$IntronExonNumber[synonym.index])
							term.search <- term.search[IntronExonHits]
							}
						if (length(term.search) > 0)
							{
							FeatureGrab <- current_record[[found.type]][term.search, ]
							FeatureIndexCheck <- unique(FeatureGrab$feat_index)
							if ((unique.features[term.index] %in% Duplicates) == FALSE)
								{
								target_gr = FeatureGrab[FeatureGrab$feat_index == FeatureIndexCheck[1]]
								Extracted.seq <- AnnotationBustR:::extract_ranges_seq(full_sequence = current_record$sequence, target_gr)
								names(Extracted.seq) <- seq.name
								
								if (synonyms$Type[synonym.index] == "CDS")
									{
									if (TranslateSeqs %in% c("Only", "Both"))
										{
										TF1 = FeatureGrab$feat_index == FeatureIndexCheck[1]
										Current.Trans <- unique(FeatureGrab[TF1]$translation)
										# This produces NA on rpl22; try the 2nd hit
										
										temp_translations = FeatureGrab$translation
										if (sum(is.na(temp_translations)) == length(temp_translations))
											{
											txt = paste0("WARNING in AnnotationBust_local(): No non-NA translations found, continuing. FeatureGrab printed below.")
											cat("\n")
											cat(txt)
											cat("\n")
											cat("term.index: ", term.index, "\n", sep="")
											cat("unique.features[term.index]: ", unique.features[term.index], "\n", sep="")
											cat("found.type: ", found.type, "\n", sep="")
											cat("synonyms$Feature[found.type]: ", synonyms$Feature[found.type], "\n", sep="")
											cat("FeatureGrab:\n")
											print(FeatureGrab)
											warning(txt)
											cat("\n")
											cat("This FeatureGrab will not be recorded. Concluding WARNING.")
											cat("\n")
											
											} else {
											# Correction, if translation is NA, look for another
											if (is.na(Current.Trans) == TRUE)
												{
												translated_TF = !is.na(FeatureGrab$translation)
												translated_TF_nums = (1:length(translated_TF))[translated_TF]
												Current.Trans <- unique(FeatureGrab[translated_TF_nums[1]]$translation)
												} # END if (is.na(Current.Trans) == TRUE)

											if (length(Current.Trans) == 1)
												{
												TransString <- Biostrings::AAStringSet(Current.Trans)
												names(TransString) <- seq.name
												Biostrings::writeXStringSet(x = TransString, filepath = paste0(File.Prefix, unique.features[term.index], "_Translation", ".fasta"), format="fasta", append=TRUE)
												Accession.Table[accession.index, grep(paste0("\\b", unique.features[term.index], "\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
												} # END if (length(Current.Trans) == 1)
											} # END if (sum(is.na(temp_translations)) == length(temp_translations))
										} # END if (TranslateSeqs %in% c("Only", "Both"))
									if (TranslateSeqs %in% c("Both", "None")) {
										Biostrings::writeXStringSet(x = Extracted.seq, 
											filepath = paste0(File.Prefix, unique.features[term.index], 
												".fasta"), format = "fasta", append = T)
										Accession.Table[accession.index, grep(paste0("\\b", 
											unique.features[term.index], "\\b"), 
											colnames(Accession.Table))] <- Accessions[accession.index]
									}
									break
								} else {
								# Continue if (synonyms$Type[synonym.index] == "CDS")
								Biostrings::writeXStringSet(x = Extracted.seq, 
									filepath = paste0(File.Prefix, unique.features[term.index], 
										".fasta"), format = "fasta", append = T)
								Accession.Table[accession.index, grep(paste0("\\b", 
									unique.features[term.index], "\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
								break
								} # END if (synonyms$Type[synonym.index] == "CDS")
							} else {
							# continue: if (unique.features[term.index] %in% Duplicates==FALSE)
							DupID <- which(unique.features[term.index] == 
								Duplicates)
							CurrentDupTargets <- DuplicateInstances[DupID]
							if (CurrentDupTargets > length(unique(FeatureGrab$feat_index)))
								{
								CurrentDupTargets <- length(FeatureGrab)
								warning(paste0("Number of duplicates specified is greater than the number in annotations. Readjusting duplicates for ", 
									unique.features[term.index], " to ", CurrentDupTargets))
								} # END if (CurrentDupTargets > length(unique(FeatureGrab$feat_index)))
							
							
							
							for (DupIndex in 1:CurrentDupTargets)
								{
								CurrentDup <- FeatureGrab[FeatureGrab$feat_index == 
									FeatureIndexCheck[DupIndex], ]
								Extracted.seq <- extract_ranges_seq(full_sequence = current_record$sequence, 
									target_gr = CurrentDup[CurrentDup$feat_index == 
										FeatureIndexCheck[DupIndex]])
								names(Extracted.seq) <- seq.name
								if (synonyms$Type[synonym.index] == "cds")
									{
									if (TranslateSeqs %in% c("Only", "Both"))
										{
										Current.Trans <- unique(CurrentDup$translation)
										if (length(Current.Trans) == 1)
											{
											TransString <- Biostrings::AAStringSet(Current.Trans)
											names(TransString) <- seq.name
											Biostrings::writeXStringSet(x = TransString, 
												filepath = paste0(File.Prefix, 
													unique.features[term.index], "_", 
													DupIndex, "_Translation", ".fasta"), 
												format = "fasta", append = T)
											Accession.Table[accession.index, 
												grep(paste0("\\b", unique.features[term.index], 
													"_", DupIndex, "\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
											} # END if (length(Current.Trans) == 1)
										} # END if (TranslateSeqs %in% c("Only", "Both"))
									if (TranslateSeqs %in% c("Both", "None")) 
										{
										Biostrings::writeXStringSet(x = Extracted.seq, 
											filepath = paste0(File.Prefix, 
												unique.features[term.index], "_", 
												DupIndex, ".fasta"), format = "fasta", 
											append = T)
										Accession.Table[accession.index, 
											grep(paste0("\\b", unique.features[term.index], 
												"_", DupIndex, "\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
										} # END if (TranslateSeqs %in% c("Both", "None")) 
								} else {
								Biostrings::writeXStringSet(x = Extracted.seq, 
									filepath = paste0(File.Prefix, unique.features[term.index], 
										"_", DupIndex, ".fasta"), format = "fasta", 
									append = T)
								Accession.Table[accession.index, grep(paste0("\\b", 
									unique.features[term.index], "_", DupIndex, 
									"\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
								} # END if (synonyms$Type[synonym.index] == "cds")
							} # END for (DupIndex in 1:CurrentDupTargets)
							break
						} # END if (unique.features[term.index] %in% Duplicates==FALSE)
					} # END if (length(term.search) > 0)
				} # END if (length(found.type) == 1)
			} # END for (synonym.index in 1:nrow(synonyms))
		} # END for (term.index in seq_along(unique.features))
	} # END for (accession.index in 1:length(Accessions))


	if (TidyAccessions == TRUE)
		{
		UniqueSpecies <- unique(Accession.Table$Species)
		Final.Accession.Table <- data.frame(data.frame(matrix(NA, 
				nrow = length(UniqueSpecies), ncol = 1 + length(file.names))))
		colnames(Final.Accession.Table) <- c("Species", file.names)
		if (Reference == TRUE)
			{
			Final.Accession.Table$Reference <- NA
			}
		for (species.index in 1:length(UniqueSpecies))
			{
			current.spec <- subset(Accession.Table, Accession.Table$Species==UniqueSpecies[species.index])
			Final.Accession.Table[species.index, 1] <- UniqueSpecies[species.index]
			for (gene.index in 2:length(Accession.Table))
				{
				current.loci <- current.spec[, gene.index]
				found.accessions <- subset(current.loci, !is.na(current.loci))
				numbers <- ifelse(length(found.accessions)==0, NA, paste(found.accessions, sep = ",", collapse = ","))
				Final.Accession.Table[species.index, gene.index] <- numbers
				Sort.Final.Accession.Table <- Final.Accession.Table[order(Final.Accession.Table$Species), ]
				}
			rownames(Sort.Final.Accession.Table) <- 1:nrow(Sort.Final.Accession.Table)
			}
		}	else {
		Final.Accession.Table <- Accession.Table
		Sort.Final.Accession.Table <- Final.Accession.Table[order(Final.Accession.Table$Species), ]
		rownames(Sort.Final.Accession.Table) <- 1:nrow(Sort.Final.Accession.Table)
		} # END if (TidyAccessions == TRUE)
	return(Final.Accession.Table)
	} # END AnnotationBust_local
