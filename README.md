# 16s_practice
A practice workflow of 16s data management and analysis

**This tutorial was taken and adapted from https://benjjneb.github.io/dada2/tutorial.html**


Download the raw reads here:
http://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip

Tutorial copy & paste:
```
#install required programs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.18")

library(dada2); packageVersion("dada2")

#once you have downloaded and unzipped your data, type in the directory here
path <- "~***/MiSeq_SOP"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])

#filter low quality reads
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#trim remaining reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Determine error rate
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#view new filtered read quality
path_2 <- "~/Desktop/16s_practice/MiSeq_SOP/filtered"
fnFs_fixed <- sort(list.files(path_2, pattern="_F_filt.fastq.gz", full.names = TRUE))
fnRs_fixed <- sort(list.files(path_2, pattern="_R_filt.fastq.gz", full.names = TRUE))
plotQualityProfile(fnFs_fixed[1:2])
```

