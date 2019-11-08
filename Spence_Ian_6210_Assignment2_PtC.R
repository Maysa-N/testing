#EXPLORATION OF POTENTIAL CANDIDATE SPECIES TO RESEARCH 16S AND COI EVOLUTION IN DAPHNIA SPECIES USING DENDROGRAM COMPARISON

#IAN SPENCE - SOFTWARE TOOLS 6210, ASSIGNMENT #2 PT. C.

##      Load Packages#####

library(muscle)
library(DECIPHER)
library(tidyverse)
library(stringr)
library(Biostrings)
library(seqinr)
library(rentrez)
library(ape)
library(dendextend)

##      Biological Question and Introduction#####


#Could the phylogenies of COI and 16S genes in closely related Daphnia species be compared to generate novel hypotheses about the evolutionary process driving the divergence or conservation of these genes?


#DNA barcodes are small sections of DNA that are conserved enough to compare, yet divergent enough to distinguish. Using different genes as DNA barcodes may produce different results in a phlyogenetic analysis, due to the many biological factors affecting the history of each gene.  Understanding evolutionary history can be essential for conserving biodiversity as well as understanding the history of these organisms.  Although COI is commonly used as the barcode of choice for animals, let's assume that we've found ourselves working in a lab that is interested in the evolutionary history of 16S in Daphnia and they're looking for which species to analyze enrich the value of the output of their program.  Let's analyze and compare the phylogenetic trees of Daphnia species created for 16S and COI genes and generate hypotheses based on this comparison to direct further research. 



#NOTE: Daphnia 16S data was originally downloaded for part B of this assignment however all of the code has been provided in this script. 

##PT 1. Download Nucleotide Data for COI and 16S Genes ######


#Load nucleotide data of 16S from NCBI and COI from BOLD.

#First read in the COI Daphnia data from BOLD
Daphnia <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Daphnia&format=tsv")

#Write the file to our computer
write_tsv(Daphnia, "Daphnia.Assignment2C.fasta")

#Read the file into our envrionment and name it "Daphnia"
Daphnia.BOLD <- read_tsv("Daphnia.Assignment2C.fasta")
#Downloaded 7:18pm Oct 28


#Now read in 16S data from NCBI. 

#We want to search the nucleotide database for 16S gene of Daphnia. Also we only want the 16S gene, not entire genomes.

(dap_16S_search <- entrez_search(db = "nucleotide", term = "(Daphnia[ORGN] AND 16S[TITL]) NOT (genome[TITL])", retmax = 328))

#Fetch the file from the online database to our R environment.
dap_16S_fetch <- entrez_fetch(db = "nucleotide", id = dap_16S_search$ids, rettype = "fasta")

#Check the class
class(dap_16S_fetch)

#Check the head of the file
head(dap_16S_fetch)

#Write it to our computer
write(dap_16S_fetch, "Daphnia_16S_fetch.fasta", sep = "\n")
#downloaded at 2:50pm October 28th, 2019

#Read the file to an object called "Dap.StringSet"
Daphnia.NCBI <- readDNAStringSet("Daphnia_16S_fetch.fasta")


##PT 2. Structure Data and Refine Daphnia Species List #####
  

#---------------16S Data Manipulation-----------------#

#We downloaded our 16S file as a DNAStringSet so let's make a dataframe to work on our data.

Dap.DF<- data.frame(Dap_16S_Title = names(Daphnia.NCBI), Dap_16S_Sequence = paste(Daphnia.NCBI))

#We want to change our title to characters and sequence to DNAStringSet.

Dap.DF$Dap_16S_Sequence <- DNAStringSet(Dap.DF$Dap_16S_Sequence)
Dap.DF$Dap_16S_Title <- as.character(Dap.DF$Dap_16S_Title)

#Check if it worked 

class(Dap.DF$Dap_16S_Title)

class(Dap.DF$Dap_16S_Sequence)


#Create new columns with Identifier, Species Name and Gene Name.  This data can be found in the Dap_16S_Title column so we can use the word function to extract the required rows for our new columns. 

Dap.DF$Identifier <- word(Dap.DF$Dap_16S_Title, 1L)
Dap.DF$Species_Name <- word(Dap.DF$Dap_16S_Title, 2L, 3L)
Dap.DF$Gene_Name <- word(Dap.DF$Dap_16S_Title, -7L)


#Rearrange the columns.

Dap.DF <- Dap.DF[, c("Identifier", "Species_Name", "Gene_Name", "Dap_16S_Sequence", "Dap_16S_Title")]

#Let's see which gene names are not labeled "16S". 
(No_16S <- c(grep("16S", Dap.DF$Gene_Name, invert = TRUE)))

#Make a dataframe of the non matches to take closer look
Dap.DF_No16S <- Dap.DF[c(No_16S), ]


#Search the sequence titles in the dataframe for "16S"
grep("16S", Dap.DF_No16S$Dap_16S_Title)

# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16

#We can see that all 16 variables have 16S in their sequence tag but we can see the top two samples (17 & 18) are not the same type of 16S sequences. We can read that sample 17 is cDNA 5', and sample 18 is BM4 gene.  

Dap.DF_No16S

#Double check it is 17 and 18 that we want to remove.
Dap.DF[17:18, 5]

#Remove them from the analysis.
Dap.16S <- Dap.DF[-c(17,18),]

#The rest of the samples are 16S gene sequences so let's overwrite the gene name column. 
Dap.16S$Gene_Name <- "16S"


#Let's do the same processes that we did for gene name for species name and search for any observations that don't have Daphnia as genus. 

No_Dap <- c(grep("Daphnia", Dap.16S$Species_Name, invert = TRUE))

No_Dap
#277 291 292 293 294 295 296 297 298 299

#These samples no not have Daphnia in their species name.

#Make a dataframe of the non matches to look closer.
Dap.DF_NoDap <- Dap.16S[c(277, 291, 292, 293, 294, 295, 296, 297, 298, 299), ]
Dap.DF_NoDap

#We can see that we have ten observations and each of them are from the Daphniopsis genus.  Let's remove these Daphniopsis samples from our dataset as our analysis is only interested in sequences from the Daphnia genus.

No_Dap
#277 291 292 293 294 295 296 297 298 299

Dap.16S.Clean <- Dap.16S[-c(277, 291, 292, 293, 294, 295, 296, 297, 298, 299),]



#Finally, let's remove the "Daphnia.sp" samples as we want our results to depict evolutionary relationships between known Daphnia species. We'll need to change the sequence data back to characters to filter. 

class(Dap.16S.Clean$Dap_16S_Sequence)
Dap.16S.Clean$Dap_16S_Sequence <- as.character(Dap.16S.Clean$Dap_16S_Sequence)
class(Dap.16S.Clean$Dap_16S_Sequence)

#Filter out all samples with the species name "Daphnia.sp"

Dap.16S.Clean.Final <- Dap.16S.Clean %>%
  filter(!Dap.16S.Clean$Species_Name == "Daphnia sp.")


#---------------COI Data Manipulation-----------------#


#Now we'll structure COI genes from BOLD. 

#Filter for all Daphnia samples sequenced with the COI markercode that have nucleotides and a species name.

Daphnia.COI <- Daphnia.BOLD %>%
  filter(!is.na(nucleotides)) %>%
  filter(markercode == "COI-5P") %>%
  filter(!is.na(species_name))


#Make a list of species that exist in COI and 16S.

#Make a table into a dataframe with all of the species name from COI data
Species.COI.counts <- as.data.frame(table(Daphnia.COI$species_name))

#Make a table into a dataframe with all of the species name from 16S data
Species.16S.counts <- as.data.frame(table(Dap.16S.Clean.Final$Species_Name))

#Which samples do not exist in both datasets? 
(Species.match <- Species.16S.counts$Var1 %in% Species.COI.counts$Var1)

#6, 35, 37 are FALSE meaning they are in one set but not the other, lets check what they are. 

check.species <- c(6, 35, 37)

Species.16S.counts$Var1[check.species]

#Daphnia cristata  Daphnia sinensis  Daphnia turbinata

#A quick search in the script revealed no 'Daphnia cristata' or 'Daphnia turbinata' matches in the COI data but there was a 'Daphnia sinensis' under a similar name 'Daphnia cf. sinensis'.  This means we will remove number 6, 'Daphnia cristata' and number 37 'Daphnia turbinata' and we will rename 35 from 'Daphnia cf. sinensis' to 'Daphnia sinensis' in the COI dataset.

#These are the species that we will use (All the species in 16S except for the two that don't have COI equivelents #6 and #37.)
Matched.Species <- Species.16S.counts[-c(6,37),]


#Which rows need their species names need fixed?
which(Daphnia.COI$species_name == "Daphnia cf. sinensis")
#179 1358

Daphnia.COI[179, "species_name"] <- "Daphnia sinensis"

#Check if that switched it correctly,
which(Daphnia.COI$species_name == "Daphnia cf. sinensis")

#Great, lets do the other sample.

Daphnia.COI[1358, "species_name"] <- "Daphnia sinensis"

#Check if that switched it correctly,
which(Daphnia.COI$species_name == "Daphnia cf. sinensis")


#This is our refined species list

Species.List <- as.vector(Matched.Species$Var1)


#Now that we have our list, lets create a dataframe with both genes for each of the species that exist in both datasets.

#First we'll refine the COI section portion of the data. Sample 1 observation from each species.

set.seed(33)
Daphnia.COI.Selected <- Daphnia.COI %>%
  group_by(species_name) %>%
  sample_n(1)


#Filter out samples that aren't represented in our character vector "Species.List"
Daphnia.COI.Sampled <- Daphnia.COI.Selected %>%
  filter(grepl(paste(Species.List, collapse="|"), species_name))


dim(Daphnia.COI.Sampled)

#We can see that our sample has 40 even though our list has 36.  After checking, we can see that some species have duplicates where one of the with variants has an extension.  We can create a logical vector of the species that match our list "Species.List".

Daphnia.COI.Sampled$species_name %in% Species.List

#[1]  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#[22]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE

#Remove samples 4, 30, 31, and 36.

Filtered.COI <- Daphnia.COI.Sampled[-c(4,30,31,36),]

#Are all matches met with only 36 observations?
Filtered.COI$species_name %in% Species.List
dim(Filtered.COI)

#Yes, great. 


#Now let's sample each of the 36 species with 16S genes.

       
set.seed(33)
Dap.16S.Selected <- Dap.16S.Clean %>%
  group_by(Species_Name) %>%
  sample_n(1)

dim(Dap.16S.Selected)

#This gives, us 39 observations, we need to remove the other two samples that are not found for COI 'Daphnia cristata' and 'Daphnia turbinata'.


Dap.16S.Selected$Species_Name %in% Species.List

#6 and #37 do not match the list.

Filtered.16S <- Dap.16S.Selected[-c(6,36, 38),]

#Are all matches met with only 36 observations?
dim(Filtered.16S)
Filtered.16S$Species_Name %in% Species.List

#--------------------COMBINING DATASETS----------------------#

#In our data we only need a few columns.  We'd like species name, gene name, ID, and sequence.

#First isolate the COI columns we want by specifying column names in a vector

columns.COI <- c("species_name", "markercode", "sampleid", "nucleotides")
columns.16S <- c("Species_Name", "Gene_Name", "Identifier", "Dap_16S_Sequence")

#Create a dataframe with the four columns. 

Filtered.COI.4columns <- Filtered.COI[,columns.COI]
Filtered.16S.4columns <- Filtered.16S[,columns.16S]


#Change the column names to the following "Species", "Gene", "ID", and "Sequence".

column.names <- c("Species", "Gene", "ID", "Sequence")
colnames(Filtered.COI.4columns) <- column.names
colnames(Filtered.16S.4columns) <- column.names


#Now lets merge the two 4 frame datasets, they should line up!


Merged.DF <- merge(Filtered.COI.4columns,Filtered.16S.4columns, by = "Species")

##When combining dataframes with the same columns either .x or .y get added to the column name.  To make it easier to recall lets change the column names to replace the .x or .y with the respective gene name. .x is COI and .y is 16S.

#Because we are adding gene names to the column titles, we can remove the gene columns. 

Merged.DF.5col <- Merged.DF[,-c(2,5)]

#Change the column names.

merged.column.names <- c("Species", "ID.COI", "Sequence.COI", "ID.16S", "Sequence.16S")
colnames(Merged.DF.5col) <- merged.column.names

Merged.DF.5col$Species

#We can see that we have only Daphnia as our genus so we can replace "Daphnia" with "D." to improve our visualization later on.

Merged.DF.5col$Species <- str_replace(Merged.DF.5col$Species, "Daphnia", "D.\\")
Merged.DF.5col$Species

#Now let's do an alignment for the COI and 16S sequences.


##PT 3. Run Alignments for COI and 16S Genes #####

#Check the str characteristics of COI and 16S sequences

mean(str_length(Merged.DF.5col$Sequence.COI))
min(str_length(Merged.DF.5col$Sequence.COI))
max(str_length(Merged.DF.5col$Sequence.COI))

mean(str_length(Merged.DF.5col$Sequence.16S))
min(str_length(Merged.DF.5col$Sequence.16S))
max(str_length(Merged.DF.5col$Sequence.16S))

#Check the classes of the sequences

class(Merged.DF.5col$Sequence.COI)
class(Merged.DF.5col$Sequence.16S)

#Change the sequence data from character to DNAStringSet. 

Merged.DF.5col$Sequence.COI <- DNAStringSet(Merged.DF.5col$Sequence.COI)
Merged.DF.5col$Sequence.16S <- DNAStringSet(Merged.DF.5col$Sequence.16S)

class(Merged.DF.5col$Sequence.COI)
class(Merged.DF.5col$Sequence.16S)


#Now let's do a quick alignment to visualize our 16S data. 

Seq.16S.Align <- DNAStringSet(muscle::muscle(Merged.DF.5col$Sequence.16S, maxiters = 2))

#Let's open it up
BrowseSeqs(Seq.16S.Align)


#And do the same for COI

Seq.COI.Align <- DNAStringSet(muscle::muscle(Merged.DF.5col$Sequence.COI, maxiters = 2))

#Let's open it up
BrowseSeqs(Seq.COI.Align)

#Two sequences, 11 and 28, are much longer than the rest.  

#Print their sequence and do a quick BLAST search. 

(Test.sequence <- as.data.frame(Merged.DF.5col[11,3]))

#TTGCGACGTTGACTCTTTTCAACTAATCATAAGGATATTGGAACCCTTTACTTCATTTTTGGAATTTGGTCCGGGATAGTCGGAACCGCTCTTAGTTTACTGATCCGGGCTGAACTTGGACAATCAGGAAGATTAATTGGAGATGACCAAATTTACAATGTAATTGTAACTGCCCACGCTTTTGTAATAATTTTTTTCATAGTTATGCCTATCATGATTGGGGGTTTTGGGAACTGACTAGTGCCCTTAATGCTTGGGGCTCCTGATATAGCTTTTCCCCGATTAAATAATCTTAGCTTCTGACTTTTACCCCCTGCTCTTACTCTCTTACTGGTAGGAGGGGCGGTGGAAAGTGGGGCAGGAACTGGGTGAACAGTTTATCCCCCCCTCTCTGCAGGGATTGCCCATGCCGGGGCCTCAGTAGATCTCAGAATTTTTTCTCTCCATCTTGCAGGAATCTCATCAATTCTTGGGGCAGTTAATTTTATTACTACTATCATTAATATACGGTCCTTAGGCATGACCTTGGATCGGATCCCTTTATTTGTATGAGCAGTCGGAATCACAGCACTCCTACTTCTCCTTAGCCTCCCCGTTCTAGCGGGGGCTATTACTATACTTCTGACCGACCGTAACCTTAATACTTCATTTTTTGACCCGGCGGGAGGGGGGGATCCAATTTTATACCAACACCTATTTTGGTTTTTTGGCCACCCGGAAGTCTATATCCTCATTCTTCCTGGATTTGGGATAATTTCCCACATTATTAGACACGAAAGAGGTAAAAAGGAAGCTTTCGGGACTTTAGGTATAATTTACGCCATACTGGCCATTGGGATCCTAGGGTTCGTAGTATGGGCCCACCACATGTTTACAGTAGGAATAGACGTAGATACCCGAGCTTACTTTACTGCCGCTACTATAATTATCGCAGTACCTACAGGTATTAAAATCTTTAGTTGACTAGGGACTCTTCATGGGACTCAATTAGTGTTTACCCCCTCTTTATTATGAGCGGCAGGATTTGTTTTCTTGTTTACAGTAGGAGGCCTAACGGGAGTTGTTTTAGCTAATTCGAGAATCGATATTATCCTTCACGATACCTACTATGTTGTTGCTCACTTTCACTATGTTCTATCTATAGGGGCAGTTTTCGCAATTTTTGCAGGAGTGGCCCACTGATATCCACTATTTACAGGACTAACTCTTCATGCTCAATGGCTTAAAATCCATTTCTTCACAATATTTCTTGGGGTTAATCTCACCTTCTTCCCTCAGCATTTTTTAGGATTAGCTGGCATACCTCGACGCTATTCTGACTACCCTGACGCCTACCTGGCATGAAATGTAGTATCATCAGTGGGCTCTATAATCTCTTTTGTAGCAACATTAGGCTTTATTTTTATTATTTGGGAGTCTCTTATTAGACTGCGGCCTGCTCTATTTGCTACTCATTTATCTACCTCCATTGAATGGCAGCATTCTTTCCCGCCTGCGGAGCATAGCTACAATGAGTTAGTTTATATCTCGCAATTC

#Daphnia galeata mitochondrial COI gene for cytochrome c oxidase subunit I, partial cds, isolate: Dgal_01


(Test.sequence2 <- as.data.frame(Merged.DF.5col[28,3]))

#GCCGAACTGGGGCAATCAGGAAGTTTAATCGGAGACGATCAGATTTATAATGTAATTGTTACCGCCCACGCGTTCGTAATAATTTTTTTTATAGTTATGCCCATCATAATTGGGGGGTTCGGGAACTGATTAGTCCCCCTTATACTAGGGGCCCCGGATATGGCCTTCCCTCGACTTAATAATTTAAGATTCTGATTTCTTCCTCCGGCACTTACACTACTTTTAGTTGGGGGGGCAGTAGAAAGAGGGGCTGGGACAGGGTGAACCGTGTACCCCCCACTCTCGGCTGGGATTGCTCATGCGGGGGCCTCTGTAGATCTAAGAATCTTCTCCCTCCATTTGGCGGGTATTTCTTCTATTCTTGGGGCTGTCAACTTTATTACTACTATTATCAATATACGATCGGCTGGGATATCATTAGACCGAATTCCTTTATTTGTATGAGCAGTAGGCATTACTGCGCTACTTCTACTATTGAGCCTACCAGTCCTAGCAGGGGCAATTACTATGCTTCTTACTGATCGTAACCTTAATACCTCCTTTTTTGACCCTGCAGGAGGGGGAGATCCGATCTTATACCAGCATTTGTTTTGGTTCTTTGGTCATCCTGAAGTATACATTTTAATTTTACCCGGGTTTGGTATAATTTCTCATATTATTAGCCACGAAAGAGGTAAAAAGGAAGCATTCGGCACGCTTGGGATAATTTACGCTATACTAGCTATCGGGGTGTTAGGCTTCGTAGTATGGGCCCACCATATATTTACTGTGGGCATGGATGTTGACACACGAGCGTATTTTACTGCAGCTACAATGATTATTGCGGTTCCCACAGGAATTAAGATCTTCAGATGACTTGGTACACTTCACGGAACCCAGCTTGTATTTACCCCTTCTCTCTTGTGAGCGGCAGGGTTTGTGTTTCTGTTTACAGTGGGGGGGCTTACTGGGGTAGTTTTAGCCAACTCAAGAATTGACATTATTCTTCATGACACGTACT

#"Daphnia pulex isolate IN-07 cytochrome c oxidase subunit I (cox1) gene, partial cds; mitochondrial"

#Looks like both of these sequences are of the correct type, just longer. Because all 16S sequences have a similar length, to minimize the noise in the comparison between these two genes, the two long COI genes will be removed. 

Merged.DF.5col.size <- Merged.DF.5col[-c(11,28),]

#Rename the rows as 1:34 (from 1:36) in our merged dataframe.

row.names(Merged.DF.5col.size) <- Merged.DF.5col.size$Species

#Check the alignments now

Seq.COI.Align <- DNAStringSet(muscle::muscle(Merged.DF.5col.size$Sequence.COI, maxiters = 2))
BrowseSeqs(Seq.COI.Align)

Seq.16S.Align <- DNAStringSet(muscle::muscle(Merged.DF.5col.size$Sequence.16S, maxiters = 2))
BrowseSeqs(Seq.COI.Align)

#Looks good, lets adjust the iterationns and gap penalty

#Run the COI alignment using muscle 
Seq.COI.Align <- DNAStringSet(muscle::muscle(Merged.DF.5col.size$Sequence.COI, maxiters = 5, gapopen = -50))
BrowseSeqs(Seq.COI.Align)

#Run the 16S alignment using muscle 
Seq.16S.Align <- DNAStringSet(muscle::muscle(Merged.DF.5col.size$Sequence.16S, maxiters = 5, gapopen = -50))
BrowseSeqs(Seq.16S.Align)

##PT 4. Create and Compare Dendrograms with Tanglegram ####

#We need to turn our alignment data into DNAbin objects
class(Seq.COI.Align)
class(Seq.16S.Align)

dnaBin.COI <- as.DNAbin(Seq.COI.Align)
dnaBin.16S <- as.DNAbin(Seq.16S.Align)

class(dnaBin.COI)
class(dnaBin.16S)


#Create a distance matrix for our alignments. The TN93 model considers evolutionary influence of transitions and transversions.

DistanceMatrix.COI <- dist.dna(dnaBin.COI, model = "TN93")
DistanceMatrix.16S <- dist.dna(dnaBin.16S, model = "TN93")


#Use our distance matrix's for each gene, create hierarchical clusters with complete linkage.
dclust1 <- hclust(DistanceMatrix.COI, method = "complete")
dendro1 <- as.dendrogram(dclust1)

dclust2 <- hclust(DistanceMatrix.16S, method = "complete")
dendro2 <- as.dendrogram(dclust2)

#Make a dendrogram list with dendrograms of both genes.
dl <- dendlist(dendro1, dendro2)


#Create a tanglegram, add titles, and set the color and highlight features. 
tanglegram(dl, sort = FALSE, main_left = "COI Dendrogram", main_right = "16S Dendrogram", common_subtrees_color_lines = TRUE, highlight_distinct_edges  = TRUE, highlight_branches_lwd = FALSE)


#Without further literature review to optimize alignment parameters, adjusting the gap penalty and the clustering method till a suitable comparison was created proved to be an effective way to infer phylogenetic relationships between genes.

#The structure of many groups of species look similar between the dendrograms. The pink lines connect 27 (D. thomsoni) and 4 (D. nivalis). The blue lines connect species 20 (D. ) and 13 (D. ). The brown lines connect species 21 and 5.  The species connected with the brown lines 




matched.pairs <- c(27, 4, 20, 13, 21, 5)


Merged.DF.5col.size[matched.pairs, 1]

#Pink matched pair
# "D. reflexa" (27)    
# "D. cephalata" (4)

#Blue matched pair
# "D. muddensis" (20)
# "D. longicephala" (13)

#Brown matched pair
# "D. neocitrina" (21)
# "D. citrina" (5)

#Of the matches, we can compare estimated distance ratios and differences between the 16S and COI genes.

  #The pink line pair:
pink.COI.dist <- 0.08
pink.16S.dist <- 0.001
pink.COI.dist/pink.16S.dist
#Pink COI:16S distance ratio = 80:1

#The blue line pair:
blue.COI.dist <- 0.05
blue.16S.dist <- 0.01
blue.COI.dist/blue.16S.dist
#Blue COI:16S distance ratio = 5:1

#The brown line pair:
brown.COI.dist <- 0.14
brown.16S.dist <- 0.01
brown.COI.dist/brown.16S.dist
#Brown COI:16S distance ratio = 14:1


#Matched pair COI:16S distance ratios:

#Brown = 14 = (0.14 / 0.01)
#Blue = 5 = (0.05 / 0.01)
#Pink = 80 = (0.08 / 0.001)


#Max and Min distance of any gene in a matched pair 

#Max = 0.14 (brown, COI)
#Min = 0.001 (pink, 16S)




##      Conclusion and References #####


#Conclusion

#This data can be used to guide the generation of hypotheses and develop statistical tests to measure if these values are likely to indicate anything of further interest. By comparing high and low levels of conservation between the COI and 16S genes between species with common ancestors,  we may construct a more informed model to predict the factors driving evolution of COI and 16S genes in Daphnia.  Further research may also show that the levels of variation in conservation between the COI and 16S genes in matched pairs are a result of more recent speciation events and comparing relationships between COI and 16S phylogenies won't yield any new information. If a common ancestor was indeed shared by the 'matched pairs' marked by pink, blue, and brown lines in the tanglegram, we can compare the rate of evolution of these genes in each of the three pairs to determine the conditions and factors which drive the evolution of these genes.  What factors were present to cause 0.14 distance in COI vs only 0.01 in 16S of brown pair? Why did the COI of the pink pair diverge 0.08 while the 16S had a distance of 0? Why is the divergence of the brown COI almost 3x higher than the blue COI but the 16S are nearly the same?  With this exploratory analysis we have provided a direction for the further research of the evolution of COI and 16S in Daphnia and the Animal Kingdom.



#-----------------References-------------------#

#R Core Team. (2019). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

#National Center for Biotechnology Information (NCBI)[Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – [cited 2019, Oct 28]. Available from: https://www.ncbi.nlm.nih.gov/

#Ratnasingham, S., & Hebert, P. D. (2007). bold: The Barcode of Life Data System (http://www.barcodinglife.org). Molecular ecology notes, 7(3), 355–364. doi:10.1111/j.1471-8286.2007.01678.x

#https://stackoverflow.com/questions/38724690/r-filter-rows-that-contain-a-string-from-a-vector/38726850

#https://www.r-bloggers.com/dendextend-a-package-for-visualizing-adjusting-and-comparing-dendrograms-based-on-a-paper-from-bioinformatics/

  #Websites acessed October 28, 2019.
 