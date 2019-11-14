      #Assignment_2_part_C
      #In this experiment, Champiaceae genus in the plant kingdom was used to identify whether the two genes (COI-5P and rbcL) yield the same or different phylogenetic hypotheses. Champianceae is a red algae which is ecologically significant as a primary producer. This algae provides structural habitat for marine organisms and also important for the maintanance of coral reefs. Also important providers of food and pharmaceutical substances. Phylogeny reveals the pattern of evolutionary relationships, recognizing the historical pattern of speciation and divergence which leads to classify life according to evolutionary scheme.Using the genus Champiaceae, identify whether the two genes (COI-5P and rbcL) yield the same or different phylogenetic hypotheses
      
      #Load packages
      library(permute)
      library(lattice)
      library(vegan)
      library(Biostrings)
      library(ape)
      library(tidyverse)
      library(RSQLite)
      library(BiocManager)
      library(DECIPHER)
      library(muscle)
      library(purrr)
      library(maps)
      library(mapdata) 
      library(dendextend)
      library(seqinr)
      library(IRanges)
      library(S4Vectors)
      library(stats)
      library(gplots)
      library(d3heatmap)
      library(circlize)
      library(corrplot)
      
      #Champia data file was obtained 20th October 2019 using the following line of code,
      Champia <- read_tsv ("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Champia&format=tsv")
      #Checking the markercode types with their counts
      table (Champia$markercode)
      #Create Champia data set without record missing data in marker code, nucleotides and country
      Champia_NArm <- Champia  %>%
        filter (!is.na(markercode)) %>%
        filter (!is.na(nucleotides)) %>%
        filter (!is.na(country))
      #For a better visualization  use ggplot to map all the marker genes of Champia genus according to the country. The data set was set in the argument data. aes means aesthetic mappings that describe how variables in the data are mapped to visual properties of geoms. stat gives the count of records. The function labs used for labels (to display the title, x axis and y axis). Flip cartesian coordinates using coord_flip, so that vertical becomes horizontal.
      ggplot(data = Champia_NArm) + 
        geom_bar(mapping = aes(x = country, colour = markercode, fill = markercode), stat = "count") +
        labs(title = "Champia BOLD Records per Country and Marker", x = "Countries", y = "Number of Records") +
        coord_flip()  
      
        #Map implies high number of species in marker code COI-5P and rbcL compared to other genes. Therefore, COI-5P and rbcL genes will be used for further downstream analysis
      #Filtering COI-5P sequences having nucleotide data and species names and assigning to new subset called Champia_COI5P. Any remaining records missing nucleotides and species names are excluded by using "!is.na".
      Champia_COI5P <- Champia %>% 
        filter(markercode=="COI-5P") %>% 
        filter(!is.na(nucleotides)) %>% 
        filter(!is.na(species_name)) 
      #Filtering rbcL sequences having nucleotide data and species names and assigning to new subset called Champia_rbcL.
      Champia_rbcL <- Champia %>% 
        filter(markercode=="rbcL") %>% 
        filter(!is.na(nucleotides)) %>% 
        filter(!is.na(species_name))
      #Checking new data sets
      dim(Champia_COI5P) # this gives a subset with 63 rows and 80 columns
      dim(Champia_rbcL) #this gives a subset with 45 rows and 80 columns
      #Check whether there is any missing value in the columns "nucleotides" and "species name"
      sum(is.na(Champia_COI5P$nucleotides))
      sum(is.na(Champia_COI5P$species_name))
      sum(is.na(Champia_rbcL$nucleotides))
      sum(is.na(Champia_rbcL$species_name))
      
      #str_count in stringr package#Count the lengths of nucleotides in Champia_COI5P and rbcLto get an idea about the lengths and to check whether it contains unusual lengths
      str_count(Champia_COI5P$nucleotides)
      str_count(Champia_rbcL$nucleotides)
      #for a nice and easy visualization use plot histogram on nucleotide lengths. rename the x axis as nucleotide lengths and the title of the histogram in the argument "main".
      hist(str_length(Champia_COI5P$nucleotides),col="violet",xlab = "COI-5P nucleotide length",main = "Histogram of COI-5P nucleotide length")
          hist(str_length(Champia_rbcL$nucleotides),col="cyan",xlab = "rbcL nucleotide length",main = "Histogram of  rbcL nucleotide length")
# first edit
# instead of check the length and then use the histogram to check unusual length we can determine the length that we need for both genes:
table(str_length(Champia_COI5P$nucleotides))
seq_length <- Champia_COI5P %>%
  filter(str_length(nucleotides) > 600)
# for rbcl gene:
table(str_length(Champia_rbcL$nucleotides))
seq_length2 <- Champia_rbcL %>%
  filter(str_length(nucleotides) > 1300)
# Now we are sure that we have length in the range we need.)
# second edid :
ggplot give better view:
ggplot(data = Champia_COI5P, aes (x = str_length(Champia_COI5P$nucleotides))) + geom_histogram(color = "black", fill = "blue", binwidth = 10) + labs(x= "COI-5P nucleotide length", y = "count", title = "Histogram of COI-5P nucleotide length")
ggplot(data = Champia_rbcL, aes (x = str_length(Champia_rbcL$nucleotides))) + geom_histogram(color = "green", fill = "red", binwidth = 30) + labs(x= "rbcl nucleotide length", y = "count", title = "Histogram of rbcl nucleotide length")
        
#According to the lengths of nucleotides there are no unusual data, therefore can proceed with this data set for the next steps
      #Combining the rows that appear in either or both Champia_COI5P and Champia_rbcL to get the count of COI5P and rbcL markers according to the country#This will use when plotting country wise gene count on the selected two markers
      dfAll_COI_rbcL <- union(Champia_COI5P,Champia_rbcL)
      #map the records for COI-5P and rbcL gene of Champia for better visualization. The data set was set in the argument data. aes means aesthetic mappings that describe how variables in the data are mapped to visual properties of geoms. stat gives the count of records. The function labs used for labels (to display the title, x axis and y axis). Flip cartesian coordinates using coord_flip, so that vertical bars becomes horizontal.
      ggplot(data = dfAll_COI_rbcL) + 
        geom_bar(mapping = aes(x = country, colour = markercode, fill = markercode), stat = "count") +
        labs(title = "Champia COI-5P and rbcL BOLD Records per Country and Marker", x = "Countries", y = "Number of Records") +
        coord_flip() 
        #change the format of the nucleotide data to be suitable for downstream analysis. Here, the sequence data were changed to a DNAStringSet object class, using a function from the Biostrings package.
      Champia_COI5P$nucleotides <- DNAStringSet(Champia_COI5P$nucleotides)
      Champia_rbcL$nucleotides <- DNAStringSet(Champia_rbcL$nucleotides)
      #MUSCLE is more theoretically reliable.It looks better for translations. It is significantly faster than the other methods  even   for large dataset. Also MUSCLE, has several different alignement algorithms to choose from (for different types of data) and has many useful extra options to work with alignments#To find out what default settings are being used for all the arguments, a log file can be written to disk using the log argument in conjunction with the verbose argument, log = "log.txt", verbose =T. This will write out the default values to the file log.txt in the current working directory of R. (verbose =T generates outputs).
      Alignment_Champia_COI5P <- DNAStringSet(muscle::muscle(Champia_COI5P$nucleotides , log = "log.tx", verbose = T))
      Alignment_Champia_rbcL <- DNAStringSet(muscle::muscle(Champia_rbcL$nucleotides , log = "log.tx", verbose = T))
      #BrowseSeqs function in DECIPHER package opens an html file in a web browser to show the sequences in a DNAStringSet.#BrowseSeqs converts the DNAStringSet into html format for viewing in a web browser. The sequences are colored and easy to recognize sequences with unusual nucleotide patterns.
      BrowseSeqs(Alignment_Champia_COI5P)
          BrowseSeqs(Alignment_Champia_rbcL)
      #Check the number of gaps across all sequences in the alignment to see if there's any unusual count. The str_count argument for the function lapply, placed after a comma, within the arguments to be passed to lapply.
      lapply(Alignment_Champia_COI5P, str_count, "-")
      lapply(Alignment_Champia_rbcL, str_count, "-")
      #as.DNAbin function(ape package) can convert the class from the package Biostrings in BioConductor for consistency within ape package. Converting this data type is essential for further downstream analysis.
      dnaBin_Alignment_Champia_COI5P <- as.DNAbin(Alignment_Champia_COI5P)  
      dnaBin_Alignment_Champia_rbcL <- as.DNAbin(Alignment_Champia_rbcL)
      #dist.dna function in ape package, computes a matrix of pairwise distances from DNA sequences using a model of DNA evolution. Tamura and Nei model was used to calculate the distances it  allows for different rates of transitions and transversions, heterogeneous base frequencies, and between site variation of the substitution rate.
      distanceMatrix_COI5P <- dist.dna(dnaBin_Alignment_Champia_COI5P, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
      class(distanceMatrix_COI5P)# check the data type of the distance matix
      distanceMatrix_rbcL<- dist.dna(dnaBin_Alignment_Champia_rbcL, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
      class(distanceMatrix_rbcL)
      # IdClusters() function from the package DECIPHER, cluster according to single linkage method and using a 0.02 (2%) divergence threshold. 0.02 is the cutoff that separating the sequences in the same cluster. The method single gave a low scale bar compared to other methods and it is good to cluster nearest neighbours similarity. To visualize a tree here "dendrogram" was used as the type      
      clusters_COI5P <- IdClusters(distanceMatrix_COI5P,
                                   method = "single",
                                   cutoff= 0.02,
                                   showPlot = TRUE,
                                   type = "dendrogram",
                                   verbose = TRUE)
      title("Dendrogram of Champia on COI-5P gene ") #Display the title of Dendrogram
        #According to the dendrogram the 7th sequence seems to   be an outlier, To check whether the sequence is an outlier, a nucleotide frquency test and a BLAST search will be performed. Because there should be a good reason to exclude unusual data from the data set otherwise this exclusion may result in missing data which leads to wrong result. 
      Champia_COI_NucFreq <- as.data.frame(letterFrequency(Champia_COI5P$nucleotides, letters = c("A", "C", "G", "T"))) %>%
        mutate(ATproportion = ((A + T) / (A + T + G + C))) %>%
        bind_cols(processid = Champia_COI5P$processid)
      
      #AT propotion of the 7th sequence in Champia_COI_NucFreq does not imply as an unusual result. Furthermore, perform a BLAST search to clarify this result.
      #get the genbank accesion number to Blast the sequence
      Champia_COI5P[7,"genbank_accession"]
      #MG894150 is the genbank accession number for the 7th sequence. Therefore using this accession number performed a BLAST search. The sequence has 100% similarity to COI gene of Champia pseudoparvula species. Therefore this is not an outlier. This data should be included in further downstream analysis.
      #Using distance matrix, cluster the rbcL gene according to single linkage method with a 0.02 (2%) divergence threshold.
      clusters_rbcL <- IdClusters(distanceMatrix_rbcL,
                                  method = "single",
                                  cutoff= 0.02,
                                  showPlot = TRUE,
                                  type = "dendrogram",
                                  verbose = TRUE)
      title("Dendrogram of Champia on rbcL gene ") #Display the title of Dendrogram
      #According to this dendrogram there's no outliers.
      
      #Now let's point out the geographic regions of Champia species on the genes COI-5P and rbcL in the world map
      #The function map uses to display the the geographic regions of Champia species. worldHires draws a world map filling in light green color. Geographic regions of Champia species on genes COI-5P and rbcL genes,
      map("worldHires", col="light green", fill = TRUE)
      
      #Distribution of Champia species on COI-5P gene in the world map according to lattitudes and longitudes. pch is an argument to select a symbol as points.cex is for expansion of symbol. col is to give a color for the symbol. 
      points(Champia_COI5P$lon,
             Champia_COI5P$lat,
             pch=12,
             col="red",
             cex=1.5)
            #red color squares show the distribution of Champia species on COI_5P gene and blue color circles represent the distribution of species on rbcL gene
      #Spread of Champia species on rbcL gene in the world map
      points(Champia_rbcL$lon,
             Champia_rbcL$lat,
             pch=21,
             col="blue",
             cex=1.5)
      title("Distibution of Champia species on COI-5P gene and rbcL gene")
      
      #This map shows the distribution of Champia species with COI-5P and rbcL genes. Champia species with COI-5P and rbcL genes are destributed in several geographic regions. 
      
      #Previously dendrograms were built with all the sequences in each gene to see whether there is any outlier. With Dendrograms, Blast results and nucleotide frequency results, it is obvious that there aren't any outlier. 
      #Now using Champia_COI5P and Champia_rbcL subsets, randomly select sequences per species for COI-5P and rbcL genes to identify whether the two genes (COI-5P and rbcL) yield the same or different phylogenetic hypotheses
      dfChampia_COI5P_Subset <- Champia_COI5P%>% 
        group_by(species_name) %>% #group by species name
        sample_n(1) #Handy function in dplyr for randomly selecting rows. Here 1 means get 1 sample from each species
      
      dfChampia_rbcL_Subset <- Champia_rbcL %>% 
        group_by(species_name) %>% 
        sample_n(1)
      
        #change the format of the nucleotide data to be suitable for downstream analysis. Here, the sequence data were changed to a DNAStringSet object class, using a function from the Biostrings package.
      dfChampia_COI5P_Subset$nucleotides <- DNAStringSet(dfChampia_COI5P_Subset$nucleotides)
      dfChampia_rbcL_Subset$nucleotides <- DNAStringSet(dfChampia_rbcL_Subset$nucleotides)
      
      #names function uses to set the names of nucleotides with the species names, which gives species name as the tip labels  in phylogenetic tree for downstream analysis
      names(dfChampia_COI5P_Subset$nucleotides) <- dfChampia_COI5P_Subset$species_name
      names(dfChampia_rbcL_Subset$nucleotides) <- dfChampia_rbcL_Subset$species_name
      
      #MUSCLE is more theoretically reliable.It looks better for translations. It is significantly faster than the other methods  even for large dataset. Also MUSCLE, has several different alignement algorithms to choose from (for different types of data) and has many useful extra options to work with alignments#To find out what default settings are being used for all the arguments, a log file can be written to disk using the log argument in conjunction with the verbose argument, log = "log.txt", verbose =T. This will write out the default values to the file log.txt in the current working directory of R. (verbose =T generates outputs).
      Alignment_Champia_COI5P_1 <- DNAStringSet(muscle::muscle(dfChampia_COI5P_Subset$nucleotides , log = "log.tx", verbose = T))
      Alignment_Champia_rbcL_1 <- DNAStringSet(muscle::muscle(dfChampia_rbcL_Subset$nucleotides , log = "log.tx", verbose = T))
 # third edit:
for better veiw we can replace the jenus name Champia with C. because all of the species have the same genus name       
names(Alignment_Champia_COI5P_1) <- gsub('Champia', "C.", names(Alignment_Champia_COI5P_1))
names(Alignment_Champia_rbcL_1) <- gsub("Champia", "C.", names(Alignment_Champia_rbcL_1))
# attached is the new dendrogram with the small species name.)      
           
#BrowseSeqs function in DECIPHER package opens an html file in a web browser to show the sequences in a DNAStringSet.#BrowseSeqs converts the DNAStringSet into html format for viewing in a web browser. The sequences are colored and easy to recognize sequences with unusual nucleotides.
      BrowseSeqs(Alignment_Champia_COI5P_1)
        BrowseSeqs(Alignment_Champia_rbcL_1)
        #Check the number of gaps across all sequences in the alignment to see if there's any unusual count. The str_count argument for the function lapply, placed after a comma, within the arguments to be passed to lapply.
      lapply(Alignment_Champia_COI5P_1, str_count, "-")
      lapply(Alignment_Champia_rbcL_1, str_count, "-")
      #as.DNAbin function(ape package) can convert the classe from the package Biostrings in BioConductor for consistency within ape package. Converting this data type is essential for further downstream analysis.
      dnaBin_Alignment_Champia_COI5P_1 <- as.DNAbin(Alignment_Champia_COI5P_1)  
      dnaBin_Alignment_Champia_rbcL_1 <- as.DNAbin(Alignment_Champia_rbcL_1)
      #dist.dna function in ape package, computes a matrix of pairwise distances from DNA sequences using a model of DNA evolution. Tamura and Nei model was used to calculate the distances it  allows for different rates of transitions and transversions, heterogeneous base frequencies, and between site variation of the substitution rate.
      distanceMatrix_COI5P_1 <- dist.dna(dnaBin_Alignment_Champia_COI5P_1, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
      distanceMatrix_rbcL_1<- dist.dna(dnaBin_Alignment_Champia_rbcL_1, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
      
      #Cluster according to unique species names on COI-5P gene. IdClusters() function from the package DECIPHER, cluster according to single linkage method and using a 0.02 (2%) divergence threshold. The method single gave a low scale bar compared to other methods. The type is set up with dendrogram as here I need to visualize the clusters in dendrogram.  
      clusters_COI5P_1 <- IdClusters(distanceMatrix_COI5P_1,
                                     method = "single",
                                     cutoff= 0.02,
                                     showPlot = TRUE,
                                     type = "dendrogram",
                                     verbose = TRUE)
      #the 8th sequence seem to be an outlier but check if the genbank accession number is same as the previously suspected outlier (The sequence which was BLAST and confirmed that it in not an outlier anymore.
      dfChampia_COI5P_Subset[8,"genbank_accession"]
      #the genbank accession number is as same as the previously BLAST sequence. Therefore no need to BLAST again as this sequence was recognized as an assential data in the data set. Exclusion of this data will result data missing and leads to incorrect result.
      #Cluster according to unique species names on rbcL gene.
      clusters_rbcL_1 <- IdClusters(distanceMatrix_rbcL_1,
                                    method = "single",
                                    cutoff= 0.02,
                                    showPlot = TRUE,
                                    type = "dendrogram",
                                    verbose = TRUE)
        
      #Give a color to the labels of species name and set the font size for a nice visualization  
      dend_1<- clusters_COI5P_1 %>% 
        set("labels_col", rainbow(7)) %>%  # change the label color 
        set("labels_cex", 0.7)   #change the font size of the label
      
      dend_2<- clusters_rbcL_1 %>% 
        set("labels_col", rainbow(15)) %>%  # change the label color 
        set("labels_cex", 0.7)   #change the font size of the label
      
      #par(mfrow = c(1, 2)) #par to generate panels with 1 row of 2 graph
      dend_1 %>% plot(main = "Champia phylogeny on COI-5P ") #Plot the dendrogram with the names of the species as the tip label
      dend_2 %>% plot(main = "Champia phylogeny on rbcL")  
      
      #par(mfrow = c(1, 1)) #par to generate panels with 1 row of 1 graph
      circlize_dendrogram(dend_1, labels_track_height = 0.3, dend_track_height = 0.5) #plot circlize dendrogram for nice visualization, adjust the height of labels and the dendrogram
      title("Champia phylogeny on COI-5P gene") #Display the title
      
      circlize_dendrogram(dend_2,  labels_track_height = 0.3, dend_track_height = 0.5) #plot circlize dendrogram for nice visualization
      title("Champia phylogeny on rbcL gene") #Display the title
      
          # To identify whether the two genes (COI-5P and rbcL) yield the same or different phylogenetic hypotheses, need to create a dataset of overlapped species.
      #First check the species names in the "dfChampia_COI5P_Subset" and "dfChampia_rbcL_Subset" datasets
      dfChampia_COI5P_Subset[ , "species_name"]
      dfChampia_rbcL_Subset[ ,"species_name"]
      #Now check the non-overlapping species names in COI-5P gene and rbcL gene
      dplyr::setdiff(dfChampia_COI5P_Subset$species_name, dfChampia_rbcL_Subset$species_name) #Species that appear in dfChampia_COI5P_Subset but not dfChampia_rbcL_Subset are "Champia inkyua", "Champia recta", "Champia sp. 1parvula"
      dplyr::setdiff(dfChampia_rbcL_Subset$species_name, dfChampia_COI5P_Subset$species_name) #Species that appear in dfChampia_rbcL_Subset but not dfChampia_COI5P_Subset are "Champia chathamensis", "Champia compressa", "Champia harveyana", "Champia japonica", "Champia lubrica", "Champia lumbricalis", "Champia puertoricensis" "Champia salicornioides", "Champia sp. CLT-219", "Champia sp. CLT-220", "Champia sp. MS-2012-1", "Champia vieillardii", "Champia viridis" 
      
      #When the overlap dataset gets both the COI-5P and rbcL has nucleotides columns. In downstream analysis this would create errors. Therefore,the columns were rename as COI-5P and rbcL nucleotides respectively 
      names(dfChampia_COI5P_Subset)# check the column number of "nucleotides"
      colnames(dfChampia_COI5P_Subset)[72]<-"COI5P_nucleotides" #change the column name
      colnames(dfChampia_rbcL_Subset)[72]<-"rbcL_nucleotides" #rename the column "nucleotides"
      names(dfChampia_COI5P_Subset) #Check the renamed column in the column number 72
      #create the overlap dataset using the column species name
      dfOverlap <- merge(dfChampia_COI5P_Subset, dfChampia_rbcL_Subset, by = "species_name", all = F)
      dim(dfOverlap) # check the number of rows and columns
      #Align sequences in the overlap dataset separately for COI-5P and rbcL genes
      #change the format of the nucleotide data to be suitable for downstream analysis. Here, the sequence data were changed to a DNAStringSet object class, using a function from the Biostrings package.
      dfOverlap$COI5P_nucleotides <- DNAStringSet(dfOverlap$COI5P_nucleotides)
      dfOverlap$rbcL_nucleotides <- DNAStringSet(dfOverlap$rbcL_nucleotides)
      #names function uses to set the names of nucleotides with the species names, which gives species name as the tip labels  in phylogenetic tree for downstream analysis
      names(dfOverlap$COI5P_nucleotides) <- dfOverlap$species_name
      names(dfOverlap$rbcL_nucleotides) <- dfOverlap$species_name
      
      #MUSCLE is more theoretically reliable.It looks better for translations. It is significantly faster than the other methods  even  for large dataset. Also MUSCLE, has several different alignement algorithms to choose from (for different types of data) and has many useful extra options to work with alignments#To find out what default settings are being used for all the arguments, a log file can be written to disk using the log argument in conjunction with the verbose argument, log = "log.txt", verbose =T. This will write out the default values to the file log.txt in the current working directory of R. (verbose =T generates outputs).
      Alignment_Champia_COI5P_2 <- DNAStringSet(muscle::muscle(dfOverlap$COI5P_nucleotides , log = "log.tx", verbose = T))
      Alignment_Champia_rbcL_2 <- DNAStringSet(muscle::muscle(dfOverlap$rbcL_nucleotides , log = "log.tx", verbose = T))
      #BrowseSeqs function in DECIPHER package opens an html file in a web browser to show the sequences in a DNAStringSet.#BrowseSeqs converts the DNAStringSet into html format for viewing in a web browser. The sequences are colored and easy to recognize sequences with unusual nucleotides.
      BrowseSeqs(Alignment_Champia_COI5P_2)
      BrowseSeqs(Alignment_Champia_rbcL_2)
      #Check the number of gaps across all sequences in the alignment to see if there's any unusual count. The str_count argument for the function lapply, placed after a comma, within the arguments to be passed to lapply.
      lapply(Alignment_Champia_COI5P_2, str_count, "-")
      lapply(Alignment_Champia_rbcL_2, str_count, "-")
      #as.DNAbin function(ape package) can convert the classe from the package Biostrings in BioConductor for consistency within ape package. Converting this data type is essential for further downstream analysis.
      dnaBin_Alignment_Champia_COI5P_2 <- as.DNAbin(Alignment_Champia_COI5P_2)  
      dnaBin_Alignment_Champia_rbcL_2 <- as.DNAbin(Alignment_Champia_rbcL_2)
      #dist.dna function in ape package, computes a matrix of pairwise distances from DNA sequences using a model of DNA evolution. Tamura and Nei model was used to calculate the distances it  allows for different rates of transitions and transversions, heterogeneous base frequencies, and between site variation of the substitution rate.
      distanceMatrix_COI5P_2 <- dist.dna(dnaBin_Alignment_Champia_COI5P_2, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
      distanceMatrix_rbcL_2<- dist.dna(dnaBin_Alignment_Champia_rbcL_2, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
      #Cluster and display dendrograms according to unique species names on COI-5P gene in overlaped dataframe 
      clusters_COI5P_2 <- IdClusters(distanceMatrix_COI5P_2,
                                     method = "single",
                                     cutoff= 0.02,
                                     showPlot = TRUE,
                                     type = "dendrogram",
                                     verbose = TRUE)
      title("Dendrogram of Champia species on COI-5P gene ") #Display the title of Dendrogram
      
      #Cluster according to unique species names on rbcL gene in overlaped dataframe 
      clusters_rbcL_2 <- IdClusters(distanceMatrix_rbcL_2,
                                    method = "single",
                                    cutoff= 0.02,
                                    showPlot = TRUE,
                                    type = "dendrogram",
                                    verbose = TRUE)
      title("Dendrogram of Champia species on rbcL gene ") #Display the title of Dendrogram
      #Display species in overlapped dataframe
      dfOverlap[ , "species_name"]
        
      #Give a color to the labels of species name and set the font size for a nice visualization 
      dend_COI5P <- clusters_COI5P_2 %>% 
        set("labels_col", rainbow(7)) %>%  # change the label color 
        set("labels_cex", 1) #set the font size of the label
      dend_rbcL <- clusters_rbcL_2 %>% 
        set("labels_col", rainbow(7)) %>%  # change the label color 
        set("labels_cex", 1) #set the font size of the label
      
      #Now using above dendrograms perform a tanglegram to compare the sequence clusters
      #A dendlist is a function in dendextend package, which produces the dendlist class. This function uses to compare two dendrograms (clusters_COI5P_2, clusters_rbcL_2) by chaining them together.
      dl <- dendlist(dend_COI5P, dend_rbcL)
      #Here tanglegram creates using dendlist. Tanglegram is a function in dendextend package, which gives two dendrogram (with the same set of labels. To get this, previously, I created an overlap dataframe of COI-5P and rbcL genes),by facing them one in front of the other and connecting their labels with lines. This is used for visually comparing two dendrograms. Nodes which contains a combination of labels  , which are not present in the other tree, can be turned on by using highlight_distinct_edges = TRUE. Also connecting lines are colored to highlight two sub-trees which are present in both dendrograms. This can be turned on by setting common_subtrees_color_lines = TRUE,I set highlight_branches_lwd = FALSE, which keep the lines thin. left and right trees named as COI-5P and rbcL genes respectively and set the margin_inner 7 and title size (cex_main = 1.3) to get an organized tanglegram.
      tanglegram(dl, common_subtrees_color_lines = TRUE, highlight_distinct_edges  = TRUE, highlight_branches_lwd = FALSE, main = "Tanglegram of COI-5P and rbcL marker",main_left = "COI-5P gene", main_right = "rbcL gene", margin_inner=7, cex_main = 1.3)
      
        #Here dendlist concatenates the two dendrograms together. Unique nodes (Distinct branches) are highlighted with  dashed lines, common nodes are in solid lines and according to this tanglegram   the species connecting lines can be seen in between the two tress.  
      #Two trees looked similar according to the tamglegram 
      #entanglement, to measure the quality of the alignment of the two trees in the tanglegram layout. lower value is better
      
      dl %>% entanglement
      #entanglement value is 0.1037517, lower the value and closer to zero means a high quality of the alignment of two trees in tanglegram layout.
      
      #all.equal function makes a global comparison of two dendrograms trees clusters_COI5P_2, clusters_rbcL_2.
      all.equal(clusters_COI5P_2, clusters_rbcL_2) #This will result the difference in branch heights.  Mean relative difference between two dendrogram is  0.4263048
      
      cor.dendlist(dl) #  create a correlation matrix using the cor.dendlist function. corelation between the two dendrogram is 0.8952775  , a strong corelation.
      corrplot(cor.dendlist(dl), "shade", "full", title="corelation plot of  clusters_COI5P_2, clusters_rbcL_2 dendrograms") # Easily shows the corelation between dendrograms by a nice visualization.  The method shade implies the color intensity and the type full shows  relationship between two dendrogram using four squares. The color intensity shows the corelation value (-1 to 1) dark blue implies higher corelation.Therefore the corelation between clusters_COI5P_2, clusters_rbcL_2 is a higher value.
      
          cor_bakers_gamma(clusters_COI5P_2, clusters_rbcL_2) #Baker’s gamma The value can range between -1 to 1. With near 0 values meaning that the two trees are not statistically similar. It shows us that the tree’s topology is not identical as the answere is closer to 0. But here the  Baker’s Gamma Index is a high value of 0.8842063, so there is a high similarity between the two trees' toplogy.
       #Further perform a permutation test for the calculation of the statistical significance of the Baker’s gamma index. Here  the distribution of Baker’s Gamma Index performs under null hypothesis
      set.seed(102)
      the_cor <- cor_bakers_gamma(clusters_COI5P_2, clusters_COI5P_2)#compare clusters_COI5P_2 itself 
      the_cor2 <- cor_bakers_gamma(clusters_COI5P_2, clusters_rbcL_2)#compare clusters_COI5P_2 with clusters_rbcL_2
      
      the_cor #display the correlation coefficient with compared to a tree it self 
      the_cor2 #display the correlation coefficient with compared to clusters_COI5P_2 tree with clusters_rbcL_2  
      R <- 100 #assign 100 to R
      cor_bakers_gamma_results <- numeric(R)
      dend_mixed <- clusters_COI5P_2 #assigning the cluster with a new name dend_mixed
      #plot 100 points iteratively to get a nice grah distributed in the range -1 to 1
      for(i in 1:R) {
        dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE) # Samples a tree without replacing, either by permuting the labels (which is usefull for a permutation test)
        cor_bakers_gamma_results[i] <- cor_bakers_gamma(clusters_COI5P_2, dend_mixed) #Calculate Baker's Gamma correlation coefficient for trees and assumes the labels in the two trees fully match. This would iterate for 100 times.
      }
      #Now plot the distribution of Baker’s Gamma Index under the null hypothesis (Both the trees are statistically similar to each other)
      plot(density(cor_bakers_gamma_results),
           main = "Baker's gamma distribution under H0", #null hypothesis : Both the trees are statistically similar to each other
           xlim = c(-1,1)) # create X axis spread in the range -1 to 1  as the values of bakers gama index range  between -1 and 1.
      
      abline(v = the_cor, lty = 2, col = 2) # draw the line for correlation coefficient with compared to a tree itself 
      abline(v = the_cor2, lty = 2, col = 4) # draw the line for correlation coefficient with compared to clusters_COI5P_2 tree with clusters_rbcL_2 
      legend("topleft", legend = c("cor", "cor2"), fill = c(2,4)) #display the details of the plot
      
      #with this plot, there are enough evidence that clusters_COI5P_2 tree and clusters_rbcL_2  are significantly “similar” with a correlation close to 1.
      
      cor_cophenetic(clusters_COI5P_2, clusters_rbcL_2)# The function cor_cophenetic is faster than cor_bakers_gamma.  The value can range between -1 to 1. With near 0 values meaning that the two trees are not statistically similar & calculating the distribution under the null hypothesis (Both the trees are statistically similar to each other) 
        #According to the results it is obvious that clusters_COI5P_2 and clusters_rbcL_2 are statistically similar with a 0.8952775 cor_cophenetic which is closer to 1. Therefore, the two genes (COI-5P and rbcL) yield statistically similar phylogenetic hypotheses
      
      
      
      
      
      
      
      
      
      
      
      
