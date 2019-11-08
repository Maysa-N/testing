####Software Tools (BINF*6210) – Assignment_1_Part_C
#Here I used specimen data from BOLD for the genus Agrilus. As I found the particular species Agrilus planipennis Fairmaire (emerald ash borer) was detected throughout southwestern Ontario, in Ottawa and nearby counties in eastern Ontario, and in Sault Ste.Marie and on Manitoulin Island in northern Ontario. Native ashes have been attacked by Emerald ash borer while the Ashes are an important component in many forest communities in Canada. Therefore, genus Agrilus was used in this project to analyse the data set as the initial step.(https://www.nrcan.gc.ca/forests/fire-insects-disturbances/top-insects/13395)

#How many BINs of the taxonomic group Agrilus have been DNA barcoded from, the first three countries with the highest count of BINs? How dissimilar is the nucleotide composition among these countries and also, between markercodes (COI5P and COI3P) in Agrilus?  How well sampled is the taxonomic group Agrilus? By calculating pairwise dissimmilarity among countries, analyse the dissimilarity to explore variability in community compossition. By using rarefaction, assess expected species richness for Agrilus based on the construction of rarefaction curve. 

#LOAD PACKAGES
#load the library to access to the functions
library(permute)

library(lattice)

library(vegan)

library(tidyverse)

library(ape)

library(RSQLite)

library(BiocManager)

library(DECIPHER)

library(muscle)

library(Biostrings)

#vignette(package="Biostrings")

##Load data.Here work with Agrilus data from BOLD
#Agrilus <- read_tsv("Agrilus_BOLD_data.tsv")
Agrilus<- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Agrilus&format=tsv")

#How many Bins of Agrilus have been DNA barcoded in the countries with the highest counts of Bins,
summary(Agrilus)#obtain a summary that some users may find more readable.  

names(Agrilus)#Getting the names and column numbers of the variables

#Number of records according to the country
Agrilus_country <- table(Agrilus$country)#assign results to Agrilus_country
table(Agrilus$country)#display the result on the screen
  
#do the same thing by piping
Agrilus_by_country <- Agrilus %>% #Name of the dataframe to work with
  group_by(country) %>% #group the data by country
  summarise(count=length(processid)) %>% #summarize the data for each country#length shows the number of records in each country
  arrange(desc(count)) %>% ##Sort the results in descending order
  print() #result print on the screen

#Number of records according to bin_uri
Agrilus_by_bin <- Agrilus %>% 
  group_by(bin_uri) %>% 
  count(bin_uri) %>% 
  print()

#How many taxonomic bins are in the dataset? 
length(unique(Agrilus$bin_uri)) #number of elements in the given vector (bin_uri)

#How many taxonomic species names are in the dataset? 
length(unique(Agrilus$species_name)) #number of species names in the given vector (species_name)

#count the number of specimens per BIN per country
Agrilus_by_bin_country <- Agrilus %>% 
  group_by(bin_uri,country) %>% 
  count(bin_uri) %>% 
  print()

#spread() from the tidyr package (within tidyverse suite) to get the data into the comm data object format. row:here all samples are treated as belongs to one big site, columns are BIN identifiers  and values in cells are the counts of individual per BIN.This function is needed for biodiversity and community ecology analysis using the vegan package.
Agrilus_by_bin_spread <- spread(Agrilus_by_bin,bin_uri,n)
par(mfrow = c(1, 1))
rarecurve(Agrilus_by_bin_spread) #The shape of the curve implies the rate of discovery of species, but it is still increasing not a plateau. prediction:There are many species remain to be sampled. Therefore by further sampling, can add species to the database. 

#now remove the rows where country and bin display NA. These are specimen records on BOLD where sequencing hasn't been performed yet, where sequencing failed, or where the sequence was too short to be assigned to a BIN.
Agrilus_by_country_bin_na_rm <- Agrilus_by_bin_country[c((!is.na(Agrilus_by_bin_country$country)) & (!is.na(Agrilus_by_bin_country$bin_uri))), ] 

#Do the same thing with piping
Agrilus_by_country_bin_na_rm1 <- Agrilus_by_bin_country %>% 
  filter((!is.na(country)) & (!is.na(bin_uri)) )

#as an extra step check whether the both commands give the same result.
all.equal(Agrilus_by_country_bin_na_rm, Agrilus_by_country_bin_na_rm1)

#spread the data accorcing to the bin and country to create a comm object
Agrilus_by_bin_spread_na_rm <- spread(Agrilus_by_country_bin_na_rm,bin_uri,n) #The first argument (Agrilus_by_country_bin_na_rm) to "spread"  Agrilus_by_country_bin_na_rm data,the  second argument (bin_uri) is the column we want to spread by and the third argument (n) is the count of the data.

#converting NAs to zeroes, as required for downstream analysis. Zeroes in Agrilus_by_country_bin_na_rm represent that a given BIN hasn't been barcoded in that country
Agrilus_by_bin_spread_na_rm[is.na(Agrilus_by_bin_spread_na_rm)] <- 0

#investigate the structure of data
str(Agrilus_by_bin_spread_na_rm)

#Typically, we want all of our data to be in the columns and cells, not as a names attribute to proceed with the data format wanted by vegan.

#Set the rownames as country, rather than having country as a data column. 
Agrilus_by_bin_spread_na_rm_row_rm <- Agrilus_by_bin_spread_na_rm %>% 
  remove_rownames %>% 
  column_to_rownames(var = "country")# Now country moves to the left margin as row names

#run species accumulation curve analysis. Resampling is performed to see how BINs accumulate as countries are added.
Agrilus_accumulation <- specaccum(Agrilus_by_bin_spread_na_rm_row_rm)
par(mfrow = c(1, 1))
Agrilus_accumulation_plot <- plot(Agrilus_accumulation) # Make a plot of the model.

#Get the number of institutes that record samples
Agrilus_institute <- Agrilus %>% 
  group_by(institution_storing) %>% 
  summarise(count=length(processid))%>% 
  arrange(desc(count)) %>% 
  print()

sum(is.na(Agrilus$country))# number of countries
sum(is.na(Agrilus$bin_uri))# number of bins   

#A specimen may have more than one row if multiple genes have been sequenced for the same specimen. So, it is important to know what genes are there in the dataset
unique(Agrilus$markercode) ##Seeing what markers are available in this dataset

#group by marker and get the count for each marker
Agrilus_count_by_marker <- Agrilus %>%
  filter(!is.na(nucleotides)) %>%
  group_by(markercode) %>%
  summarize(count= length(processid)) %>%
  arrange(desc(count)) %>%
  print()
#Filtering the dataset first to retain those records having a COI-5P markercode. Use "==" for exactly equal.And also want to filter out records lacking a nucleotide sequence. We are looking to see if there is an A, C, G, or T in the nucleotides column.

Agrilus_COI_5P <- Agrilus %>%
  filter((markercode == "COI-5P") & (str_detect(nucleotides, "[ACGT]")))
#Filtering the dataset first to retain those records having a COI-5P markercode.

Agrilus_COI_3P <- Agrilus %>%
  filter((markercode == "COI-3P") & (str_detect(nucleotides, "[ACGT]")))
unique(Agrilus_COI_5P$markercode) # Checking whether the results have only COI5P marker only 
unique(Agrilus_COI_3P$markercode) # Checking whether the results have only COI5P marker only 
class(Agrilus_COI_5P$nucleotides) #Checking the data type of the nucleotides.
class(Agrilus_COI_3P$nucleotides) #Checking the data type of the nucleotides.

#Changing the class of nucleotide data to a DNAStringSet.
Agrilus_COI_5P$nucleotides <- DNAStringSet(Agrilus_COI_5P$nucleotides)
Agrilus_COI_3P$nucleotides <- DNAStringSet(Agrilus_COI_3P$nucleotides)

class(Agrilus_COI_5P$nucleotides) #Checking the data type of the new nucleotide column
class(Agrilus_COI_3P$nucleotides) #Checking the data type of the new nucleotide column

#Calculates the frequency of each nucleotide separately.
Agrilus_5P_nuc_seq <- as.data.frame(letterFrequency(Agrilus_COI_5P$nucleotides,letters = c( "A","T","C","G")))
# Calculate the AT and GC frequency (as a proportion of the total). The function mutate() add a column onto the end of Agrilus_5P_nuc_seq.
Agrilus_5P_nuc_seq_AT_GC_frq <- Agrilus_5P_nuc_seq %>% 
  mutate(ATpropotion= (A+T)/(A+T+C+G),GCpropotion= (G+C)/(A+T+C+G))

# Calculate the AT and GC frequency (as a proportion of the total). The function mutate() add a column onto the end of Agrilus_3P_nuc_seq.
Agrilus_3P_nuc_seq <- as.data.frame(letterFrequency(Agrilus_COI_3P$nucleotides,letters = c( "A","T","C","G")))
Agrilus_3P_nuc_seq_AT_GC_frq <- Agrilus_3P_nuc_seq %>% 
  mutate(ATpropotion= (A+T)/(A+T+C+G),GCpropotion=(G+C)/(A+T+C+G))


par(mfrow = c(2, 2))#  use par to generate panels with 2 row of 4 graphs
hist(Agrilus_5P_nuc_seq_AT_GC_frq$ATpropotion,col = rainbow(10)) #viewing histogram of Agrilus_5P_nuc_seq_AT_GC_frq$ATpropotion
hist(Agrilus_3P_nuc_seq_AT_GC_frq$ATpropotion,col = rainbow(11)) #viewing histogram of Agrilus_3P_nuc_seq_AT_GC_frq$ATpropotion
hist(Agrilus_5P_nuc_seq_AT_GC_frq$GCpropotion,col = rainbow(12)) #viewing histogram of Agrilus_5P_nuc_seq_AT_GC_frq$GCpropotion
hist(Agrilus_3P_nuc_seq_AT_GC_frq$GCpropotion,col = rainbow(13)) #viewing histogram of Agrilus_3P_nuc_seq_AT_GC_frq$GCpropotion

mean(Agrilus_5P_nuc_seq_AT_GC_frq$ATpropotion) #calculating mean of 5P
mean(Agrilus_3P_nuc_seq_AT_GC_frq$ATpropotion) #calculating mean of 3P 

#Check the mean AT frequencies significantly different between COI5P and COI3P
t.test(Agrilus_5P_nuc_seq_AT_GC_frq$ATpropotion,Agrilus_3P_nuc_seq_AT_GC_frq$ATpropotion)
#Check the mean GC frequencies significantly different between COI5P and COI3P
t.test(Agrilus_5P_nuc_seq_AT_GC_frq$GCpropotion,Agrilus_3P_nuc_seq_AT_GC_frq$GCpropotion)
# reject null hypothesis, There's a significance difference between the nuceotide frequencies (AT and GC) between COI5P and COI3P markercodes. 

#Check whether the results of nucleotide frequencies significantly different among the top most three countries.
a <- Agrilus %>% 
  group_by(country,bin_uri,markercode,nucleotides) %>% 
  summarise(n=length(processid)) %>% 
  arrange(desc(n)) %>% 
  print()
#get a subset of Agrilus with 5 columns country, bin_uri, markercode, nucleotide sequence and count(n).

#remove all the NAs in the columns bin_uri, markercode, country and nucleotides.
a_NA_rm <- a %>% 
  filter(!is.na(bin_uri) & (!is.na(markercode) & (!is.na(country) & (!is.na(nucleotides))))) %>% 
  print()

#number of unique bin_uri according to the country
group_country <- a_NA_rm %>% 
  group_by(country) %>% 
  count(country) %>% 
  print()

# Get the top 3 countries with the highest counts of Bins
Top_three_countries<- sort(table(Agrilus$country), decreasing = TRUE)[1:3] %>% 
print

#Sub set of Canada samples with bin,marker code and nucleotides with the count
Agrilus_Canada <- a_NA_rm %>% 
  filter(country=="Canada" ) %>% 
  print()
sum(Agrilus_Canada$n)
Agrilus_Canada$nucleotides <- DNAStringSet(Agrilus_Canada$nucleotides)
class(Agrilus_Canada$nucleotides)

Agrilus_5P_Canada_nuc_seq <- as.data.frame(letterFrequency(Agrilus_Canada$nucleotides,letters = c( "A","T","C","G")))
Agrilus_5P_Canada_nuc_seq_AT_GC_frq <- Agrilus_5P_Canada_nuc_seq %>% 
  mutate(ATpropotion= (A+T)/(A+T+C+G),GCpropotion= (G+C)/(A+T+C+G))


#Sub set of Slovakia samples with bin,marker code and nucleotides with the count
Agrilus_Slovakia <- a_NA_rm %>% 
  filter(country=="Slovakia" ) %>% 
  print()
sum(Agrilus_Slovakia$n)
Agrilus_Slovakia$nucleotides <- DNAStringSet(Agrilus_Slovakia$nucleotides)
class(Agrilus_Slovakia$nucleotides)

Agrilus_5P_Slovakia_nuc_seq <- as.data.frame(letterFrequency(Agrilus_Slovakia$nucleotides,letters = c( "A","T","C","G")))
Agrilus_5P_Slovakia_nuc_seq_AT_GC_frq <- Agrilus_5P_Slovakia_nuc_seq %>% 
  mutate(ATpropotion= (A+T)/(A+T+C+G),GCpropotion= (G+C)/(A+T+C+G))


#Sub set of Russia samples with bin,marker code and nucleotides with the count
Agrilus_Russia <- a_NA_rm %>% 
  filter(country=="Russia" ) %>% 
  print()
sum(Agrilus_Russia$n)

Agrilus_Russia$nucleotides <- DNAStringSet(Agrilus_Russia$nucleotides)
class(Agrilus_Russia$nucleotides)

Agrilus_5P_Russia_nuc_seq <- as.data.frame(letterFrequency(Agrilus_Russia$nucleotides,letters = c( "A","T","C","G")))
Agrilus_5P_Russia_nuc_seq_AT_GC_frq <- Agrilus_5P_Russia_nuc_seq %>% 
  mutate(ATpropotion= (A+T)/(A+T+C+G),GCpropotion= (G+C)/(A+T+C+G))


par(mfrow = c(2, 3)) # use par to generate panels with 2 rows of 6 graphs,in each row 3 graphs. col=blues, display histogram with dens colors 
hist(Agrilus_5P_Canada_nuc_seq_AT_GC_frq$ATpropotion,col = blues9) 
hist(Agrilus_5P_Slovakia_nuc_seq_AT_GC_frq$ATpropotion,col = blues9)
hist(Agrilus_5P_Russia_nuc_seq_AT_GC_frq$ATpropotion,col = blues9)
hist(Agrilus_5P_Canada_nuc_seq_AT_GC_frq$GCpropotion,col = blues9) 
hist(Agrilus_5P_Slovakia_nuc_seq_AT_GC_frq$GCpropotion,col = blues9)
hist(Agrilus_5P_Russia_nuc_seq_AT_GC_frq$GCpropotion,col = blues9)

#Compare the AT frequences of COI5P in Canada and Slovakia using a t.test. to check whether there is a significant different 
t.test(Agrilus_5P_Canada_nuc_seq_AT_GC_frq$ATpropotion,Agrilus_5P_Slovakia_nuc_seq_AT_GC_frq$ATpropotion)#fail to reject null hypothesis, There's no significance difference of nucleotide frequencies (AT) of COI5P in Canada and Slovakia.

#Compare the GC frequences of COI5P in Canada and Slovakia using a t.test. to check whether there is a significant different
t.test(Agrilus_5P_Canada_nuc_seq_AT_GC_frq$GCpropotion,Agrilus_5P_Slovakia_nuc_seq_AT_GC_frq$GCpropotion)#fail to reject null hypothesis, There's no significance difference of nucleotide frequencies (CG) of COI5P in Canada and Slovakia.

#Compare the AT frequences of COI5P in Canada and Russia using a t.test. to check whether there is a significant different
t.test(Agrilus_5P_Canada_nuc_seq_AT_GC_frq$ATpropotion,Agrilus_5P_Russia_nuc_seq_AT_GC_frq$ATpropotion)# reject null hypothesis, There's a significance difference of nucleotide frequencies (AT) of COI5P in Canada and Russia.

#Compare the GC frequences of COI5P in Canada and Russia using a t.test. to check whether there is a significant different
t.test(Agrilus_5P_Canada_nuc_seq_AT_GC_frq$GCpropotion,Agrilus_5P_Russia_nuc_seq_AT_GC_frq$GCpropotion)# reject null hypothesis, There's a significance difference of nucleotide frequencies (GC) of COI5P in Canada and Russia.

#Compare the AT frequences of COI5P in Slovakia and Russia using a t.test. to check whether there is a significant different
t.test(Agrilus_5P_Slovakia_nuc_seq_AT_GC_frq$ATpropotion,Agrilus_5P_Russia_nuc_seq_AT_GC_frq$ATpropotion)# reject null hypothesis, There's a significance difference of nucleotide frequencies (AT) of COI5P in Slovakia and Russia.

#Compare the GC frequences of COI5P in Slovakia and Russia using a t.test. to check whether there is a significant different
t.test(Agrilus_5P_Slovakia_nuc_seq_AT_GC_frq$GCpropotion,Agrilus_5P_Russia_nuc_seq_AT_GC_frq$GCpropotion)# reject null hypothesis, There's a significance difference of nucleotide frequencies (GC) of COI5P in Slovakia and Russia.

# Check the diversity of the group Agrilus,
#calculate Simpson's 1-D Index of Diversity for each country. (closer to 1 = greater diversity)
Agrilus_simpson_dist <- diversity(Agrilus_by_bin_spread_na_rm_row_rm,"simpson") %>% 
  print()

#Shannon's is default(ypically ranges from 1.5 - 3.4, higher = more diverse )
Agrilus_shannon_dist <- diversity(Agrilus_by_bin_spread_na_rm_row_rm) %>% 
  print()

#par to generate panels with 1 row of 2 graphs
par(mfrow = c(1, 2))
hist(Agrilus_simpson_dist,col=rainbow(10))
hist(Agrilus_shannon_dist,col = rainbow(10))

#According to the simpson's diversity index Canada,China, Costa Rica, Finland, France, Germany, Greece, Italy, Russia, Slovakia, Spain,United States, Vietnam have numerical values closer to 1, which means a high diversity.#According to the Shannon diversity index Canada,China, Costa Rica, Finland, France, Germany, Greece, Italy, Russia, Slovakia, Spain,United States, Vietnam have numerical values lie between the range 1.5- 3.4 , which means a high diversity.Therefore both the methods implies same result.

#calcuate all of the pair-wise dissimilarity (distance) measures between countries based on their bin_uri composition using the function vegdist. Vegdist computes dissimilarity indices. Here uses bray-curtis which is good in detecting underlying ecological gradients.Bray dissimilaarity is used to quantify the compositional dissimilarity between two different countries. They are bounded between 0 and 1, where 0 = same composition, 1 = maximally dissimilar. Dissimilarity analysis is a good way to explore variability in community composition.
Agrilus_bray_dis <- vegdist(Agrilus_by_bin_spread_na_rm_row_rm, "bray") %>% 
  print()

#Agrilus_gower_dis <- vegdist(Agrilus_by_bin_spread_na_rm_row_rm, "gower") %>% 
#print()

#par to generate panels with 1 row of 1 graph
par(mfrow = c(1, 1))
hist(Agrilus_bray_dis, xlim = range(0.0,1.0),col ="pink")
#hist(Agrilus_gower_dis, xlim = range(0.0,1.0),col = "purple")
#According to the histogram, almost all the countries imply maximum dissimilarity as all the frequencies plot around bray distance 1.

#By calculating pairwise dissimmilarity among countries, could identify a maximum dissimilarity which leads to the variability in community compossition.


#Rarefaction is a technique to assess expected species richness. Rarefaction allows the calculation of species richness for a given number of individual samples, based on the construction of rarefaction curves.
rowSums(Agrilus_by_bin_spread_na_rm_row_rm) #gives the number of samples found in each country

#par to generate panels with 1 row of 1 graph
par(mfrow = c(1, 1))
rarecurve(Agrilus_by_bin_spread_na_rm_row_rm, col = "Red") # produces rarefaction curves # squares are the countries 

#The shape of the rarecurve implies the rate of discovery of species,  Rarefaction curves generally grow rapidly at first, as the most common species are found.If the curves plateau then only the rarest species remain to be sampled. But here it is still increasing, not a plateau. The prediction is there are many species remain to be sampled. Therefore species can added to the database by further sampling.
#The  Agrilus diversity is high in Canada,China, Costa Rica, Finland, France, Germany, Greece, Italy, Russia, Slovakia, Spain,United States, and Vietnam. The community composition is highly variable when comparing pairwise dissimilarity. Further, the nucleotide frequencies (AT and GC)are significantly different between the markercodes, COI5P and COI3P. When comparing the AT and GC frequencies of COI5P in Canada, slovaki and Russia, there is no significant different between Canada and Slovaki. while there are significant differences between Canada-Russia, and Russia- Slovaki.
  