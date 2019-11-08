#####
# Part C 

#####
# Question 1

# The biodiversity questions that I am interested in asking are listed below: 
# How well sampled is the taxonomic group as a whole? 
# What is the specis diversity among each region, in each country sampled. 
# How dissimilar is the BIN composition between the two most well-sampled countries? 
# How many regions have been sampled for each country that has been sampled? 

#####
# Question 2

# Use BOLD API tool to download data directly into R. The data were also downloaded and are attached in this assignment. File was downloaded Sunday, September 29, 2019.
library("tidyverse")

country_df <- function(df, country.list,i) {
  output <- data.frame
  output <- df %>% 
    filter(country == country.list[i]) %>% 
    remove_rownames %>% 
    column_to_rownames(var="region") %>% 
    select(-country)
  output[is.na(output)] <- 0
  dfName <- country.names[i]
  assign(dfName, output, env=.GlobalEnv)
}

Porifera <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Porifera&format=tsv")

# Write file to hard disk
write_tsv(Porifera, "Porifera_BOLD_data.tsv")

# Get summary of data 
summary(Porifera)

# Which column names do we have? 
names(Porifera)

# What countries have been sampled? 
unique(Porifera$country)

# What markercodes are represented in our sample? 
unique(Porifera$markercode)

# Which families are represented in the sample? 
unique(Porifera$family_name)

# Which orders are represented in the sample? 
unique(Porifera$order_name)

# What is the highest sampled country? How many BINs do we have for each country?
count.by.country <- Porifera %>%
  group_by(country, bin_uri) %>%
  filter(!is.na(country)) %>% 
  filter(!is.na(bin_uri)) %>% 
  count(bin_uri)

Porifera.spread <- spread(count.by.country, bin_uri, n)
Porifera.spread[is.na(Porifera.spread)] <- 0

country.names <- unique(Porifera.spread$country)

# What is the country with the highest 

family.count.by.country <- Porifera %>% 
  group_by(country, family_name) %>% 
  filter(!is.na(country)) %>% 
  filter(!is.na(family_name)) %>% 
  count(family_name)

family.count.region.country <- Porifera %>% 
  group_by(country, family_name, region) %>% 
  filter(!is.na(country)) %>% 
  filter(!is.na(family_name)) %>%
  filter(!is.na(region)) %>% 
  count(region) %>% 
  arrange(desc(n))

print(family.count.region.country)

species.count.by.country <- Porifera %>% 
  group_by(country,species_name) %>% 
  filter(!is.na(country)) %>% 
  filter(!is.na(species_name)) %>% 
  count(country)

species.count.region.country <- Porifera %>% 
  group_by(country, species_name, region) %>% 
  filter(!is.na(country)) %>% 
  filter(!is.na(species_name)) %>%
  filter(!is.na(region)) %>% 
  count(region) %>% 
  arrange(desc(country))

# Subset each country into its own df and 

for (i in seq_along(Porifera.spread)) {
  country_df(Porifera.spread,country.names,i)
}


# Subset each country into its own df

Porifera.by.country <- Porifera %>%
  group_by(country, bin_uri) %>%
  filter(!is.na(country)) %>%
  filter(!is.na(bin_uri)) %>%
  count(bin_uri)

Porifera.spread.by.country <- Porifera.by.country %>%
  remove_rownames %>%
  column_to_rownames(var="country")

Porifera.spread.by.country[is.na(Porifera.spread.by.country)] <- 0


# Calculate species accumulation for each country

commAssign1.spec <- specaccum(Porifera.spread.by.country)
plot(commAssign1.spec)

par(mar = c(1,1,1,1))
Brazil.spec <- specaccum(Brazil)
plot(Brazil.spec)
Canada.spec <- specaccum(Canada)
plot(Canada.spec, col = "2")
