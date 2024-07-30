library(ggplot2)
library(MASS)
library(squire)

countrydatafile <- '../../data/country_data.csv'

p2data <- read.csv(countrydatafile)
p2countries <- p2data$country

cmcols <- which(grepl('CM',colnames(p2data)))
Npopcols <- which(grepl('Npop',colnames(p2data)))

# prepare squire
populations <- squire::population
countries <- unique(populations$country)

# rename p2 countries to match squire
p2countries[p2countries=='Palestine'] <- 'State of Palestine'
p2countries[p2countries=='Laos'] <- 'Lao PDR'
p2countries[p2countries=='Taiwan'] <- 'China, Taiwan Province of China'
p2countries[p2countries=='Congo'] <- 'Republic of the Congo'
p2countries[p2countries=='Brunei'] <- 'Brunei Darussalam'
p2countries[p2countries=='Kyrgyzstan'] <- 'Kyrgyz Republic'
p2countries[p2countries=='Hong Kong'] <- 'Hong Kong SAR, China'
p2countries[p2countries=='Macao'] <- 'Macao SAR, China'
p2countries[p2countries=='Virgin Islands US'] <- 'United States Virgin Islands'
p2countries <- gsub('St ','St. ',p2countries)

# check what we are still missing
p2countries[!p2countries%in%countries]
countries[!countries%in%p2countries]

# the countries we will replace
rowindices <- p2countries%in%countries
countries_to_get <- p2countries[rowindices]

# get all matrices
matrices <- lapply(countries_to_get,function(cn)
  c(squire:::process_contact_matrix(get_mixing_matrix(cn), get_population(cn,simple_SEIR = T)$n)))
names(matrices) <- countries_to_get


# update populations and matrices
p2data[!rowindices,cmcols] <- NA
for(i in 1:length(matrices)){
  # which country
  cn <- names(matrices)[i]
  rowindex <- which(p2countries==cn)
  # replace matrix
  p2data[rowindex,cmcols] <- matrices[[i]]
  # get saved population
  oldpop <- p2data[rowindex,Npopcols]
  # get new population
  newpop <- subset(populations,country==cn)$n
  # trim to fit age groups
  if(sum(oldpop[18:21],na.rm=T)>newpop[17])
    print(c(sum(oldpop[17:21]),newpop[17]))
  oldpop[1:16] <- newpop[1:16]
  oldpop[17:21] <- as.numeric(oldpop[17:21])*newpop[17]/sum(oldpop[17:21])
  oldpop[is.na(oldpop)] <- 0
  # overwrite
  p2data[rowindex,Npopcols] <- oldpop
}
# View(p2data[,c(1,cmcols)])

write.csv(p2data,countrydatafile,row.names = F, na="")


