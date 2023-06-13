
get_most_recent <- function(dataset){
  setDT(dataset)
  dataset[,mostrecent:=max(date),by=country]
  subset(dataset,date==mostrecent)
}


indicators[grepl('television',indicators$indicator),]
# IT.TVS.HOUS.ZS
indicators[grepl('radio',indicators$indicator),]
# IT.RAD.HOUS.ZS
get_most_recent(wb_data('IT.TVS.HOUS.ZS') %>%
                  filter(!is.na(IT.TVS.HOUS.ZS)&date<2020))
get_most_recent(wb_data('IT.RAD.HOUS.ZS') %>%
                  filter(!is.na(IT.RAD.HOUS.ZS)&date<2020))

indicators <- wb_indicators()
View(indicators[grepl('BMI',indicators$indicator),])

meanbmiadultsid <- 'HF.STA.BM18'
meanbmiadults <- wb_data(meanbmiadultsid) %>%
  filter(!is.na(HF.STA.BM18)&date<2020)

get_most_recent(meanbmiadults)

View(indicators[grepl('tourism',indicators$indicator),])

#ST.INT.RCPT.CD International tourism, receipts (current US$)
tourism <- wb_data(c('ST.INT.RCPT.CD','NY.GDP.MKTP.CD')) %>%
  filter(!is.na(ST.INT.RCPT.CD)&date<2020)
tourism <- get_most_recent(tourism)
tourism[,tourpc:=ST.INT.RCPT.CD/NY.GDP.MKTP.CD*100]
subset(tourism,tourpc>30)

tourism <- wb_data(c('ST.INT.RCPT.CD','NY.GDP.MKTP.CD','NV.SRV.TRAD.CD')) %>%
  filter(!is.na(ST.INT.RCPT.CD)&!is.na(NV.SRV.TRAD.CD)&date<2020)
tourism <- get_most_recent(tourism)
tourism[,tourpc:=ST.INT.RCPT.CD/NY.GDP.MKTP.CD*100]
tourism[,serpc:=NV.SRV.TRAD.CD/NY.GDP.MKTP.CD*100]
subset(tourism,tourpc>20)
ggplot(tourism) + 
  annotate('segment',x=0,xend=40,y=0,yend=40) +
  geom_label(aes(x=serpc,y=tourpc,label=country)) +
  theme_bw(base_size=15) +
  labs(x='Services, % of GDP',y='International tourism, % of GDP')

# NV.SRV.TOTL.CD Services correspond to ISIC divisions 50-99. 
# They include value added in wholesale and retail trade (including hotels and restaurants), 
# transport, and government, financial, professional, and personal services such as education, 
# health care, and real estate services. Also included are imputed bank service charges and import duties.

# NV.SRV.OTHR.CD Other services is a subset of services, comprising real estate, renting and business activities 
# (excluding services of owner occupied dwellings), education, health and social work, other community, 
# social and personal service activities, private households with employed persons, and extra territorial 
# organizations (ISIC 70-74, 80-99).

# NV.SRV.DWEL.CD	Dwellings is a subset of services (part of ISIC 70)

# NV.SRV.ADMN.CD  Public administration and defense is a subset of services (ISIC 75).

# NV.SRV.TRAN.CD Transport is a subset of services, comprising transport, storage and communications (ISIC 60-64).

# NV.SRV.TRAD.CD Trade is a subset of services, comprising wholesale and retail trade, and hotel and restaurants (ISIC 50-55)

View(indicators[grepl('services',indicators$indicator),])


# faa from oecd, 2017, for albania, georgia, north macedonia, senegal, serbia:
c(33738.691,1437.469,9709.805,161527,65733.01)/c(1355126.08,35347.67,535726.218,11025109,3952669.83)*100

# suriname
1702452/31482516*100
1170110/17947089*100

#cayman islands
120628/2433964*100
129245/2541559*100
