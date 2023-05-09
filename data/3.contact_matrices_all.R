options(java.parameters = "-Xmx8000m")
library("readxl")
library(xlsx)

countries <- read_xlsx("3.contact_matrices_list.xlsx")

load("contact_work.rdata")#
codes <- names(contact_work)#
print(all(codes==countries$Code))
df <- data.frame(matrix(unlist(contact_work), nrow=16, byrow=FALSE))#

for (i in 1:177){
  a <- 16*(i-1)+1
  b <- 16*i
  write.xlsx(df[a:b], file="4.contact_matrices_work.xlsx", sheetName=countries$Country[i], append=TRUE, row.names=FALSE, col.names=FALSE)#
  gc()
}

extra = list('Australia','Japan','Taiwan','Antigua and Barbuda','Seychelles','Monaco','Andorra','Lebanon','Kiribati','Haiti')

for (i in 1:length(extra)){
  add <- read.xlsx("_4.contact_matrices_work.xlsx", extra[[i]], header=FALSE)#
  write.xlsx(add, file="4.contact_matrices_work.xlsx", sheetName=extra[[i]], append=TRUE, row.names=FALSE, col.names=FALSE)#
  gc()
}