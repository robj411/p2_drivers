
fn <- function(){}
folder <- utils::getSrcDirectory(fn)
setwd(ifelse(folder=='','.',folder))


print('1 pathogenprofile.R')
source('pathogenprofile.R')

print('2 rewrite_cm.R')
source('rewrite_cm.R')

print('3 process_cmix.R')
source('process_cmix.R')

print('4 get_distributions.R')
source('get_distributions.R')

print('5 behaviour.R')
source('behaviour.R')

