
if(!file.exists('process_results.R')){
  setwd(getSrcDirectory(function(){})[1])
}

# get costs and deliveries

setwd('costing')
source('costing_script.R')
rmarkdown::render('README.Rmd', 'bookdown::github_document2')
rmarkdown::render('README.Rmd', 'bookdown::pdf_document2', clean=F)
setwd('..')

# run matlab impact model
require(matlabr)
matlabr::run_matlab_script('impact_script.m')

# aggregate impacts to pandemics
source('aggregation_script.R')

# build results markdown
require(rmarkdown)
setwd('..')
rm(params)
rmarkdown::render('outputs.Rmd', 'bookdown::pdf_document2', params=list(lbfile ='data/vaccine_delivery.xlsx'), clean=F)

