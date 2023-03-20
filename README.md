# Daedalus-P2-Dashboard

This is the README file for the Daedalus P2 Dashboard.

Country data is provided in .mat files according to the country name

Disease parameters are provided in p2Params_XXX.m

The main function is p2Sim.m, which takes as inputs the country name, disease name, mitigation strategy, level of preparedness and two additional parameters to scale R0 and severity (set these to 1)

The other functions called by the main code include:
#p2Params.m to define other country, disease and preparedness specific parameters
#p2MakeDs.m to create the contact matrices
#p2Run.m is the main solver
#p2Cost.m calculates the cost of a given epidemic scenario
#p2Plot.m plots the epidemic trajectory and costs
