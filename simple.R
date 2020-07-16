source('chick_pipp.R')
nreps = 19
expt.col = 'blue'
try_field("DN1", "Blue", 5, 9, 5)  ## should work

try_field("DN1", "Blue", 10, 9, 5) ## should not converge
