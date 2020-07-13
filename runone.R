## This script allows us to run one field 
source('chick_pipp.R')

nreps = 99

type = "Blue"; field = "DN6"



## Nothing to edit below here...

stopifnot(is.element(field, fields))
stopifnot(is.element(type, types))


expt.col = "blue"

deltas = seq(from=5, to=10, by=1)
sigma = seq(from=1, to=10, by=2)
kappa = seq(from=1, to=6, by=2)

grid = as.matrix(expand.grid(deltas, sigma, kappa))

ngrid = nrow(grid)



options(mc.cores=30)
results = mclapply(1:ngrid, run.one)

results_df = do.call(rbind, results)

write.csv(results_df, file.path("res", sprintf("%s_%s.csv", field, type)))



