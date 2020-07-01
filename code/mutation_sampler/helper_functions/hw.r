library(doMC)  # Multicore support 
no.cpus = 4
registerDoMC(no.cpus)  # No. of cores to use

do.in.parallel = TRUE