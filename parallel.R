nsim = 50
#myseeds = floor(1e4 * runif(nsim))
#save(myseeds, file = "myseed.RData")
load("~/fMRI_data/myseed.RData")

for (i in 1:nsim){
  system(paste0('sbatch -o main.out -t 10:00:00 -n 4 --mem-per-cpu=4g -N 1-1 --wrap="Rscript group_shen.R  myseed=', myseeds[i], '"'))
}
