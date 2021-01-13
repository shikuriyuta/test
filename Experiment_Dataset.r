dir_ = "C:/Users/shiku/Work/Paper/ICML2021/" # directory
data = "asia" # dataset 
N = 100 # sample size
iter_n = 10 # number of samplings

library(bnlearn)
bn.fit = load(paste(dir_, data, "/", data, ".rda", sep=""))
for (num in 1:iter_n) {
    samples = rbn(bn, n=N)
    df <- data.frame(samples)
    write.csv(df, file=paste(dir_, data, "/", data, "_", as.character(N), "_", as.character(num), ".csv", sep=""), row.names=FALSE)
}