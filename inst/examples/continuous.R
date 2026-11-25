library(MixedDeepMM)
data(iris)
x <- iris[1:4]  # species labels are withheld from fitting
fit <- mixeddeepmm(x, family="gaussian", heads="single",
  layers=list(k=c(2L,2L),r=c(2L,1L)),
  control=list(max_iter=40L,nstart=2L,seed=20L))
print(summary(fit))
cat("ARI:",mixeddeepmm_ari(iris$Species,predict(fit,type="cluster")),"\n")
print(mixeddeepmm_diagnostics(fit))
