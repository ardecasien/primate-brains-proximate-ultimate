library(mvMORPH)
library(phylopath)
library(phytools)
library(caper)
library(ggplot2)
library(ppcor)
library(igraph)

`%!in%` <- Negate(`%in%`)

## load dN/dS data

dnds = read.csv('dnds_data.csv')
columns_with_inf = sapply(dnds, function(col) any(is.infinite(col)))
table(columns_with_inf)
dnds = dnds[, !columns_with_inf]

## load phenotypic data

pheno = read.csv('phenotype-data.csv')
dim(pheno)

## merge data

data = merge(pheno, dnds, by = 'Species_names')
data = data[complete.cases(data[,c('treename_10k','logbrain','logbody','DQI')]),]
dim(data) # 50 species

## load tree

tree = read.nexus('consensusTree_10kTrees_Primates_Version3.nex')
tree2 = drop.tip(tree,setdiff(tree$tip.label,data$treename_10k))

######################
# run PGLS models
# abs/rel brain size as a function of per gene dNdS
######################

pgls.results = data.frame()

dwarves = 'yes'
dwarves = 'no'

brain = 'absolute'
brain = 'relative'

for(i in 1:8327){
  
  print(i)
  now = data[,c('logbrain','logbody','DQI')]
  now$gene = data[,i+6]
  now$Taxon = data$treename_10k
  rownames(now) = now$Taxon
  
  now = subset(now, gene > 0 & gene < 1)
  now = now[complete.cases(now),]
  
  # use if removing species that are phyletic dwarves
  
  if(dwarves == 'yes') {now = now} else {
  now = subset(now, Taxon %!in% c('Callithrix_jacchus','Microcebus_murinus','Saguinus_midas'))}
  
  if(dim(now)[1] < 10) {print('not enough species')} else {
  
  treenow = drop.tip(tree2,setdiff(tree2$tip.label,rownames(now)))
  brainsphylo = comparative.data(phy = treenow, data = now, names.col = Taxon, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
  
  if(brain == 'absolute') {
    mod1 = pgls(logbrain ~ log(gene), data = brainsphylo, lambda = 'ML', bounds=list(lambda=c(0.02 , 0.98)))} else {
      mod1 = pgls(logbrain ~ logbody + log(gene), data = brainsphylo, lambda = 'ML', bounds=list(lambda=c(0.02 , 0.98)))  
    }
  s = summary(mod1)
  o = data.frame(s$coefficients)
  o$variable = rownames(o)
  o$lamba = s$param[2]
  rownames(o) = NULL
  o$gene = colnames(data)[i+6]
  o$n = dim(now)[1]
  o$model = brain
  o$dwarves = dwarves
  pgls.results = rbind(pgls.results,o)
  
}}

colnames(pgls.results) = c('estimate','se','t','p','predictor','lambda','gene','n','model','dwarves_incl')
pgls.results$padj = p.adjust(pgls.results$p, method =  "BH")

write.csv(pgls.results, file = paste(brain,dwarves,'pgls_results.csv',sep="_"))

######################
# partial correlations
# abs/rel brain size + dNdS + DQI
######################

model.results = data.frame()
partCor.results = data.frame()

dwarves = 'yes'
dwarves = 'no'

for(i in 1:8327){
  
  print(colnames(data)[i+6])
  now = data[,c('logbrain','logbody','DQI')]
  now$gene = data[,i+6]
  now$Taxon = data$treename_10k
  rownames(now) = now$Taxon
  
  now = subset(now, gene > 0 & gene < 1)
  now = now[complete.cases(now),]
  now$logGene = log(now$gene)
  now$logDQI = log(now$DQI)
  
  if(dwarves == 'yes') {now = now} else {
    now = subset(now, Taxon %!in% c('Callithrix_jacchus','Microcebus_murinus','Saguinus_midas'))}
  
  if(dim(now)[1] < 10) {print('not enough species')} else {
    
  treenow = drop.tip(tree2,setdiff(tree2$tip.label,rownames(now)))
  
  Y = as.matrix(now[, c("logGene", "logbrain" ,"logbody", "logDQI")])
  fit1 <- tryCatch(mvgls(Y~1, tree = treenow, model = "OU",  penalty = "RidgeArch", upper = Inf), error=function(err) NA)
  fit2 <- tryCatch(mvgls(Y~1, tree = treenow, model = "BM", penalty = "RidgeArch", upper = Inf), error=function(err) NA)
  fit3 <- tryCatch(mvgls(Y~1, tree = treenow, model = "EB",  penalty = "RidgeArch", upper = Inf), error=function(err) NA)

  Gfit1 <- tryCatch(GIC(fit1)$GIC, error=function(err) NA)
  Gfit2 <- tryCatch(GIC(fit2)$GIC, error=function(err) NA)
  Gfit3 <- tryCatch(GIC(fit3)$GIC, error=function(err) NA)

  model_data = data.frame(Model = c("Model-OU", "Model-BM", "Model-EB"), 
                             GIC = c(Gfit1, Gfit2, Gfit3))
  min_index = which.min(model_data$GIC)
  best_model = model_data$Model[min_index]
  
  new_data = data.frame(colnames(data)[i+6], model_data,best_model)
  model.results = rbind(model.results,new_data,best_model)
  
  if (best_model == "Model-OU") {
    best_fit = fit1  } else if (best_model == "Model-BM") {
    best_fit = fit2 } else if (best_model == "Model-EB") {
    best_fit = fit3 }
  
  # partial correlation 
  N = cov2pcor(best_fit$sigma$Pinv)
  rownames(N) = colnames(N) = c("logGene", "logbrain" ,"logbody", "logDQI")
  N
  
  # significance edge exclusion test
  edge_results <- edge_exclusion_test(C = best_fit$sigma$Pinv, n = dim(now)[1])
  edge_results
  
  p_values <- edge_results$p_values
  significance_levels <- edge_results$matrix
 
  # format
  N_flat <- N[upper.tri(N)]
  edge_labels <- significance_levels[upper.tri(N)]
  edge_labels <- ifelse(edge_labels == "*", paste0(round(N_flat, 3), " *"),"")
  edge_labels
  
  # save results
  
  new_data = data.frame(r = N_flat, p = p_values)
  new_data$variables = c('gene-brain','gene-body','brain-body','gene-dqi','brain-dqi','body-dqi')
  new_data$gene = colnames(data)[i+6]
  partCor.results = rbind(partCor.results, new_data)
  
  # plots
  
  pdf(paste(colnames(data)[i+6],"_diagram.pdf",sep = ""), width = 6, height = 6)
  my_plot_pcor(M = best_fit$sigma$Pinv, 
                row_names = rownames(N), 
               edge_labels = edge_labels,
               0, 0, "Partial correlations diagram")
  dev.off()

  pdf(paste(colnames(data)[i+6],"_corrplot.pdf",sep = ""), width = 6, height = 6)
  corrplot(N, method = "circle", tl.col = "black", tl.cex = 0.8)
  dev.off()

  }
}

saveRDS(model.results, file = paste(dwarves,'partcor_GIC_results.rds',sep="_"))
write.csv(model.results, file = paste(dwarves,'partcor_GIC_results.csv',sep="_"))

partCor.results$padj = p.adjust(partCor.results$p, method = 'BH')
signcomp = data.frame(gene = unique(partCor.results$gene),
                      r_sign_gene_brain = sign(partCor.results$r[partCor.results$variables == "gene-brain"]),
                      r_sign_gene_dqi = sign(partCor.results$r[partCor.results$variables == "gene-dqi"]))
partCor.results = merge(partCor.results, signcomp, by = 'gene', all = T)
sigcomp = subset(partCor.results, variables %in% c('gene-brain','gene-dqi')) 
sigcomp$sig = ifelse(sigcomp$p < 0.1, 'sig', 'ns')
sigcomp = sigcomp %>% group_by(gene,sig) %>% summarise(n=n())
sigcomp = subset(sigcomp, sig == 'sig' & n == 2)
partCor.results = merge(partCor.results, sigcomp[,c('gene','sig')], by = 'gene', all = T)
partCor.results$code <- ifelse(partCor.results$sig == 'sig' &
                                 partCor.results$r_sign_gene_brain == partCor.results$r_sign_gene_dqi,
                               "yes", "no")

saveRDS(partCor.results,  file = paste(dwarves,'partcor_results.rds',sep="_"))
write.csv(partCor.results,  file = paste(dwarves,'partcor_results.csv',sep="_"))

#################
## path analyses
#################

path.results = data.frame()

dwarves = 'yes'
dwarves = 'no'

for(i in 1:8327){
  
  print(i)
  now = data[,c('logbrain','logbody','DQI')]
  now$gene = data[,i+6]
  now$Taxon = data$treename_10k
  rownames(now) = now$Taxon
  
  now = subset(now, gene > 0 & gene < 1)
  now = now[complete.cases(now),]
  now$logGene = log(now$gene)
  now$logDQI = log(now$DQI)
  
  if(dwarves == 'yes') {now = now} else {
    now = subset(now, Taxon %!in% c('Callithrix_jacchus','Microcebus_murinus','Saguinus_midas'))}
  
  if(dim(now)[1] < 10) {print('not enough species')} else {
    
  treenow = drop.tip(tree2,setdiff(tree2$tip.label,rownames(now)))
    
  # define models to be compared
  models = define_model_set(
    A = c(logbrain ~ logbody +logDQI, logbody ~logDQI),
    B = c(logbrain ~ logbody +logDQI, logbody ~logDQI, logbrain ~ logGene),
    C = c(logbrain ~ logbody +logDQI, logbody ~logDQI, logbrain ~ logGene, logGene ~logDQI),
    D = c(logbrain ~ logbody +logDQI, logbody ~logDQI, logbrain ~ logGene, logGene ~ logbody),
    E = c(logbrain ~ logbody +logDQI, logbody ~logDQI, logGene ~logDQI),
    F = c(logbrain ~ logbody +logDQI, logbody ~logDQI, logGene ~ logbody))
  
  # plot path models
  plot_model_set(models)
  
  # run models
  result = phylo_path(models, data = now, tree = treenow, model = 'lambda') # another round with 3“OUfixedRoot”
  s = summary(result)
  rownames(s) = NULL
  s$gene = colnames(data)[i+6]
  s$best_model = s[which(s$CICc == min(s$CICc)),]$model
  s$method = 'lambda'
  path.results = rbind(path.results, s)
  
  result = phylo_path(models, data = now, tree = treenow, model = 'OUfixedRoot') 
  s = summary(result)
  rownames(s) = NULL
  s$gene = colnames(data)[i+6]
  s$best_model = s[which(s$CICc == min(s$CICc)),]$model
  s$method = 'OUfixedRoot'
  path.results = rbind(path.results, s)
  
  result = phylo_path(models, data = now, tree = treenow, model = 'EB') 
  s = summary(result)
  rownames(s) = NULL
  s$gene = colnames(data)[i+6]
  s$best_model = s[which(s$CICc == min(s$CICc)),]$model
  s$method = 'EB'
  path.results = rbind(path.results, s)
  }
}
write.csv(pgls.results, file = paste(brain,dwarves,'path_results.csv',sep="_"))

View(path.results)
