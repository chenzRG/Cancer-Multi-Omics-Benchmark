########################################
############## SIMULATION ##############
########################################

get.results.dir = function() {
  return('RESULTS_DIR_PATH')
}

get.simulation.missing.data.range <- function() {
  return(seq(0, 0.8, 0.1))
}

get.num.neighbors.sequence <- function() {
  return(seq(25, 105, by=10))
}

get.missing.values.range <- function() {
  return(c(0, seq(0.1, 0.7, by=0.2)))
}

load.nemo.results.libraries <- function() {
  library('SNFtool')
  library('bnstruct')
  library('fossil')
  library('parallel')
  library('mixtools')
}

get.full.data.results.dir <- function() {
  return('path_to_results_on_full_data')
}

create.all.robustness.results <- function() {
  robustness.results.dir = file.path(get.results.dir(), 'robustness')
  neighbors.range = get.num.neighbors.sequence()
  num.clusters.results.dir = file.path(get.results.dir(), 'robustness_num_clusters')
  num.clusters.range = 2:15
  
  benchmark.ret = analyze.benchmark()
  nemo.surv = benchmark.omics.surv(benchmark.ret)['nemo',]
  nemo.clin = benchmark.omics.clinical(benchmark.ret)['nemo',]
  nemo.num.clusters = benchmark.omics.num.clusters(benchmark.ret)['nemo',]
  full.data.nemo.results = list(survival=nemo.surv, clinical=nemo.clin, num_clusters=nemo.num.clusters)
  
  create.robustness.results(robustness.results.dir, T, neighbors.range, full.data.nemo.results, 'number of neighbors', NULL)
  create.robustness.results(num.clusters.results.dir, F, num.clusters.range, full.data.nemo.results, 'number of clusters', nemo.num.clusters)
}

create.robustness.results <- function(results.dir, is.neighbors, arg.sequence, full.data.results, x.title, nemo.params) {
  set.seed(42)
  
  robustness.results = list()
  add.num.samples = is.null(nemo.params)

  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype.name = current.subtype.data$name
    subtype.display.name = current.subtype.data$display.name
    
    pvalues = c()
    num.enriched = c()
    all.num.clusters = c()
    print(subtype.name)
    
    subtype.raw.data = get.raw.data(subtype.name, 
                                    only.primary = current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, 
                                           current.subtype.data)
					   
    if (add.num.samples) {
      nemo.params = c(nemo.params, round(ncol(subtype.raw.data[[1]]) / NUM.NEIGHBORS.RATIO))
    }
    
    for (argum in arg.sequence) {
      cur.results.file.name = file.path(results.dir, paste(subtype.name, argum, sep='_'))
      if (file.exists(cur.results.file.name)) {
        next
      }
      
      if (is.neighbors) {
        clustering.ret = run.nemo(subtype.raw.data, subtype.name, num.neighbors = argum)
      } else {
        clustering.ret = run.nemo(subtype.raw.data, subtype.name, num.clusters = argum)
      }
      
      clustering = clustering.ret$clustering
      save(clustering, file=cur.results.file.name)
    }
  }
  
  # calculate survival and clinical enrichment
  surv.robustness.results = matrix(0, nrow=length(SUBTYPES.DATA), ncol=length(arg.sequence))
  clin.robustness.results = matrix(0, nrow=length(SUBTYPES.DATA), ncol=length(arg.sequence))
  num.clusters.robustness.results = matrix(0, nrow=length(SUBTYPES.DATA), ncol=length(arg.sequence))
  
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype.name = current.subtype.data$name
    subtype.display.name = current.subtype.data$display.name
    
    
    pvalues = c()
    num.enriched = c()
    all.num.clusters = c()
    print(subtype.name)
    
    subtype.raw.data = get.raw.data(subtype.name, 
                                    only.primary = current.subtype.data$only.primary)
    
    
    for (j in 1:length(arg.sequence)) {
      argum = arg.sequence[j]
      cur.results.file.name = file.path(results.dir, paste(subtype.name, argum, sep='_'))
      load(cur.results.file.name)

      clin.file.path = paste(cur.results.file.name, 'clin', sep='_')
      surv.file.path = paste(cur.results.file.name, 'surv', sep='_')
      if (!file.exists(clin.file.path)) {
	  empirical.clin = check.clinical.enrichment(clustering, subtype.name)
	  save(empirical.clin, file=clin.file.path)  
      } else {
	  # loads empirical clin
	  load(clin.file.path)
      }

      if (!file.exists(surv.file.path)) {
	  empirical.surv = get.empirical.surv(clustering, subtype.name)
	  save(empirical.surv, file=surv.file.path)  
      } else {
	  # loads empirical surv
	  load(surv.file.path)
      }
      surv.robustness.results[i, j] = -log10(empirical.surv$pvalue)
      clin.robustness.results[i, j] = sum(empirical.clin * length(empirical.clin) < 0.05)
      num.clusters.robustness.results[i, j] = max(clustering)
    }
  }
  robustness.results = list(surv.robustness.results=surv.robustness.results, clin.robustness.results=clin.robustness.results, 
                            num.clusters.robustness.results=num.clusters.robustness.results)
  print(robustness.results)
  plot.robustness.results(surv.robustness.results, -log10(0.05), 0.5, '-log10(logrank p-value)', file.path(results.dir, 'survival_robustness.png'), arg.sequence, x.title, list(params=nemo.params, values=full.data.results$survival))
  plot.robustness.results(clin.robustness.results, 1, 1, '# enriched clinical parameters', file.path(results.dir, 'clinical_robustness.png'), arg.sequence, x.title, list(params=nemo.params, values=full.data.results$clinical))
  plot.robustness.results(num.clusters.robustness.results, 1, 1, 'number of clusters', file.path(results.dir, 'num_clusters_robustness.png'), arg.sequence, x.title, list(params=nemo.params, values=full.data.results$num_clusters))
}

analyze.all.missing.real.data <- function() {
  all.subtypes.results = list()
  all.subtypes.with.third.results = list()
  for (i in 1:length(SUBTYPES.DATA)) {
    # we set the seed for each subtype so that we can run them separately
    set.seed(42)
  
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    print(subtype)
    subtype.raw.data = get.raw.data(subtype, 
                                    only.primary=current.subtype.data$only.primary)
   
    num.iters = 5
    missing.probs = get.missing.values.range()
    subtype.with.third = paste(subtype, 'with_third', sep='_')
    
    all.subtypes.results[[i]] = analyze.missing.real.data(subtype, subtype, missing.probs, num.iters, subtype.raw.data[1:2])
    all.subtypes.with.third.results[[i]] = analyze.missing.real.data(subtype.with.third, subtype, missing.probs, num.iters, subtype.raw.data)
  }
  plot.all.missing.real.data(all.subtypes.results, 'two_omics')
  plot.all.missing.real.data(all.subtypes.with.third.results, 'three_omics')
  
  return(list(all.subtypes.results, all.subtypes.with.third.results))
}

plot.all.missing.real.data <- function(all.subtypes.results, plot.prefix) {
  subtype.surv.data = list()
  subtype.clin.data = list()
  
  # survival plot
  png(file.path(get.results.dir(), paste0(plot.prefix, '_all_missing_data_surv.png')), height=600, width=1350)
  par(mfrow=c(2, 5))
  for (i in 1:length(SUBTYPES.DATA)) {
    subtype.surv.data[[i]] = lapply(all.subtypes.results[[i]], function(x) -log10(x$surv.mat))
    plot.partial.results(subtype.surv.data[[i]], get.missing.values.range(), '-log10 logrank pvalue', main=SUBTYPES.DATA[[i]]$display.name, include.labels=(i==1))
  }
  dev.off()
  par(mfrow=c(1, 1))

  alg.scores = reformat.subtype.results(all.subtypes.results, subtype.surv.data)
  png(file.path(get.results.dir(), paste0(plot.prefix, '_mean_missing_data_surv.png')))
  plot.partial.results(alg.scores, get.missing.values.range(), '-log10 logrank pvalue')
  dev.off()
  
  # clinical parameters plot
  png(file.path(get.results.dir(), paste0(plot.prefix, '_all_missing_data_clin.png')), height=600, width=1350)
  par(mfrow=c(2, 5))
  for (i in 1:length(SUBTYPES.DATA)) {
    subtype.clin.data[[i]] = lapply(all.subtypes.results[[i]], function(x) x$clin.mat)
    plot.partial.results(subtype.clin.data[[i]], get.missing.values.range(), '# enriched clinical parameters', main=SUBTYPES.DATA[[i]]$display.name, include.labels=(i==1))
  }
  dev.off()
  par(mfrow=c(1, 1))
  
  alg.scores = reformat.subtype.results(all.subtypes.results, subtype.clin.data)
  png(file.path(get.results.dir(), paste0(plot.prefix, '_mean_missing_data_clin.png')))
  plot.partial.results(alg.scores, get.missing.values.range(), '# enriched clinical parameters')
  dev.off()
}

reformat.subtype.results <- function(all.subtypes.results, criterion.results) {
  alg.scores = list()
  for (i in 1:length(all.subtypes.results[[1]])) {
    for (j in 1:length(SUBTYPES.DATA)) {
      if (j == 1) {
        method.sum = colMeans(criterion.results[[j]][[i]])
      } else {
        method.sum = method.sum + colMeans(criterion.results[[j]][[i]])
      }
    }
    alg.scores[[i]] = method.sum / length(SUBTYPES.DATA)
  }
  names(alg.scores) = names(all.subtypes.results[[1]])
  return(alg.scores)
}



analyze.missing.real.data <- function(subtype.dir, subtype.name, missing.probs, num.iters, subtype.raw.data) {
  alg.names = c('pvc', 'nemo', 'mkl', 'cca', 'cca2', 'nemo_after')
  subtype.results = list()
  for (alg.index in 1:length(alg.names)) {
    alg.name = alg.names[alg.index]
    
    surv.mat = matrix(NA, ncol=length(missing.probs), nrow=num.iters)
    clin.mat = matrix(NA, ncol=length(missing.probs), nrow=num.iters)
    num.cluster.mat = matrix(NA, ncol=length(missing.probs), nrow=num.iters)
    timing.mat = matrix(NA, ncol=length(missing.probs), nrow=num.iters)
    
    res.dir = file.path(get.results.dir(), paste0(alg.name, '_results'), subtype.dir)
    
    for (prob.index in 1:length(missing.probs)) {
      print(prob.index)
      miss.prob = missing.probs[prob.index]
      for (iter in 1:num.iters) {
        file.name = paste0(prob.index, "_", iter)
        res.file.path = file.path(res.dir, file.name)
	if (alg.name == 'pvc') {
	  indices.dir = file.path(get.results.dir(), paste0(alg.name, '_indices'), subtype.dir)
	  if (!is.pvc.results.calculated(indices.dir, subtype.raw.data, prob.index, iter)) next # since pvc does not run on 3 omics
	  cur.clustering = get.pvc.results(indices.dir, subtype.raw.data, prob.index, iter)  
	} else {
	  var.name = load(res.file.path)
	  cur.clustering = get(var.name)
	}
	names(cur.clustering$clustering) = colnames(subtype.raw.data[[1]])
	
	surv.file.path = file.path(get.results.dir(), 'surv', paste(alg.name, subtype.dir, prob.index, iter, sep='_'))
	if (!file.exists(surv.file.path)) {
	  empirical.surv = get.empirical.surv(cur.clustering$clustering, subtype.name)
	  save(empirical.surv, file=surv.file.path)
	} else {
	  # loads empirical surv
	  load(surv.file.path)
	}
	
	clin.file.path = file.path(get.results.dir(), 'clinical', paste(alg.name, subtype.dir, prob.index, iter, sep='_'))
	if (!file.exists(clin.file.path)) {
	  empirical.clin = check.clinical.enrichment(cur.clustering$clustering, subtype.name)
	  save(empirical.clin, file=clin.file.path)
	} else {
	  # loads empirical clin
	  load(clin.file.path)
	}

	num.clusters = max(cur.clustering$clustering)
	timing = cur.clustering$timing
	
	surv.mat[iter, prob.index] = empirical.surv$pvalue
	clin.mat[iter, prob.index] = sum(empirical.clin * length(empirical.clin) < 0.05)
	num.cluster.mat[iter, prob.index] = num.clusters
	timing.mat[iter, prob.index] = timing

      }
    }
    alg.mats = list(surv.mat=surv.mat, clin.mat=clin.mat, num.cluster.mat=num.cluster.mat, timing.mat=timing.mat)
    subtype.results[[alg.name]] = alg.mats
  }
  return(subtype.results)
}

run.missing.real.data <- function() {
  for (i in 1:length(SUBTYPES.DATA)) {
    # we set the seed for each subtype so that we can run them separately
    set.seed(42)
  
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    print(subtype)
    subtype.raw.data = get.raw.data(subtype, 
                                    only.primary=current.subtype.data$only.primary)
    
    for (j in 1:length(subtype.raw.data)) {
      subtype.raw.data[[j]] = subtype.raw.data[[j]][apply(subtype.raw.data[[j]], 1, var) > 0,]
    }
    
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, 
                                           current.subtype.data)
    # rearrange omics
    subtype.raw.data = subtype.raw.data[c(2, 1 ,3)]
    
    num.iters = 5
    missing.probs = get.missing.values.range()
    run.partial.analysis(subtype, subtype.raw.data[1:2], missing.probs, num.iters)
    subtype.with.third = paste(subtype, 'with_third', sep='_')
    run.partial.analysis(subtype.with.third, subtype.raw.data, missing.probs, num.iters)

  }
}

analyze.missing.simulation.data <- function(sim.dir, missing.probs, num.iters, expected.clustering) {
  #alg.names = c('nemo', 'mkl', 'cca', 'pvc')
  alg.names = c('nemo', 'mkl', 'cca', 'cca2', 'nemo_after', 'pvc')
  all.aris = list()
  for (alg.index in 1:length(alg.names)) {
    ari.mat = matrix(NA, ncol=length(missing.probs), nrow=num.iters)
    alg.name = alg.names[alg.index]
    res.dir = file.path(get.results.dir(), paste0(alg.name, '_results'), sim.dir)
    for (prob.index in 1:length(missing.probs)) {
      print(prob.index)
      miss.prob = missing.probs[prob.index]
      for (iter in 1:num.iters) {
        file.name = paste0(prob.index, "_", iter)
        res.file.path = file.path(res.dir, file.name)
        if (alg.name == 'pvc') {
          indices.dir = file.path(get.results.dir(), 'pvc_indices', sim.dir)
	  temp.mock.data = matrix(nrow=1, ncol=length(expected.clustering))
	  colnames(temp.mock.data) = 1:ncol(temp.mock.data)
	  if (!is.pvc.results.calculated(indices.dir, list(temp.mock.data), prob.index, iter)) next
          cur.clustering = get.pvc.results(indices.dir, list(temp.mock.data), prob.index, iter)  
        } else {
	  if (!file.exists(res.file.path)) next # This happens since mkl crashed on two of the simulated datasets
          var.name = load(res.file.path)
    	  cur.clustering = get(var.name)
        }
	ari = adj.rand.index(cur.clustering$clustering, expected.clustering)
	ari.mat[iter, prob.index] = ari
       }
     }
     all.aris[[alg.name]] = ari.mat
   }   
   return(all.aris)
}

create.simulation.results <- function() {
  set.seed(42)
  num.simulation.repeats = 10
  missing.probs = get.simulation.missing.data.range()
  num.iters = 10

  for (i in 1:num.simulation.repeats) {
    run.missing.data.simulation(i)
  }
  
  png(file.path(get.results.dir(), 'simulation_plot.png'), height=600, width=1250)
  par(mfrow=c(1, 2))
  for (sim.group.index in 1:2) {
    sim.name.suffix = c('', '_with_third')[sim.group.index]
    sim.results = NULL
    for (i in 1:num.simulation.repeats) {
      cur.results = analyze.missing.simulation.data(paste0(i, sim.name.suffix), missing.probs, num.iters, c(rep(1, 150), rep(2, 150)))
      if (is.null(sim.results)) {
        sim.results = cur.results
      } else {
        for (j in 1:length(sim.results)) {
          sim.results[[j]] = rbind(sim.results[[j]], cur.results[[j]])
        }
      }
    }
    plot.partial.results(sim.results, missing.probs, 'Adjusted Rand Index', y.min=0.2)
  }
  dev.off()
  par(mfrow=c(1, 1))
  dev.off()
  
  return(sim.results)
}

run.missing.data.simulation <- function(sim.number) {
  # The function runs one simulation, and sim.number is the index of that simulation run.
  
  # generate data
  num.samples = 300
  dim = 10
  first.cluster = rmvnorm(num.samples / 2, c(rep(2, dim / 2), rep(0, dim / 2)))
  second.cluster = rmvnorm(num.samples / 2, rep(0, dim))
  
  cor.matrix = 2 * diag(1, dim, dim)
  first.cluster.first.omic = first.cluster + rmvnorm(num.samples / 2, rep(0, dim), cor.matrix)
  second.cluster.first.omic = second.cluster + rmvnorm(num.samples / 2, rep(0, dim), cor.matrix)
  
  
  first.cluster.second.omic = first.cluster + rmvnorm(num.samples / 2, rep(0, dim), cor.matrix)
  second.cluster.second.omic = second.cluster + rmvnorm(num.samples / 2, rep(0, dim), cor.matrix)
  
  # third omic with no signal at all
  first.cluster.third.omic = rmvnorm(num.samples / 2, rep(0, dim), 5 * cor.matrix)
  second.cluster.third.omic = rmvnorm(num.samples / 2, rep(0, dim), 5 * cor.matrix)
  data = list(t(rbind(first.cluster.first.omic, second.cluster.first.omic)),
              t(rbind(first.cluster.second.omic, second.cluster.second.omic)),
              t(rbind(first.cluster.third.omic, second.cluster.third.omic)))
  
  # transform data to be non negative.
  for (i in 1:length(data)) {
    data[[i]] = data[[i]] - min(data[[i]])
  }
  
  expected = c(rep(1, num.samples / 2), rep(2, num.samples / 2))
  
  missing.probs = get.simulation.missing.data.range()
  num.iters = 10
  
  for (i in 1:length(data)) {
    data[[i]] = as.data.frame(data[[i]])
    colnames(data[[i]]) = paste0('patient', 1:num.samples)
    attr(data[[i]], 'is.seq') = F
  }
  
  sim.number = as.character(sim.number)
  
  run.partial.analysis(sim.number, data[1:2], missing.probs, num.iters, num.clusters=2, mcca.penalty=sqrt(dim), impute.transpose=F)
  run.partial.analysis(paste(sim.number, 'with_third', sep='_'), data, missing.probs, num.iters, num.clusters=2, mcca.penalty=sqrt(dim), impute.transpose=F)

}

partial.analysis.one.dataset <- function(subtype.name, keep.sample, omics.list, prob.index, iter, orig.is.seq, num.clusters=NULL, mcca.penalty=NULL,
                                 impute.transpose=T) {

  num.samples = ncol(omics.list[[1]])
  results.dir.name = get.results.dir()
  indices.dir = file.path(results.dir.name, 'pvc_indices', subtype.name)
  current.data = omics.list
  current.data[[2]] = omics.list[[2]][,keep.sample]
  for (i in 1:length(current.data)) {
    attr(current.data[[i]], 'is.seq') = orig.is.seq[i]
  }

  mkl.dir = file.path(results.dir.name, 'mkl_results', subtype.name)
  dir.create(mkl.dir)
  cca.dir = file.path(results.dir.name, 'cca_results', subtype.name)
  dir.create(cca.dir)
  cca2.dir = file.path(results.dir.name, 'cca2_results', subtype.name)
  dir.create(cca2.dir)
  nemo.dir = file.path(results.dir.name, 'nemo_results', subtype.name)
  dir.create(nemo.dir)
  nemo.after.dir = file.path(results.dir.name, 'nemo_after_results', subtype.name)
  dir.create(nemo.after.dir)
  
  dir.create(indices.dir)
  file.name = paste0(prob.index, "_", iter)
  
  cca.results.path = file.path(cca.dir, file.name)
  cca2.results.path = file.path(cca2.dir, file.name)
  mkl.results.path = file.path(mkl.dir, file.name)
  nemo.results.path = file.path(nemo.dir, file.name)
  nemo.after.imp.results.path = file.path(nemo.after.dir, file.name)
  
  nemo.data = current.data
  if (!file.exists(nemo.results.path)) {
    nemo.clustering = run.nemo(nemo.data, list(name=subtype.name), num.clusters)
    save(nemo.clustering, file=nemo.results.path)
  } else {
    # loads nemo.clustering variable
    load(nemo.results.path)
  }
  nemo.num.clusters = max(nemo.clustering$clustering)
  
  write.partial.index.file(indices.dir, paste(file.name, nemo.num.clusters, sep='_'),
                           (1:num.samples)[keep.sample])
  
  data.before.imputation = omics.list
  data.before.imputation[[2]][,!keep.sample] = NA
  
  imputation.start.time = Sys.time()
  all.num.features = sapply(data.before.imputation, nrow)
  
  if (length(data.before.imputation) == 3) {
    
    if (impute.transpose) {
      imputed.data.binded = t(knn.impute(t(as.matrix(rbind(data.before.imputation[[1]], data.before.imputation[[2]], data.before.imputation[[3]]))), cat.var=c()))
    } else {
      imputed.data.binded = knn.impute(as.matrix(rbind(data.before.imputation[[1]], data.before.imputation[[2]], data.before.imputation[[3]])), cat.var=c())
    }
    imputed.data = list(imputed.data.binded[1:all.num.features[1],], imputed.data.binded[(all.num.features[1] + 1):(all.num.features[1] + all.num.features[2]),], 
                    imputed.data.binded[(all.num.features[1] + all.num.features[2] + 1):sum(all.num.features),])
  } else {
    if (impute.transpose) {
      imputed.data.binded = t(knn.impute(t(as.matrix(rbind(data.before.imputation[[1]], data.before.imputation[[2]]))), cat.var=c()))
    } else {
      imputed.data.binded = knn.impute(as.matrix(rbind(data.before.imputation[[1]], data.before.imputation[[2]])), cat.var=c())
    }
    imputed.data = list(imputed.data.binded[1:all.num.features[1],], imputed.data.binded[(all.num.features[1] + 1):(all.num.features[1] + all.num.features[2]),])
  }
  
  imputed.data = lapply(imputed.data, as.data.frame)
  imputation.time.taken = as.numeric(Sys.time() - imputation.start.time, units='secs')
  imputation.dir = file.path(results.dir.name, 'imputation_times', subtype.name)
  dir.create(imputation.dir)
  save(imputation.time.taken, file=file.path(imputation.dir, file.name))
  
  for (i in 1:length(imputed.data)) {
    attr(imputed.data[[i]], 'is.seq') = orig.is.seq[i]
  }
  
  if (!file.exists(cca.results.path)) {
    # get solution for MCCA after imputation
    cca.clustering = run.mcca(imputed.data, list(name=subtype.name), num.clusters, rep(mcca.penalty, length(imputed.data)))
    save(cca.clustering, file=cca.results.path)
  } else {
    # loads cca.clustering variable
    load(cca.results.path)
  }
  
  if (!file.exists(cca2.results.path)) {
    cca2.clustering = run.mcca(imputed.data, list(name=subtype.name), num.clusters, rep(mcca.penalty, length(imputed.data)), rep.omic=2)
    save(cca2.clustering, file=cca2.results.path)
  } else {
    # loads cca2.clustering variable
    load(cca2.results.path)
  }
  
  if (!file.exists(nemo.after.imp.results.path)) {
    nemo.after.clustering = run.nemo(imputed.data, list(name=subtype.name), num.clusters)
    save(nemo.after.clustering, file=nemo.after.imp.results.path)
  } else {
    # loads nemo.after.clustering variable
    load(nemo.after.imp.results.path)
  }
        
  # get solution for mkl after imputation
  if (!file.exists(mkl.results.path)) {
    mkl.clustering = run.mkl(imputed.data, list(name=subtype.name), num.clusters, run.id=file.name)
    save(mkl.clustering, file=mkl.results.path)
  } else {
    # loads mkl.clustering variable
    load(mkl.results.path)
  }
  
  # Now run the pvc command line.
  if (!is.pvc.results.calculated(indices.dir, current.data, prob.index, iter) & (length(data.before.imputation) == 2)) {
    system.cmd = paste0('matlab -nodisplay -wait -nosplash -nodesktop -r "expected_prefix=\'', file.name, '_\';subtype_name=string(\'', subtype.name, '\');run(\'', pvc.script.path(), '\');exit;"')
    command.ret = system(system.cmd)
    stopifnot(command.ret == 0)
  }
}

run.partial.analysis <- function(subtype.name, omics.list, missing.probs, num.iters, num.clusters=NULL, mcca.penalty=NULL,
                                 impute.transpose=T) {
  results.dir.name = get.results.dir()
  indices.dir = file.path(results.dir.name, 'pvc_indices', subtype.name)
  num.samples = ncol(omics.list[[1]])
  results.dir.name = get.results.dir()
  # save the generated data, to use the same data in pvc.
  dir.name = file.path(results.dir.name, subtype.name)
  dir.create(dir.name)
  
  
  orig.is.seq = sapply(omics.list, function(omic) attr(omic, 'is.seq'))
  
  omics.for.pvc = log.and.normalize(omics.list, 'UNUSED', normalize=F, filter.var=T)
  for (i in 1:2) {
    save.matlab.format(omics.for.pvc[[i]], dir.name = dir.name, file.name = i)
  }
  
  indices = vector('list', 4)
  for (i in 1:length(indices)) {
    indices[[i]] = matrix(0, num.iters, length(missing.probs))
  }
  
  all.keep.samples = list()
  for (prob.index in 1:length(missing.probs)) {
    miss.prob = missing.probs[prob.index]
    for (iter in 1:num.iters) {
      
      # simulate data loss
      rand.permutation = sample(1:num.samples, num.samples)
      keep.sample = rep(T, num.samples)
      if (miss.prob != 0) {
        keep.sample[1:round(miss.prob*num.samples)] = F  
      }
      keep.sample = keep.sample[rand.permutation]
      
      all.keep.samples[[paste(prob.index, iter, sep='_')]] = keep.sample
    }
  }
  #lapply(names(all.keep.samples), function(dataset.id) partial.analysis.one.dataset(subtype.name, all.keep.samples[[dataset.id]], omics.list, strsplit(dataset.id, '_')[[1]][1], 
  #       strsplit(dataset.id, '_')[[1]][2], orig.is.seq, num.clusters, mcca.penalty, impute.transpose))

  mclapply(names(all.keep.samples), function(dataset.id) partial.analysis.one.dataset(subtype.name, all.keep.samples[[dataset.id]], omics.list, strsplit(dataset.id, '_')[[1]][1], 
         strsplit(dataset.id, '_')[[1]][2], orig.is.seq, num.clusters, mcca.penalty, impute.transpose), mc.cores=60)
}

create.subtype.analysis.results <- function() {
  results.dir.name = file.path(get.results.dir(), 'aml_results')
  aml.all.data = get.raw.data('aml', intersect.patients=F, only.primary = F)
  
  patient.original.names = Reduce(union, lapply(aml.all.data, colnames))
  aml.data.with.fixed.names = fix.patient.names(aml.all.data)
  
  attr(aml.data.with.fixed.names[[1]], 'is.seq') = T
  attr(aml.data.with.fixed.names[[2]], 'is.seq') = F
  attr(aml.data.with.fixed.names[[3]], 'is.seq') = T
  
  aml.sub.omics = list(mRNA=1, methylation=2, miRNA=3, multiomics=c(1, 2, 3))
  for (i in 1:length(aml.sub.omics)) {
    aml.sub.omic.name = names(aml.sub.omics)[[i]]
    aml.sub.omic = aml.sub.omics[[i]]
    print(aml.sub.omic)
    clustering = run.nemo(aml.data.with.fixed.names[aml.sub.omic], 'UNUSED', is.missing.data=T)$clustering
    print(get.empirical.surv(clustering, 'aml'))
    save(clustering, file=file.path(results.dir.name, paste0(aml.sub.omic.name, '_clustering')))
    
    if (length(aml.sub.omic) > 1) {
      write.promo.cluster.file(patient.original.names, clustering,
                               file.path(results.dir.name, 'fab_clustering'))
      
      clinical.params = get.clinical.params('aml')
      fab.calls = clinical.params[names(clustering), 'leukemia_french_american_british_morphology_code']
      names(fab.calls) = names(clustering)
      print(get.empirical.surv(as.factor(fab.calls), 'aml'))
      
      print(table(clustering, fab.calls))
    }
  }
}

write.promo.cluster.file <- function(patient.names, clustering, file.path) {
  patient.clustering = rep(0, length(patient.names))
  for (i in 1:length(patient.clustering)) {
    cluster = clustering[startsWith(names(clustering), 
                                    substr(patient.names[i], 1, 12))][1]
    patient.clustering[i] = cluster
  }
  
  perm = sort.list(patient.clustering)
  patient.names = patient.names[perm]
  patient.names = gsub('\\.', '-', patient.names)
  patient.clustering = patient.clustering[perm]
  promo.clusters = data.frame(patient.names, patient.clustering)
  write.table(promo.clusters, file=file.path,
              row.names = F, col.names = F, quote = F, sep = '\t', eol='\n')
}

##########################
###### Code for PVC ######
##########################

save.matlab.format <- function(mat, dir.name, file.name) {
  full.file.name = file.path(dir.name, file.name)
  write.table(mat, full.file.name)
}

get.pvc.results <- function(dir.name, subtype.raw.data, prob.index, iter) {
  file.prefix = paste0(prob.index, "_", iter, '_')
  file.names = list.files(dir.name)
  
  pvalues = c()
  for (file.name in file.names) {
    if (startsWith(file.name, file.prefix) & endsWith(file.name, 'clustering')) {
      full.file.name = file.path(dir.name, file.name)
      clustering = as.numeric(read.table(full.file.name, sep=','))
      names(clustering) = colnames(subtype.raw.data[[1]])  
      timing = as.numeric(strsplit(file.name, '_')[[1]][4])
      return(list(clustering=clustering, timing=timing))
    }
  }
}

is.pvc.results.calculated <- function(dir.name, subtype.raw.data, prob.index, iter) {
  file.prefix = paste0(prob.index, "_", iter, '_')
  file.names = list.files(dir.name)
  
  pvalues = c()
  for (file.name in file.names) {
    if (startsWith(file.name, file.prefix) & endsWith(file.name, 'clustering')) {
      return(T)
    }
  }
  return(F)
}

write.partial.index.file <- function(dir.name, file.name, indices) {
  missing.samples.path = file.path(dir.name, file.name)
  write.table(indices, file=missing.samples.path,
              col.names = F, row.names = F)
}

pvc.script.path <- function() {
  return('PATH_TO_PVC/PVC_script.m')
}

##########################
####### Plots code #######
##########################

plot.robustness.results <- function(robustness.results, y.max, ticks.breaks, y.lab, plot.name, argum.range, x.title, full.data.results) {
  png(plot.name, width=1500, height=750)
  
  subtype.plots = list()
  plot.titles = sapply(1:length(SUBTYPES.DATA), function(i) SUBTYPES.DATA[[i]]$display.name)
  for (i in 1:length(SUBTYPES.DATA)) {
  
    param = full.data.results$params[i]
    param.value = full.data.results$values[i]
    current.subtype.nemo.values = robustness.results[i,]
    
    # insert in an ordered manner
    if (!(param %in% argum.range)) {
      index.to.insert = max(which(argum.range<param))
      cur.argum.range = append(argum.range, param, index.to.insert)
      current.subtype.nemo.values = append(current.subtype.nemo.values, param.value, index.to.insert)
    } else {
      cur.argum.range = argum.range
    }
    
    
    results.df = data.frame(xrange=cur.argum.range, 
                            nemo.values = current.subtype.nemo.values)
    size.vector = results.df$xrange == param

   ylim.max = ceiling(max(results.df$nemo.values, y.max))
    print(ticks.breaks)
    
    if (i == 1) {
      
      subtype.plot = (ggplot(data=results.df, aes(x=xrange)) + 
                        labs(title=plot.titles[[i]], x=x.title, y=y.lab) + 
                        geom_line(aes(y=nemo.values, color='red'), size=1) +
			#geom_point(data=results.df[which(cur.argum.range == param),], aes(y=force(param.value), size=1.5)) + 
			#geom_point(aes(x=xrange, y=nemo.values, size=ifelse(size.vector, 1.5, 0))) + 
                        scale_y_continuous(breaks=seq(0, ylim.max, by=ticks.breaks), limits=c(0, ylim.max)) +
                        theme(aspect.ratio = 1/1, plot.title = element_text(size=14),
                              
                              axis.title.x = element_text(size=12, color='black'),
                              axis.title.y = element_text(size=12, color='black'),
                              axis.text.x = element_text(size=10, color='black'),
                              axis.text.y = element_text(size=10, color='black'),
                              axis.line.x = element_line(color='black'),
                              axis.line.y = element_line(color='black'),
                              
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              legend.title = element_blank(),
                              legend.background = element_rect(fill='transparent', color = 'transparent'),
                              legend.box.background = element_rect(fill='transparent', linetype = 'blank'),
                              legend.key = element_blank(),
                              legend.text = element_text(size=14),
                              legend.position="none"
                        ))
    } else {
      subtype.plot = (ggplot(data=results.df, aes(x=xrange)) + 
                        labs(title=plot.titles[[i]], x='', y='') + 
                        geom_line(aes(y=nemo.values, color='red'), size=1) +
			#geom_point(data=results.df[which(cur.argum.range == param),], aes(y=param.value, size=1.5)) + 
			#geom_point(aes(x=xrange, y=nemo.values, size=ifelse(size.vector, 1.5, 0))) + 
                        scale_y_continuous(breaks=seq(0, ylim.max, by=ticks.breaks), limits=c(0, ylim.max)) +
                        theme(legend.position="none", aspect.ratio = 1/1, plot.title = element_text(size=14),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              axis.title.x = element_text(size=12, color='black'),
                              axis.title.y = element_text(size=12, color='black'),
                              axis.text.x = element_text(size=10, color='black'),
                              axis.text.y = element_text(size=10, color='black'),
                              axis.line.x = element_line(color='black'),
                              axis.line.y = element_line(color='black'),))
    }
    
    subtype.plots[[length(subtype.plots) + 1]] = subtype.plot
  }
  print(grid.arrange(subtype.plots[[1]] + theme(legend.position = 'none'), subtype.plots[[2]], subtype.plots[[3]],
                     subtype.plots[[4]], subtype.plots[[5]], subtype.plots[[6]], subtype.plots[[7]], subtype.plots[[8]], subtype.plots[[9]], subtype.plots[[10]], nrow=2))
  
  dev.off()
}


plot.missing.real.data.results <- function(all.partial.results) {
  plot.path = get.results.dir()
  tiff(paste(plot.path, 'missing_real_data_results', '.tiff', sep=''), width=1200, height=300)
  
  palette = c("#000000",  "#56B4E9", "#E69F00")
  values.list = palette
  names(values.list) = c('NEMO', 'PVC', 'significance')
  
  subtype.plots = list()
  plot.titles = names(all.partial.results)
  for (i in 1:length(all.partial.results)) {
    current.subtype.results = all.partial.results[[i]]
    current.subtype.num.clusters = current.subtype.results[[1]]
    current.subtype.nemo.pvalues = colMeans(-log10(current.subtype.results[[2]]))
    current.subtype.pvc.pvalues = colMeans(-log10(current.subtype.results[[3]]))
    results.df = data.frame(range=get.missing.values.range(), 
                            nemo.pvalues = (current.subtype.nemo.pvalues),
                            pvc.pvalues = (current.subtype.pvc.pvalues), 
                            sig.threshold=rep(-log10(0.05), length(current.subtype.nemo.pvalues)))
    ylim.max = ceiling(2 * max(results.df$nemo.pvalues, 
                           results.df$pvc.pvalues,
                           -log10(0.05))) / 2
    ylim.min = floor(2 * min(results.df$nemo.pvalues, 
                           results.df$pvc.pvalues,
                           -log10(0.05))) / 2
    
    if (i == 5) {
      ticks.break = 1
    } else {
      ticks.break = 0.5
    }
    
    if (i == 1) {
      subtype.plot = (ggplot(data=results.df, aes(x=range)) + 
                        labs(title=plot.titles[[i]], x='missing data fraction', y='-log10 logrank pvalue') + 
                        geom_line(aes(y=nemo.pvalues, color='NEMO', linetype='NEMO'), size=1) +
                        geom_line(aes(y=pvc.pvalues, color='PVC', linetype='PVC'), size=1) + 
                        geom_line(aes(y=sig.threshold, color='significance', linetype='significane'), size=1) + 
                        scale_y_continuous(breaks=seq(ylim.min, ylim.max, by=ticks.break), limits=c(ylim.min, ylim.max)) +
                        theme(legend.position='none', aspect.ratio = 1/1, plot.title = element_text(size=10),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                              panel.background = element_blank(),
                              axis.title.x = element_text(size=12, color='black'),
                              axis.title.y = element_text(size=12, color='black'),
                              axis.text.x = element_text(size=10, color='black'),
                              axis.text.y = element_text(size=10, color='black'),
                              axis.line.x = element_line(color='black'),
                              axis.line.y = element_line(color='black'),
                        ) + scale_color_manual(values=values.list) +
                        scale_linetype_manual(values=c('solid', 'dashed', 'dotted')))
    } else {
      if (i == 2) {
        plot.theme = theme(aspect.ratio = 1/1, plot.title = element_text(size=10),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              
              axis.title.x = element_text(size=12, color='black'),
              axis.title.y = element_text(size=12, color='black'),
              axis.text.x = element_text(size=10, color='black'),
              axis.text.y = element_text(size=10, color='black'),
              axis.line.x = element_line(color='black'),
              axis.line.y = element_line(color='black'),
              
              legend.title = element_blank(),
              legend.position = c(0.4, -0.15),
              legend.direction = 'horizontal',
              legend.background = element_rect(fill='transparent', color = 'transparent'),
              legend.box.background = element_rect(fill='transparent', linetype = 'blank'),
              legend.key = element_blank(),
              legend.text = element_text(size=10))
      } else {
        plot.theme = theme(legend.position="none", aspect.ratio = 1/1, plot.title = element_text(size=10),
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(),
                           
                           axis.title.x = element_text(size=12, color='black'),
                           axis.title.y = element_text(size=12, color='black'),
                           axis.text.x = element_text(size=10, color='black'),
                           axis.text.y = element_text(size=10, color='black'),
                           axis.line.x = element_line(color='black'),
                           axis.line.y = element_line(color='black'))
      }
      subtype.plot = (ggplot(data=results.df, aes(x=range)) + 
                        labs(title=plot.titles[[i]], x='', y='') + 
                        geom_line(aes(y=nemo.pvalues, color='NEMO', linetype='NEMO'), size=1) +
                        geom_line(aes(y=pvc.pvalues, color='PVC', linetype='PVC'), size=1) + 
                        geom_line(aes(y=sig.threshold, color='significance', linetype='significance'), size=1) + 
                        scale_y_continuous(breaks=seq(ylim.min, ylim.max, by=ticks.break), limits=c(ylim.min, ylim.max)) +
                        plot.theme +
                        scale_color_manual(values=values.list, labels=names(values.list), name='') +
                        scale_linetype_manual(values=c('solid', 'dashed', 'dotted'), labels=names(values.list), name='')+
                        guides(color = guide_legend(override.aes=list(linetype=c('solid', 'dashed', 'dotted')))))
    }
    
    subtype.plots[[length(subtype.plots) + 1]] = subtype.plot
  }
  plot.legend = g_legend(subtype.plots[[2]])
  print(grid.arrange(subtype.plots[[1]], subtype.plots[[2]], subtype.plots[[3]],
                     subtype.plots[[4]], plot.legend, nrow=1))
  dev.off()
}

plot.partial.results <- function(simulation.results, x.values, ylabel, y.min=0, main='', include.labels=T) {
  result.colors = c('black', 'brown', 'blue', 'purple', 'red', 'forestgreen')
  lines.styles = c(1, 1, 2, 2, 3, 3)
  result.names = c('nemo', 'nemo_after', 'cca', 'cca2', 'mkl', 'pvc')
  
  max.lim = max(sapply(simulation.results, function(x) max(x, na.rm=T)))
  
  if (include.labels) {
    xlab.val = 'Fraction of missing data'
    ylab.val = ylabel
  } else {
    xlab.val = ''
    ylab.val = ''
  }
  
  if (is.vector(simulation.results[[result.names[1]]])) {
    func = identity
  } else {
    func = function(x) colMeans(x, na.rm=T)
  }
  
  plot(x.values, func(simulation.results[[result.names[1]]]), col=result.colors[1], type='l', ylim=c(y.min, max.lim),
       xlab=xlab.val, ylab=ylab.val, lty=lines.styles[1], lwd=2, cex.lab=2.4, cex.axis=2, main=main)
  for (i in 2:length(result.names)) {
    lines(x.values, func(simulation.results[[result.names[i]]]), col=result.colors[i], type='l', lty=lines.styles[i], lwd=2)
  }

  #legend(0.05, 0.3, legend=c('NEMO without imputation', 'NEMO with imputation', 'MCCA using first omic representation', 'MCCA using second omic representation', 'rMKL-LPP', 'PVC'),
  #            col=result.colors, box.lty=0, lty=lines.styles, lwd=2, cex=1.1, horiz=T)
}

run.handwritten.data <- function() {
  set.seed(42)
  handwritten.data.and.labs = get.handwritten()
  handwritten.data = handwritten.data.and.labs[[1]]
  handwritten.labels = handwritten.data.and.labs[[2]]
  handwritten.name = 'handwritten'
  
  all.sel.samples = c()
  all.sel.samples = sample(1:2000, 500)
  handwritten.labels = handwritten.labels[all.sel.samples]
  handwritten.data = lapply(handwritten.data, function(x) x[,all.sel.samples])
  
  attr(handwritten.data[[1]], 'is.seq') = F
  attr(handwritten.data[[2]], 'is.seq') = F
  
  num.iters = 10
  missing.probs = c(0, 0.5)
  set.seed(42)
  run.partial.analysis(handwritten.name, handwritten.data, missing.probs, num.iters, num.clusters=10, impute.transpose=F)
  
  # and now analyze the results
  handwritten.analysis = analyze.missing.simulation.data(handwritten.name, missing.probs, num.iters, handwritten.labels)
  print(handwritten.analysis)
  return(handwritten.analysis)
}

get.handwritten <- function() {
  file1 = file.path("path_to_pixel_omic")
  file2 = file.path("path_to_fourier_omic")
  all.data = lapply(c(file1, file2), function(fpath) {
    data.mat = read.table(fpath)
    rownames(data.mat) = paste0('web', 1:nrow(data.mat))
    colnames(data.mat) = paste0('feat', 1:ncol(data.mat))
    ret = as.data.frame(t(data.mat))
    attr(ret, 'is.seq') = F
    return(ret)
  })
  labs = unlist(lapply(1:10, function(x) rep(x, 200)))
  return(list(all.data, labs))
}
