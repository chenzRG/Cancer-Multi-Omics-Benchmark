SUBTYPES.DATA = list(
  list(name='aml', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='AML'),
  list(name='breast', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='BIC'),
  list(name='colon', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='COAD'),
  list(name='gbm', only.primary=T, is.rna.seq=F, is.mirna.seq=F, display.name='GBM'),
  list(name='kidney', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='KIRC'),
  list(name='liver', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='LIHC'),
  list(name='lung', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='LUSC'),
  list(name='melanoma', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='SKCM'),
  list(name='ovarian', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='OV'),
  list(name='sarcoma', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SARC'))

MAX.NUM.CLUSTERS = 15

OMIC.SUBSETS = list('multi_omics', 'exp', 'methy', 'mirna')
names(OMIC.SUBSETS) = c('all', '1', '2', '3')

get.clustering.results.dir.path <- function() {
  return('RESULTS_PATH')
}

get.plots.dir.path <- function() {
  results.dir.path = get.clustering.results.dir.path()
  return(file.path(results.dir.path, 'plots'))
}

get.tables.dir.path <- function() {
  results.dir.path = get.clustering.results.dir.path()
  return(file.path(results.dir.path, 'tables'))
}

subtype.to.display.name <- function(subtype) {
  for (i in 1:length(SUBTYPES.DATA)) {
    if (SUBTYPES.DATA[[i]]$name == subtype) {
      return(SUBTYPES.DATA[[i]]$display.name)
    }
  }
}

run.nemo <- function(omics.list, subtype.data, num.clusters=NULL, is.missing.data=F, num.neighbors=NA) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
  clustering = nemo.clustering(omics.list, num.clusters, num.neighbors)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

set.omics.list.attr <- function(subtype.raw.data, subtype.data) {
  attr(subtype.raw.data[[1]], 'is.seq') = subtype.data$is.rna.seq
  attr(subtype.raw.data[[2]], 'is.seq') = F
  attr(subtype.raw.data[[3]], 'is.seq') = subtype.data$is.mirna.seq
  return(subtype.raw.data)
}

ALGORITHM.NAMES = c('kmeans', 'spectral', 'lracluster', 'pins', 'snf', 'mkl', 
                     'mcca', 'nmf', 'iCluster', 'nemo')
ALGORITHM.DISPLAY.NAMES = as.list(c('K-means', 'Spectral', 'LRAcluster', 'PINS', 
                           'SNF', 'rMKL-LPP', 'MCCA', 'MultiNMF', 'iClusterBayes', 'NEMO'))
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES
			
print.matrix.latex.format <- function(mat) {
  print(do.call(paste, as.list(c(colnames(mat), sep=' & '))))
  for (i in 1:nrow(mat)) {
    print(do.call(paste, as.list(c(rownames(mat)[i], round(mat[i,], digits=2), sep=' & '))))
  }
}


perform.all.analyses <- function(benchmark.ret) {

  for (i in 1:4) {
    
    cur.func = list(benchmark.omics.time, benchmark.omics.num.clusters,
      benchmark.omics.surv, benchmark.omics.clinical)[[i]]
    for (omic.subset in names(OMIC.SUBSETS)) {
      print(paste('current omic subset ', omic.subset))
      benchmark.data = cur.func(benchmark.ret, omic.subset)
      
      displayed.benchmark.data = benchmark.data
      colnames(displayed.benchmark.data)[1:ncol(displayed.benchmark.data)] = 
        sapply(as.list(colnames(displayed.benchmark.data)[1:ncol(displayed.benchmark.data)]), 
	subtype.to.display.name)
      rownames(displayed.benchmark.data) = unlist(ALGORITHM.DISPLAY.NAMES[rownames(displayed.benchmark.data)])
      print.matrix.latex.format(displayed.benchmark.data)
      
      table.name = c('runtime', 'num_cluster', 'survival', 'clinical')[i]
      write.csv(displayed.benchmark.data, file=file.path(get.tables.dir.path(), paste0(table.name, '_', OMIC.SUBSETS[[omic.subset]], '.csv')))
    }
    
    
    print('------------------------')
  }
  
  
  # plots that include all datasets
  for (i in 1:4) {
    omic.subset = names(OMIC.SUBSETS)[[i]]
    benchmark.surv = benchmark.omics.surv(benchmark.ret, omic.subset)
    benchmark.clinical = benchmark.omics.clinical(benchmark.ret, omic.subset)
    plot.name = list('multi_omics_surv_clinical.png', 'exp_surv_clinical.png', 
      'methy_surv_clinical.png', 'mirna_surv_clinical.png')[[i]]
    create.clinical.survival.plots(benchmark.surv, benchmark.clinical, plot.name)
  }
  
  
  # plots for the mean behaviour
  create.mean.clinical.survival.plot(benchmark.ret)
}

get.best.single.omic.mat <- function(single.omic.list1, single.omic.list2) {
  best.mat1 = matrix(0, ncol=ncol(single.omic.list1[[1]]), 
                            nrow=nrow(single.omic.list1[[1]]))
  best.mat2 = matrix(0, ncol=ncol(single.omic.list1[[1]]), 
                            nrow=nrow(single.omic.list1[[1]]))
  for (i in 1:nrow(best.mat1)) {
    for (j in 1:ncol(best.mat1)) {
      values1 = sapply(single.omic.list1, function(mat) mat[i, j])
      values2 = sapply(single.omic.list2, function(mat) mat[i, j])
      if (any(is.na(values1))) {
        best.mat1[i, j] = NA
	best.mat2[i, j] = NA
      } else {
        best.omic = which.max(values1)
        best.mat1[i, j] = values1[best.omic]
	best.mat2[i, j] = values2[best.omic]
      }
    }
  }
  return(list(best.mat1, best.mat2))
}

create.mean.clinical.survival.plot <- function(benchmark.ret) {
  png(file.path(get.plots.dir.path(), 'mean_surv_clinical.png'), width=1500, height=750)
  layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5), nrow=2, byrow=T))
  single.omic.surv.benchmarks = list()
  single.omic.clin.benchmarks = list()
  for (i in 2:6) {
    omic.subset = c(names(OMIC.SUBSETS), 'best_surv', 'best_clin')[i]
    alg.cols = c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#842121")
    pch = c(rep(15, 3), rep(3, 3), rep(19, 3))
    
    if (i <= 4) {
      benchmark.surv = benchmark.omics.surv(benchmark.ret, omic.subset)
      benchmark.clinical = benchmark.omics.clinical(benchmark.ret, omic.subset)
    
      if (i %in% 2:4) {
        single.omic.surv.benchmarks[[i - 1]] = benchmark.surv
        single.omic.clin.benchmarks[[i - 1]] = benchmark.clinical
      }
    
    } else if (i == 5) {
      best.mats = get.best.single.omic.mat(single.omic.surv.benchmarks, single.omic.clin.benchmarks)
      benchmark.surv = best.mats[[1]]
      benchmark.clinical = best.mats[[2]]

    } else {
      best.mats = get.best.single.omic.mat(single.omic.clin.benchmarks, single.omic.surv.benchmarks)
      benchmark.surv = best.mats[[2]]
      benchmark.clinical = best.mats[[1]]
    }
    
    surv.sum = rowSums(benchmark.surv, na.rm=T)
    clin.sum = rowSums(benchmark.clinical, na.rm=T)
    
    # for mcca, the sum is for an empty vector, and is equal 0, so we remove this
    available.indices = surv.sum != 0
    surv.sum = surv.sum[available.indices]
    clin.sum = clin.sum[available.indices]
    alg.cols = alg.cols[available.indices]
    pch = pch[available.indices]
    
    
    if (omic.subset == '1') {
      xlab = '-log10(logrank pvalue)'
      ylab = '# enriched clinical parameters'
    } else {
      xlab = ''
      ylab = ''
    }
    
    if (i %in% c(1, 5, 6)) {
      y.min = min(clin.sum) - 1
      x.min = min(surv.sum) - 1
    } else {
      y.min = 0
      x.min = 0
    }
    
    
    subplot.name = c('Multi-omics', 'Gene Expression', 'DNA Methylation', 'miRNA Expression',
                     'Best single-omic (survival)', 'Best single-omic (clinical)')[i]
    plot(surv.sum, clin.sum, main=subplot.name, xlab=xlab, ylab=ylab,
         xlim=c(x.min, max(surv.sum) + 1), ylim=c(y.min, max(clin.sum + 1)), col=alg.cols, pch=pch, cex.lab=1.8, cex=4, cex.axis=1.5, cex.main=2.5, lwd=4)
  }
  dev.off()
}

create.clinical.survival.plots <- function(benchmark.surv, benchmark.clinical, plot.name) {
  num.subtypes = ncol(benchmark.surv)
  
  png(file.path(get.plots.dir.path(), plot.name), width=1200)
  alg.cols = c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#842121")
  pch = c(rep(15, 3), rep(3, 3), rep(19, 3))

  par(mfrow=c(2, num.subtypes / 2)) 
  for (i in 1:num.subtypes) {
    subtype = colnames(benchmark.surv)[i]
    subtype.surv = benchmark.surv[,subtype]
    subtype.clinical = benchmark.clinical[,subtype]
    
    available.indices = !is.na(subtype.surv)
    subtype.surv = subtype.surv[available.indices]
    subtype.clinical = subtype.clinical[available.indices]
    current.cols = alg.cols[available.indices]
    current.pch = pch[available.indices]
    
    surv.significance = -log10(0.05)
    
    if (i == 1) {
      xlab = '-log10(logrank pvalue)'
      ylab = '# enriched clinical parameters'
    } else {
      xlab = ''
      ylab = ''
    }
    
    plot(subtype.surv, subtype.clinical, main=subtype.to.display.name(subtype), xlab=xlab, ylab=ylab,
         xlim=c(0, max(subtype.surv, surv.significance) + 0.2), ylim=c(0, max(subtype.clinical + 1)), col=current.cols, pch=pch, cex.lab=1.8, cex=3, cex.axis=1.5, cex.main=2, lwd=3)
    abline(v=surv.significance, col='red')
  }
  
  dev.off()
}

plot.legend <- function() {
  alg.cols = c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#842121")
  pch = c(rep(15, 3), rep(3, 3), rep(19, 3))
  algorithm.names = ALGORITHM.DISPLAY.NAMES
  
  if (T) {
    alg.cols = alg.cols[-7]
    pch = pch[-7]
    algorithm.names = algorithm.names[-7]
  }
  
  png(file.path(get.plots.dir.path(), 'legend2.png'), width=2000, height=400)
  par(font=2)
  plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
  legend(0.7, 1, legend=algorithm.names, col=alg.cols, pch=pch, horiz=T, pt.cex=1.8, cex=1.2, bty='n')
  dev.off()
}

analyze.benchmark <- function() {
  all.clusterings = list()
  all.timings = list()
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, 
                                    only.primary=current.subtype.data$only.primary)
				    
    all.clusterings[[subtype]] = list()
    all.timings[[subtype]] = list()
    
    for (algorithm.name in ALGORITHM.NAMES) {
      all.clusterings[[subtype]][[algorithm.name]] = list()
      all.timings[[subtype]][[algorithm.name]] = list()
      for (j in c('all', '1', '2', '3')) {
        clustering.path = file.path(get.clustering.results.dir.path(),
                                  paste(subtype, algorithm.name, j, sep='_'))
        timing.path = file.path(get.clustering.results.dir.path(),
                                  paste(subtype, algorithm.name, j, 'timing', sep='_'))
	load(clustering.path)
	load(timing.path)
	if (!any(is.na(clustering))) {
	  names(clustering) = colnames(subtype.raw.data[[1]])
	}
	
	all.clusterings[[subtype]][[algorithm.name]][[j]] = clustering
	all.timings[[subtype]][[algorithm.name]][[j]] = timing
      }
    }
  }
  return(list(all.clusterings=all.clusterings, all.timings=all.timings))
}

check.empirical.surv <- function(old.pvals, new.pvals) {
  is.in.conf = matrix(0, ncol=ncol(old.pvals), nrow=nrow(old.pvals))
  for (i in 1:nrow(old.pvals)) {
    for (j in 1:ncol(old.pvals)) {
      old.pval = old.pvals[i, j]
      if (is.na(old.pval)) {
        is.in.conf[i, j] = NA
	next
      }
      if (old.pval == 1e-10) old.pval = 1e-5
      num.runs = floor(30 / old.pval)
      new.pval = new.pvals[i, j]
      num.success = round(new.pval * num.runs)
      print(c(num.success, num.runs))
      
      conf.int = binom.test(num.success, num.runs)$conf.int
      in.conf.int = old.pval >= conf.int[1]  & old.pval <= conf.int[2]
      is.in.conf[i, j] = in.conf.int
    }
  }
  return(is.in.conf)
}

get.empirical.surv <- function(clustering, subtype) {
  set.seed(42)
  surv.ret = check.survival(clustering, subtype)
  orig.chisq = surv.ret$chisq
  orig.pvalue = get.logrank.pvalue(surv.ret)
  # The initial number of permutations to run
  num.perms = round(min(max(10 / orig.pvalue, 1000), 1e6))
  should.continue = T
  
  total.num.perms = 0
  total.num.extreme.chisq = 0
  
  while (should.continue) {
    print('Another iteration in empirical survival calculation')
    print(num.perms)
    perm.chisq = as.numeric(mclapply(1:num.perms, function(i) {
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
      cur.chisq = check.survival(cur.clustering, subtype)$chisq
      return(cur.chisq)
    }, mc.cores=60))
    
    total.num.perms = total.num.perms + num.perms
    total.num.extreme.chisq = total.num.extreme.chisq + sum(perm.chisq >= orig.chisq)
    
    binom.ret = binom.test(total.num.extreme.chisq, total.num.perms)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
    
    print(c(total.num.extreme.chisq, total.num.perms))
    print(cur.pvalue)
    print(cur.conf.int)
    
    sig.threshold = 0.05
    is.conf.small = ((cur.conf.int[2] - cur.pvalue) < min(cur.pvalue / 10, 0.01)) & ((cur.pvalue - cur.conf.int[1]) < min(cur.pvalue / 10, 0.01))
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if ((is.conf.small & !is.threshold.in.conf) | (total.num.perms > 2e7)) {
      should.continue = F
    } else {
      num.perms = 1e5
    }
  }
  
  return(list(pvalue = cur.pvalue, conf.int = cur.conf.int, total.num.perms=total.num.perms, 
              total.num.extreme.chisq=total.num.extreme.chisq))
}

benchmark.omics.surv <- function(benchmark.results, omics='all') {
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      surv.path = file.path(get.clustering.results.dir.path(),
                                  paste(subtype, ALGORITHM.NAMES[j], omics, 'surv', sep='_'))
      if (file.exists(surv.path)) {
        load(surv.path)
	pvalue = empirical.surv.ret$pvalue
      } else {
        clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]][[omics]]
        if (length(table(clustering)) > 1) {
          #pvalue = -log10(get.logrank.pvalue(check.survival(clustering, subtype)))
	  empirical.surv.ret = get.empirical.surv(clustering, subtype)
	  save(empirical.surv.ret, file=surv.path)
	  pvalue = empirical.surv.ret$pvalue

        } else {
          pvalue = NA
        }  
      }
      all.surv.pvalues[j, i] = -log10(pvalue)
    }
  }
  return(all.surv.pvalues)
}

benchmark.omics.num.clusters <- function(benchmark.results, omics='all') {
  num.clusters = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(num.clusters) = ALGORITHM.NAMES
  colnames(num.clusters) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(num.clusters)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]][[omics]]
      num.clusters[j, i] = max(clustering)
    }
  }
  return(num.clusters)
}

benchmark.omics.clinical <- function(benchmark.results, omics='all') {
  num.clinical.enrich = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(num.clinical.enrich) = ALGORITHM.NAMES
  colnames(num.clinical.enrich) = sapply(SUBTYPES.DATA, function(x) x$name)
  total.num.tested.parameters = 0
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(num.clinical.enrich)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      print('checking clinical enrichment')
      print(c(i, j))
      
      
      clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]][[omics]]
      if (any(is.na(clustering))) {
        num.clinical.enrich[j, i] = NA
      } else {
        clin.path = file.path(get.clustering.results.dir.path(),
                                  paste(subtype, ALGORITHM.NAMES[j], omics, 'clin', sep='_'))
				  
        if (file.exists(clin.path)) {
	  # load enrichment.pvalues
	  load(clin.path)
	} else {
	  enrichment.pvalues = check.clinical.enrichment(clustering, subtype)
	  save(enrichment.pvalues, file=clin.path)
	}
        
	if (j == 1) {
	  total.num.tested.parameters = total.num.tested.parameters + length(enrichment.pvalues)
	}
        num.clinical.enrich[j, i] = sum(enrichment.pvalues * length(enrichment.pvalues) < 0.05)
      }
    }
  }
  print(paste0('Total number of parameters tested:', total.num.tested.parameters))
  return(num.clinical.enrich)
}

benchmark.omics.time <- function(benchmark.results, omics='all') {
  all.alg.times = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.alg.times) = ALGORITHM.NAMES
  colnames(all.alg.times) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.timings = benchmark.results$all.timings
  for (i in 1:length(all.timings)) {
    subtype = colnames(all.alg.times)[i]
    subtype.timings = all.timings[[subtype]]
    for (j in 1:length(subtype.timings)) {
      timing = subtype.timings[[ALGORITHM.NAMES[j]]][[omics]]
      all.alg.times[j, i] = timing
    }
  }
  return(all.alg.times)
}

benchmark.alg.agreement <- function(benchmark.results, omics='all') {
  all.clusterings = benchmark.results$all.clusterings
  rand.matrix = matrix(NA, ncol=length(ALGORITHM.NAMES), nrow=length(ALGORITHM.NAMES))
  colnames(rand.matrix) = ALGORITHM.NAMES
  rownames(rand.matrix) = ALGORITHM.NAMES
  for (alg1.index in 1:length(ALGORITHM.NAMES)) {
    alg1 = ALGORITHM.NAMES[alg1.index]
    for (alg2.index in 1:length(ALGORITHM.NAMES)) {
      alg2 = ALGORITHM.NAMES[alg2.index]
      rands = c()
      for (i in 1:length(SUBTYPES.DATA)) {
        clustering1 = all.clusterings[[i]][[alg1]][[omics]]
	clustering2 = all.clusterings[[i]][[alg2]][[omics]]
	rands = c(rands, adj.rand.index(clustering1, clustering2))
      }
      mean.rand = mean(rands)
      rand.matrix[alg1.index, alg2.index] = mean.rand
    }
  }
  return(rand.matrix)
}

get.logrank.pvalue <- function(survdiff.res) {
  1 - pchisq(survdiff.res$chisq, length(survdiff.res$n) - 1)  
}

run.benchmark <- function() {
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, 
                                    only.primary=current.subtype.data$only.primary)
    
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, 
                                           current.subtype.data)
    
    for (algorithm.name in ALGORITHM.NAMES) {
      for (j in c('all', '1', '2', '3')) {
        set.seed(42)
        print(paste('subtype', subtype, 'running algorithm', algorithm.name, j))
        clustering.path = file.path(get.clustering.results.dir.path(),
                                  paste(subtype, algorithm.name, j, sep='_'))
        timing.path = file.path(get.clustering.results.dir.path(),
                                  paste(subtype, algorithm.name, j, 'timing', sep='_'))
	
	
        if (!file.exists(clustering.path)) {
          algorithm.func.name = paste0('run.', algorithm.name)
          algorithm.func = get(algorithm.func.name)
	  if (j == 'all') {
	    cur.iteration.data = subtype.raw.data
	  } else {
	    cur.iteration.data = subtype.raw.data[as.numeric(j)]
	  }
          algorithm.ret = algorithm.func(cur.iteration.data, current.subtype.data)
          clustering = algorithm.ret$clustering
          timing = algorithm.ret$timing
	  print('before saving')
          save(clustering, file = clustering.path)
          save(timing, file = timing.path)
        }
      }
    }
  }
}

get.dataset.dir.path <- function() {
  return('DATASET DIR PATH')
}

log.and.normalize <- function(omics.data, subtype.data, normalize=T,
                              filter.var=F) {
  # filter features with no variance at all
  for (i in 1:length(omics.data)) {
    omics.data[[i]] = omics.data[[i]][apply(omics.data[[i]], 1, var) > 0,]
  }
			      
  for (i in 1:length(omics.data)) {
    if (attr(omics.data[[i]], 'is.seq')) {
      omics.data[[i]] = log(1+omics.data[[i]])
    }
  }
  
  if (filter.var) {
    omics.data = lapply(omics.data, keep.high.var.features)
  }
  
  if (normalize) {
    omics.data = lapply(omics.data, normalize.matrix)    
  }
  
  return(omics.data)
}

normalize.matrix <- function(data.matrix) {
  temp = data.matrix - rowMeans(data.matrix)
  should.keep = (apply(temp, 1, sd) != 0)
  return ((temp / apply(temp, 1, sd))[should.keep, ])
}

filter.non.tumor.samples <- function(raw.datum, only.primary=only.primary) {
  # 01 is primary, 06 is metastatic, 03 is blood derived cancer
  if (!only.primary)
    return(raw.datum[,substring(colnames(raw.datum), 14, 15) %in% c('01', '03', '06')])
  else
    return(raw.datum[,substring(colnames(raw.datum), 14, 15) %in% c('01')])
}

get.fixed.names <- function(patient.names, include.type=F) {
  # fix the TCGA names to only include the patient ids
  if (include.type) {
    return(gsub('-', '\\.', toupper(substring(patient.names, 1, 15))))
  } else {
    return(gsub('-', '\\.', toupper(substring(patient.names, 1, 12))))  
  }
}

fix.patient.names <- function(subtype.raw.data, include.type=F) {
  for (i in 1:length(subtype.raw.data)) {
    colnames(subtype.raw.data[[i]]) = get.fixed.names(colnames(subtype.raw.data[[i]]),
                                                      include.type)
  }
  return(subtype.raw.data)
}

get.raw.data <- function(subtype.name,
                         datasets.path = get.dataset.dir.path(),
                         only.primary=NA,
			 intersect.patients=T) {
  omics.dir = file.path(datasets.path, subtype.name)
  omics.files = list.files(omics.dir)
  omics.files = setdiff(omics.files, c('survival'))  
  raw.data = lapply(file.path(omics.dir, omics.files), read.table)
  
  if (!is.na(only.primary)) {
    raw.data = lapply(raw.data, function(x) filter.non.tumor.samples(x, only.primary = only.primary))
  }
  if (intersect.patients) {
    name.corrected.data = fix.patient.names(raw.data)
    patients.intersection = Reduce(intersect, lapply(name.corrected.data, colnames))
    ret.data = lapply(name.corrected.data, function(datum) datum[,patients.intersection])  
  } else {
    ret.data = raw.data
  }
  
  return(ret.data)
}

get.elbow <- function(values, is.max) {
  second.derivatives = c()
  for (i in 2:(length(values) - 1)) {
    second.derivative = values[i + 1] + values[i - 1] - 2 * values[i]
    second.derivatives = c(second.derivatives, second.derivative)
  }
  print(second.derivatives)
  if (is.max) {
    return(which.max(second.derivatives) + 1)
  } else {
    return(which.min(second.derivatives) + 1)
  }
}

# Does not support a single omic dataset
# Does not support a single omic dataset
run.mcca <- function(omics.list, subtype.data, known.num.clusters=NULL, penalty=NULL, rep.omic=1) {

  if (length(omics.list) == 1) {
    return(list(clustering=rep(1, ncol(omics.list[[1]])), timing=1))
  }
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, 
                               normalize = T,
                               filter.var = T)

  
  max.dim = min(MAX.NUM.CLUSTERS, nrow(omics.list[[1]]))
  
  subtype = subtype.data$name
  omics.transposed = lapply(omics.list, t)
  
  cca.ret = PMA::MultiCCA(omics.transposed, 
                          ncomponents = max.dim, penalty=penalty)
  sample.rep = omics.transposed[[rep.omic]] %*% cca.ret$ws[[rep.omic]]
  
  explained.vars = sapply(1:max.dim, 
              function(i) sum(unlist(apply(sample.rep[1:i,,drop=F], 2, var))))
  
  dimension = get.elbow(explained.vars, is.max=F)
  print(dimension)
  sample.rep = sample.rep[,1:dimension]
  sils = c()
  clustering.per.num.clusters = list()
  if (is.null(known.num.clusters)) {
    for (num.clusters in 2:max.dim) {
      cur.clustering = kmeans(sample.rep, num.clusters, iter.max=100, nstart=30)$cluster  
      sil = get.clustering.silhouette(list(t(sample.rep)), cur.clustering)
      sils = c(sils, sil)
      clustering.per.num.clusters[[num.clusters - 1]] = cur.clustering
    }
    cca.clustering = clustering.per.num.clusters[[which.min(sils)]]
  } else {
    cca.clustering = kmeans(sample.rep, known.num.clusters, iter.max=100, nstart=30)$cluster  
  }
  
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=cca.clustering, timing=time.taken))
}

run.snf <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
  subtype = subtype.data$name
  alpha=0.5
  T.val=30
  num.neighbors = round(ncol(omics.list[[1]]) / 10)
  similarity.data = lapply(omics.list, function(x) {affinityMatrix(dist2(as.matrix(t(x)),as.matrix(t(x))), 
                                                                   num.neighbors, alpha)})
  if (length(similarity.data) == 1) {
    W = similarity.data[[1]]
  } else {
    W = SNF(similarity.data, num.neighbors, T.val)  
  }
  
  num.clusters = estimateNumberOfClustersGivenGraph(W, 2:MAX.NUM.CLUSTERS)[[3]]  
  clustering = spectralClustering(W, num.clusters)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.iCluster <- function(omics.list, subtype.data) {
  omics.list = log.and.normalize(omics.list, subtype.data, normalize = F)

  start = Sys.time()
  subtype = subtype.data$name
  dev.ratios = c()
  icluster.rets = list()

  if (length(omics.list) == 1) {
    icluster.ret = iClusterPlus::tune.iClusterBayes(cpus=(MAX.NUM.CLUSTERS - 1), t(omics.list[[1]]), 
                              K=1:(MAX.NUM.CLUSTERS - 1), type=c('gaussian'))$fit
  } else {
    icluster.ret = iClusterPlus::tune.iClusterBayes(cpus=(MAX.NUM.CLUSTERS - 1), t(omics.list[[1]]), 
                              t(omics.list[[2]]), 
                              t(omics.list[[3]]), 
                              K=1:(MAX.NUM.CLUSTERS - 1), type=rep('gaussian', 3))$fit
  }
  dev.ratios = lapply(1:(MAX.NUM.CLUSTERS - 1), function(i) icluster.ret[[i]]$dev.ratio)

  print('dev.ratios are:')
  print(dev.ratios)
  
  optimal.solution = icluster.ret[[which.max(dev.ratios)]]
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=optimal.solution$clusters, 
              timing=time.taken))
}

get.mkl.binary.path = function() {
  return('MKL_BINARY_PATH')
}

get.mkl.arguments.path = function() {
  return('MKL_ARGS_PATH')
}


run.mkl <- function(omics.list, subtype.data, known.num.clusters=NULL, run.id=NULL) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
  subtype = subtype.data$name
  time.taken = as.numeric(Sys.time() - start, units='secs')
  export.subtype.to.mkl(omics.list, subtype, run.id=run.id)
  
  start = Sys.time()
  bin.path = get.mkl.binary.path()
  subtype.dir = paste0(file.path(get.mkl.arguments.path(), subtype, run.id), '/')
  paste0(subtype.dir, 'kernels')
  command = paste(bin.path, paste0(subtype.dir, 'kernels'),
                  paste0(subtype.dir, 'ids'),
                  paste0(subtype.dir, 'output'), 
                  '9', '5')
  command.return = system(command)
  stopifnot(command.return == 0)
  
  time.taken2 = as.numeric(Sys.time() - start, units='secs')
  clustering = get.mkl.clustering(subtype, known.num.clusters, run.id)
  return(list(clustering=clustering, 
              timing=time.taken + time.taken2))
}

run.nmf <- function(omics.list, subtype.data) {
  total.time.taken = 0
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data,
                                 filter.var = T, normalize = F)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  total.time.taken = total.time.taken + time.taken
  
  save.subtype.matlab.format(omics.list)
  subtype = subtype.data$name
  if (length(omics.list) > 1) {
    command.ret = system('MULTI_NMF_COMMAND')
    stopifnot(command.ret == 0)
    nmf.timing = read.csv('SOME_TEMP_PATH_TIMING', header=F)[1, 1]
    total.time.taken = total.time.taken + nmf.timing
  } else {
   for (k in 1:MAX.NUM.CLUSTERS) {
     start = Sys.time()
     file.name = paste0('SOME_TEMP_PATH/', k, '_consensus')
     nmf.ret = nmf(omics.list[[1]], k, method='lee')
     coef.mat = t(coef(nmf.ret))
     time.taken = as.numeric(Sys.time() - start, units='secs')
     total.time.taken = total.time.taken + time.taken
     write.table(coef.mat, file=file.name, quote=F, row.names=F, col.names=F, sep=',')
   }
  }
  
  explained.vars = c()
  clustering.per.num.clusters = list()
  for (k in 1:MAX.NUM.CLUSTERS) {
    file.name = paste0('SOME_TEMP_PATH/', k, '_consensus')
    consensus.mat = read.csv(file.name, header=F)
    
    start = Sys.time()
    cur.clustering = apply(consensus.mat, 1, which.max)
    explained.var = sum(unlist(apply(consensus.mat, 2, var)))
    explained.vars = c(explained.vars, explained.var)
    clustering.per.num.clusters[[k]] = cur.clustering
    time.taken = as.numeric(Sys.time() - start, units='secs')
    total.time.taken = total.time.taken + time.taken
  }
  
  dimension = get.elbow(explained.vars, is.max=F)
  nmf.clustering = clustering.per.num.clusters[[dimension]]
  return(list(clustering=nmf.clustering, timing=total.time.taken))  
}

run.pins <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, normalize = F)
  subtype = subtype.data$name
  omics.transposed = lapply(omics.list, t)
  if (length(omics.list) == 1) {
    pins.ret = PINSPlus::PerturbationClustering(data=omics.transposed[[1]],
                                            kMax = MAX.NUM.CLUSTERS)
    clustering = pins.ret$cluster
    
  } else {
    pins.ret = PINSPlus::SubtypingOmicsData(dataList=omics.transposed,
                                            kMax = MAX.NUM.CLUSTERS)
    clustering = pins.ret$cluster2
  }
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.lracluster <- function(omics.list, subtype.data) {
  omics.list = log.and.normalize(omics.list, subtype.data, normalize = F)
  
  subtype = subtype.data$name
  start = Sys.time()
  
  dim.range = 1:MAX.NUM.CLUSTERS
  all.clustering.results = list()
  
  omics.matrix.list = lapply(omics.list, as.matrix)
  for (dimension in dim.range) {
    print(paste('running lra cluster for dimension', dimension))
    data.names = c('gene expression', 'methylation', 'miRNA expression')
    clustering.results = LRAcluster(omics.matrix.list, 
                                    rep('gaussian', length(omics.list)), 
                                    dimension=dimension, data.names)
    all.clustering.results[[dimension]] = clustering.results
  }
  explained.var = sapply(all.clustering.results, function(x) x$potential)
  print(explained.var)
  dimension = get.elbow(explained.var, is.max=F)
  print(dimension)
  solution = all.clustering.results[[dimension]]$coordinate
  
  sils = c()
  clustering.per.num.clusters = list()
  for (num.clusters in 2:MAX.NUM.CLUSTERS) {
    print(paste('running kmeans in lra cluster for num clusters', num.clusters))
    cur.clustering = kmeans(t(solution), num.clusters, iter.max=100, nstart=60)$cluster
    sil = get.clustering.silhouette(list(solution), cur.clustering)
    sils = c(sils, sil)
    clustering.per.num.clusters[[num.clusters - 1]] = cur.clustering
  }
  print(sils)
  chosen.clustering = clustering.per.num.clusters[[which.min(sils)]]
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=chosen.clustering, timing=time.taken))
}

run.kmeans <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, 
                                 filter.var = T)
  
  subtype = subtype.data$name
  all.withinss = c()
  all.clusterings = list()
  k.range = 1:MAX.NUM.CLUSTERS
  for (k in k.range) {
    concat.omics = do.call(rbind, omics.list)
    kmeans.ret = kmeans(t(concat.omics), k, iter.max=100, nstart=60)
    all.withinss = c(all.withinss, kmeans.ret$tot.withinss)
    all.clusterings[[k]] = kmeans.ret$cluster
  }
  
  best.k = get.elbow(all.withinss, is.max=T)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=all.clusterings[[best.k]], 
              timing=time.taken))
}

run.spectral <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, 
                                 filter.var = T)
  subtype = subtype.data$name
  concat.omics = do.call(rbind, omics.list)
  
  similarity.data = affinityMatrix(dist2(as.matrix(t(concat.omics)),
                                         as.matrix(t(concat.omics))), 
                                   20, 0.5)
  num.clusters = estimateNumberOfClustersGivenGraph(similarity.data, 
                                      2:MAX.NUM.CLUSTERS)[[3]]  
  clustering = spectralClustering(similarity.data, num.clusters)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

load.libraries <- function() {
  library('PMA')
  library('R.matlab')
  library('SNFtool')
  library('PINSPlus')
  #library('LRAcluster')
  library('kernlab')
  library('survival')
  library('NMF')
  
  # bioconductor packages
  source("https://bioconductor.org/biocLite.R")
  biocLite("impute")
  #biocLite("iClusterPlus")
}



########################################
###############   MKL    ###############
########################################

radial.basis <- function(mat, gamma) {
  if (missing(gamma)) {
    gamma = 1 / (2*nrow(mat)**2)
  }
  npatients = ncol(mat)
  output.mat = matrix(0, ncol=npatients, nrow=npatients)
  for (i in 1:npatients) {
    for (j in 1:npatients) {
      output.mat[i, j] = exp(-norm(as.matrix(mat[,i] - mat[,j]), type = 'F')**2 * gamma)
    }
  }
  
  D = apply(output.mat, 2, sum) / npatients
  E = sum(D) / npatients
  J = matrix(1, nrow=npatients, ncol=1) %*% D
  ret = output.mat - J - t(J) + E * matrix(1, ncol=npatients, nrow=npatients)
  ret = diag(1/sqrt(diag(ret))) %*% ret %*% diag(1/sqrt(diag(ret)))
  return(ret)
}

clear.dir <- function(dir.path) {
  files.in.dir = list.files(dir.path)
  for (file.in.dir in files.in.dir) {
    full.file.path = file.path(dir.path, file.in.dir)
    file.remove(full.file.path)
  }
}

export.subtype.to.mkl <- function(omics.list, dir.name, run.id=NULL) {
  folder.path = file.path(get.mkl.arguments.path(), dir.name, run.id)
  
  if (!dir.exists(folder.path)) {
    dir.create(folder.path)
  }
  
  kernels.path = file.path(folder.path, 'kernels')
  
  if (!dir.exists(kernels.path)) {
    dir.create(kernels.path)
  }
  clear.dir(kernels.path)
  
  gammas = 10 ** seq(-6, 6, by=3)
  for (i in 1:length(omics.list)) {
    for (j in 1:length(gammas)) {
      datum = omics.list[[i]]
      gamma = gammas[[j]] / (2*nrow(datum)**2)
      mat = radial.basis(datum, gamma)
      R.matlab::writeMat(file.path(kernels.path, paste(i, '_', j, sep='')), mat=mat)
    }
  }
  
  output.path = file.path(folder.path, 'output')
  if (!dir.exists(output.path)) {
    dir.create(output.path)
  }
  clear.dir(output.path)
  
  write.table(colnames(omics.list[[1]]), file=file.path(folder.path, 'ids'),
              quote=F, row.names = F, col.names = F)
  
}

get.mkl.clustering <- function(dir.name, num.clusters=NULL, run.id=NULL) {
  folder.path = file.path(get.mkl.arguments.path(), dir.name, run.id)
  output.path = file.path(folder.path, 'output')
  output.files = list.files(output.path)
  if (is.null(num.clusters)) {
    clustering = read.csv(file.path(output.path, output.files[grep('clusters', output.files)]))[,2]
  } else {
    dr.file.path = file.path(output.path, output.files[grep('dimensions', output.files)])
    dr.rep = read.table(dr.file.path, sep=',', header=T, row.names=1)
    clustering = kmeans(dr.rep, num.clusters, iter.max=100, nstart=30)$cluster
  }
  return(clustering)
}


check.survival <- function(groups, subtype, survival.file.path) {
  if (missing(survival.file.path)) {
    survival.file.path = get.subtype.survival.path(subtype)
  }
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))

  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  return(survdiff(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
  
}

get.subtype.survival.path <- function(subtype) {
  datasets.path = get.dataset.dir.path()
  survival.file.path = file.path(datasets.path, subtype, 'survival')
  return(survival.file.path)
}

get.clustering.silhouette <- function(raw.data, clustering) {
  sils = c()
  for (i in 1:length(raw.data)) {
    x = raw.data[[i]]
    distmatrix = dist2(as.matrix(t(x)),as.matrix(t(x)))
    sil = silhouette(clustering, dmatrix = distmatrix)[,3]
    sils = c(sils, mean(sil))
  }
  return(mean(sils))
}

get.clinical.params.dir <- function() {
  return('CLINICAL PARAMS DIR')
}

get.clinical.params <- function(subtype.name) {
  clinical.data.path = paste(get.clinical.params.dir(), subtype.name, sep = '')
  clinical.params = read.table(clinical.data.path,
                               sep='\t', header=T, row.names = 1, stringsAsFactors = F)
  rownames.with.duplicates = get.fixed.names(rownames(clinical.params))  
  clinical.params = clinical.params[!duplicated(rownames.with.duplicates),]
  rownames(clinical.params) = rownames.with.duplicates[!duplicated(rownames.with.duplicates)]
  return(clinical.params)
}

check.clinical.enrichment <- function(clustering, subtype.name) {
  clinical.params = get.clinical.params(subtype.name)  
  
  clinical.metadata = list(gender='DISCRETE', age_at_initial_pathologic_diagnosis='NUMERIC',
    pathologic_M='DISCRETE', pathologic_N='DISCRETE', pathologic_T='DISCRETE', pathologic_stage='DISCRETE')
  
  pvalues = c()
  
  params.being.tested = c()
  
  for (clinical.param in names(clinical.metadata)) {
    
    if (!(clinical.param %in% colnames(clinical.params))) {
      #print(paste0('WARNING: ', clinical.param, ' does not appear for subtype ', subtype.name))
      next
    }
    
    clinical.values = clinical.params[names(clustering),clinical.param]
    is.discrete.param = clinical.metadata[clinical.param] == 'DISCRETE'
    is.numeric.param = clinical.metadata[clinical.param] == 'NUMERIC'
    stopifnot(is.discrete.param | is.numeric.param)
    
    # skip parameter if many missing values
    
    if (is.numeric.param) {
      numeric.entries = !is.na(as.numeric(clinical.values))
      if (2 * sum(numeric.entries) < length(clinical.values)) {
        #print(paste0('WARNING: skipping on ', clinical.param, ' for subtype ', subtype.name))
        next
      }
    } else {
      not.na.entries = !is.na(clinical.values)
      should.skip = F
      if (2 * sum(not.na.entries) < length(clinical.values)) {
        should.skip = T
      } else if (length(table(clinical.values[not.na.entries])) == 1) {
        should.skip = T
      }
      if (should.skip) {
        #print(paste0('WARNING: skipping on ', clinical.param, ' for subtype ', subtype.name))
        next
      }
    }
    
    params.being.tested = c(params.being.tested, clinical.param)
    
    if (is.discrete.param) {
      #clustering.with.clinical = cbind(clustering, clinical.values)
      #tbl = table(as.data.frame(clustering.with.clinical[!is.na(clinical.values),]))
      #test.res = chisq.test(tbl)
      #pvalue = test.res$p.value
      pvalue = get.empirical.clinical(clustering[!is.na(clinical.values)], clinical.values[!is.na(clinical.values)], T)
      
    } else if (is.numeric.param) {
      #test.res = kruskal.test(as.numeric(clinical.values[numeric.entries]),
      #				clustering[numeric.entries])
      #pvalue = test.res$p.value
      pvalue = get.empirical.clinical(clustering[numeric.entries], as.numeric(clinical.values[numeric.entries]), F)
    }
    
    pvalues = c(pvalues, pvalue)
    
  }
  names(pvalues) = params.being.tested
  return(pvalues)
}

get.empirical.clinical <- function(clustering, clinical.values, is.chisq) {
  set.seed(42)
  if (is.chisq) {
      clustering.with.clinical = cbind(clustering, clinical.values)
      tbl = table(as.data.frame(clustering.with.clinical))
      test.res = chisq.test(tbl)
  } else {
    test.res = kruskal.test(as.numeric(clinical.values), clustering)
  }
  orig.pvalue = test.res$p.value
  num.iter = 1000
  total.num.iters = 0
  total.num.extreme = 0
  should.continue = T
  
  while (should.continue) {
    print('another iteration in empirical clinical')
    perm.pvalues = as.numeric(mclapply(1:num.iter, function(i) {
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
    
      if (is.chisq) {
        clustering.with.clinical = cbind(cur.clustering, clinical.values)
        tbl = table(as.data.frame(clustering.with.clinical))
        test.res = chisq.test(tbl)
      } else {
        test.res = kruskal.test(as.numeric(clinical.values), cur.clustering)
      }
      cur.pvalue = test.res$p.value
      return(cur.pvalue)
    }, mc.cores=60))
    total.num.iters = total.num.iters + num.iter
    total.num.extreme = total.num.extreme + sum(perm.pvalues <= orig.pvalue)
    
    binom.ret = binom.test(total.num.extreme, total.num.iters)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
    
    sig.threshold = 0.05
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if (!is.threshold.in.conf | total.num.iters > 1e5) {
      should.continue = F
    }
  }
  return(cur.pvalue)
}
  

get.nmf.datasets.dir <- function() {
  return('NMF_DATASETS_PATH')
}

save.subtype.matlab.format <- function(subtype.raw.data) {
  data.names = c('1', '2', '3')
  
  full.dir.name = get.nmf.datasets.dir()
  dir.create(full.dir.name, showWarnings = F)
  
  for (i in 1:length(subtype.raw.data)) {
    full.file.name = file.path(full.dir.name, data.names[i])
    write.table(subtype.raw.data[[i]], full.file.name)
  }
}

keep.high.var.features <- function(omic, num.features=2000) {
  if (nrow(omic) < num.features) {
    return(omic)
  } else {
    feature.vars = apply(omic, 1, var)
    threshold = feature.vars[order(feature.vars, decreasing = T)][num.features]
    return(omic[feature.vars >= threshold,])    
  }
}
