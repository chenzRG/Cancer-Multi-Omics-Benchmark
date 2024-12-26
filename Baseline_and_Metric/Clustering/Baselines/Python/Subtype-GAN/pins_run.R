require(optparse)
require(PINSPlus)

option_list = list(
  make_option(c("-i", "--input"), action = "store", default = '',
              type = 'character', help = "all_data"),
  make_option(c("-t", "--type"), action = "store", default = 'BRCA',
              type = 'character', help = "all_data"),
  make_option(c("-n", "--cluster_num"), action = "store", default = 5,
              type = 'integer', help = "cluster_num")
)
opt = parse_args(OptionParser(option_list = option_list))
flist <- scan(opt$input, what = "", sep = "\n")
cluster_num = opt$cluster_num
nb_file <- length(flist)
l <- list()
for (i in 1:nb_file) {
  print(flist[i])
  if (i == 1) {
    ids <- colnames(read.table(flist[i], check.names = FALSE, row.names = 1, sep = ",", header = TRUE))
    print(ids)
  }
  l[[i]] <- t(data.matrix(read.table(flist[i], check.names = FALSE, row.names = 1, sep = ",", header = TRUE)))
}
print(l[[1]][1:10, 1:10])
start <- Sys.time()
pins.ret = PINSPlus::SubtypingOmicsData(dataList = l, k = opt$cluster_num)
end <- as.numeric(Sys.time() - start, units = 'secs')
clustering = pins.ret$cluster2
print(clustering)

df <- as.data.frame(clustering, col.names = c("pins"))
fileout <- paste("./results/", opt$type, ".pins", sep = "")
write.table(df, file = fileout, quote = FALSE, sep = "\t")
df <- as.data.frame(c(end))
colnames(df) <- c('time')
fileout <- paste("./results/", opt$type, ".pins.time", sep = "")
write.table(df, file = fileout, quote = FALSE, sep = "\t")

