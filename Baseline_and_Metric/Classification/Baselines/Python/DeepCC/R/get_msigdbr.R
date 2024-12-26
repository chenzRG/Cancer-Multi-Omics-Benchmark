#' get MSigDBr from R package `msigdbr`
#'
#' This function defines MSigDBr from a R package, with 25, 724 gene sets
#'
#' @param cores a integer indicating cpu cores used in parallel computing (default = all cores -2 )
#'
#' @return a list containing 25, 724 gene sets, each sets contains multiple entrez_gene
#' @examples
#' MSigDBr <- get_msigdbr()

get_msigdbr <- function(cores){
  m_df <- msigdbr(species = "Homo sapiens")
  set_name <- unique(m_df$gs_name)

  m_df_2 <- m_df %>%
    select(gs_name, entrez_gene)

  get_list <- function(g_name){
    tmp <- m_df_2 %>%
      filter(gs_name == g_name)  %>%
      select(entrez_gene) %>%
      sapply(as.character) %>%
      as.character()
    tmp
  }

  doParallel::registerDoParallel(cores)

  gene_sets <- foreach(idx = 1:length(set_name)) %dopar% {
    get_list(set_name[idx])
  }
  names(gene_sets) <- set_name

  gene_sets
}
