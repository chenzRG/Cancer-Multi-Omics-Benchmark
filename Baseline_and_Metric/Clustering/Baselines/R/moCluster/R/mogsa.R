mogsa <- function(x, sup, nf=NULL, factors = NULL, proc.row=NULL, w.data=NULL, 
  w.row=NULL, statis=FALSE, ks.stat=FALSE, ks.B = 1000, ks.cores = NULL, p.adjust.method = "fdr") {
  
  # if sup is NULL .....
  # extract data and moa
  if (inherits(x, "list")) {
    if (is.null(nf) & is.null(factors))
      stop("x is an object of \"list\", nf or factors need to be set.")
    if (is.null(proc.row))
      stop("x is an object of \"list\", proc.row need to be set.")
    if (is.null(w.data))
      stop("x is an object of \"list\", w.data need to be set.")
    if (is.null(statis))
      stop("x is an object of \"list\", statis need to be set.")
    data <- x
    r <- moa(data=data, proc.row=proc.row, 
      w.data=w.data, w.row=w.row, statis=statis)
    
    # sup data
    if (inherits(sup, "list")) {
      supr <- sup.moa (X=r, sup=sup, nf=nf, factors = factors, 
        ks.stat=ks.stat, ks.B = ks.B, ks.cores = ks.cores, p.adjust.method = p.adjust.method)
    } else if (inherits(sup, "moa.sup")) {
      stop("sup cannot be an object of class moa.sup if x is an object of class list.")
    } 
  } else if (inherits(x, "moa")) {
    if (!is.null(proc.row))
      cat("x is an object of \"moa\", proc.row is not used")
    if (!is.null(w.data))
      cat("x is an object of \"moa\", w.data is not used")
    if (!is.null(w.row))
      cat("x is an object of \"moa\", w.row is not used")
    if (!is.null(statis))
      cat("x is an object of \"moa\", statis is not used")
    data <- x@data
    r <- x

    # sup data
    if (inherits(sup, "list")) {
      supr <- sup.moa (X=r, sup=sup, nf=nf, factors = factors, 
        ks.stat=ks.stat, ks.B = ks.B, ks.cores = ks.cores, p.adjust.method = p.adjust.method)
    } else if (inherits(sup, "moa.sup")) {
      if (!is.null(nf))
        cat("x is an object of \"moa\" and sup is an object of class \"sup.moa\", nf is not used")
      supr <- sup
    }
  }
  
  new("mgsa",
      call=match.call(),
      moa=r,
      sup=supr)
}