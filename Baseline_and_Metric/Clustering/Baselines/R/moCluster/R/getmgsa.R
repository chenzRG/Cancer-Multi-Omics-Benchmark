getmgsa <- function(mgsa, value) {
  
  if (!inherits(mgsa, "mgsa"))
    stop("the first input argument should be an object of class 'mgsa'.")
  if (value %in% c("call", "moa", "sup"))
    r <- slot(mgsa, value) else 
      if (value %in% c("eig", "tau", "partial.eig", "eig.vec", "loading", 
                       "fac.scr", "partial.fs", "ctr.obs", "ctr.var", "ctr.tab", "RV"))
        r <- slot(mgsa@moa, value) else 
          if (value %in% c("data", "coord.sep", "coord.comb", "score", 
                           "score.data", "score.pc", "score.sep", "p.val")) 
            r <- slot(mgsa@sup, value) else
              stop("unknown value selected.")
  return(r)
}
