setMethod(f = "combine", 
          signature = c("mgsa", "mgsa"), 
          definition = function(x, y, ...) {
            
            x@moa@call <- y@moa@call <- "merged_moa"
            j1 <- identical(x@moa, y@moa)
            j2 <- length(x@sup@score.pc) == length(x@sup@score.pc)
            j <- j1 & j2
            
            if (!j)
              stop("different moa or different number of PCs retained in x and y")
            
            r <- x
            r@call <- "merged mogsa"
            r@sup@sup <- mapply(cbind, r@sup@sup, y@sup@sup, SIMPLIFY = FALSE)
            r@sup@coord.sep <- mapply(rbind, r@sup@coord.sep, y@sup@coord.sep, SIMPLIFY = FALSE)
            r@sup@coord.comb <- rbind(r@sup@coord.comb, y@sup@coord.comb)
            r@sup@score <- rbind(r@sup@score, y@sup@score)
            r@sup@score.data <- mapply(rbind, r@sup@score.data, y@sup@score.data, SIMPLIFY = FALSE)
            r@sup@score.pc <- mapply(rbind, r@sup@score.pc, y@sup@score.pc, SIMPLIFY = FALSE)
            r@sup@score.sep <- mapply(function(a, b) {
              mapply(c, a, b, SIMPLIFY = FALSE)
            }, x@sup@score.sep, y@sup@score.sep, SIMPLIFY = FALSE)
            r@sup@p.val <- rbind(r@sup@p.val, y@sup@p.val)
            
            return(r)
          })