setMethod("summary", signature = "moa", 
  definition = function(object) {
  nc <- ncol(object@fac.scr)
  list("@call" = object@call,
       "@tab.dim" = object@tab.dim,
       "@ids"  = lapply(object@data, function(x) head(rownames(x))),
       "@RV" = object@RV,
       "@eig" = head(object@eig),
       "@fac.scr" = summary(object@fac.scr[, 1:min(6, nc)]),
       "@loading" = summary(object@loading[, 1:min(6, nc)]))
})

setMethod("print", signature = "moa", 
  definition = function(x) {
  sx <- summary(x)
  print(sx)
  invisible(sx)
})



setMethod("summary", signature = "moa.sup", 
  definition=  function(object) {
  list(
    "@sup.dim" = sapply(object@sup, dim),
    "@score" = summary(object@score)[, 1:min(6, ncol(object@score))], 
    "@gsets" = head(rownames(object@score))
  )
})

setMethod("print", signature = "moa.sup",
  definition = function(x) {
  sx <- summary(x)
  print(sx)
  invisible(sx)
})


setMethod("summary", signature = "mgsa", 
  definition = function(object) {
  list("@call" = object@call, 
       "@moa" = summary(object@moa), 
       "@sup"  = summary(object@sup)
  )
})

setMethod("print", signature = "mgsa", 
  definition = function(x) {
  sx <- summary(x)
  print(sx)
  invisible(sx)
})


setMethod(f = "show", 
          signature = "mgsa", 
          definition = function(object) {
            print(object)
          })

setMethod(f = "show", 
          signature = "moa", 
          definition = function(object) {
            print(object)
          })

setMethod(f = "show", 
          signature = "moa.sup", 
          definition = function(object) {
            print(object)
          })