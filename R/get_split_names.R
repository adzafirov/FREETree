#' Method for extracting names of splitting features used in a tree.
#'
#' @param x a tree object.
#' @param y train or test set.
#' @return names of splitting features extracted from tree object.
#' @import pre
get_split_names = function(tree,data){

  # path: the string that contains all the node information
  paths <- pre:::list.rules(tree, removecomplements = FALSE)
  vnames = names(data)
  split_names = vnames[sapply(sapply(vnames, FUN = function(var) grep(paste(paste(var,"<="),"|",paste(var,">"),sep=""), paths)), length) > 0]
  return (split_names)
}
