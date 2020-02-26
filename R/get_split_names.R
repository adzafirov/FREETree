#' Method for extracting names of splitting features used in a tree.
#'
#' @param tree a tree object.
#' @param data train or test set.
#' @import pre
#' @return names of splitting features extracted from tree object.
get_split_names = function(tree,data){

  # path: the string that contains all the node information
  paths <- eval(parse(text = "pre:::list.rules(tree, removecomplements = FALSE)"))
  vnames = names(data)
  split_names = vnames[sapply(sapply(vnames, FUN = function(var) grep(paste(paste(var,"<="),"|",paste(var,">"),sep=""), paths)), length) > 0]
  return (split_names)
}
