# Methods for extracting names of splitting features used in a tree
# tree: a tree object; data: the train or test set
get_split_names = function(tree,data){
  # path: the string that contains all the node information
  paths <- pre:::list.rules(tree, removecomplements = FALSE)
  vnames = names(data)
  # the regex for a variable
  # tomatch = paste(paste(var,"<="),"|",paste(var,">"),sep="")
  # match to tomatch in path
  tmp = vnames[sapply(sapply(vnames, FUN = function(var) grep(paste(paste(var,"<="),"|",paste(var,">"),sep=""), paths)), length) > 0]
  return (tmp)
}