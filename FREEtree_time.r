
# FREEtree_time: used when var_select (list of features to screen & select) and fixed_regress (eg: time) are non-empty
# FREEtree is equivalent to this FREEtree_time in this case
FREEtree_time = function(data,fixed_regress,fixed_split,var_select, power, minModuleSize, cluster,
                         maxdepth_factor_screen,maxdepth_factor_select,
                         Fuzzy,minsize_multiplier,alpha_screen,alpha_select,
                         alpha_predict){
  # Cluster var_select
  data_WGCNA = data[var_select]
  # Must set numericLabels = FALSE so that it uses actual colors like "grey"
  # "This function performs automatic network construction and module detection 
  #  on large expression datasets in a block-wise manner."
  # !!! attention needs to be paid to minModuleSize
  net = blockwiseModules(data_WGCNA, power = power,TOMType = "unsigned", 
                         minModuleSize = minModuleSize, reassignThreshold = 0, 
                         mergeCutHeight = 0.25,numericLabels = FALSE, 
                         pamRespectsDendro = FALSE,verbose = 0)
  
  
  # the correspondance between feature names and colors
  colors = net$colors # it is a string vector with names (that is the name is V1)
  print(colors)
  module_names = unique(colors) # all color names
  #"dictionary" with keys=name of color,values=names of features of that color
  module_dic = list() 
  for (i in 1:length(module_names)){
    module_dic[[module_names[i]]] = names(colors[colors==module_names[i]])
  }
  
  imp_var = list() # used to store the names of important features
  
  if(Fuzzy==TRUE){
    # Do the selection step just like Fuzzy Forest:
    # for each module (including grey), use them as split_var and select
    # finally, use all selected ones as split_var and select
    
    for (name in module_names){
      # in the formula, add fixed_split as split_var, also include the module features
      split_var = c(module_dic[[name]],fixed_split)
      maxdepth = ceiling(maxdepth_factor_screen*length(split_var))
      # use fixed_regress as regressor
      regress_var = fixed_regress
      
      # Formula for lmtree
      Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                 "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
      
      # fit the tree
      mytree = lmertree(Formula, data = data,alpha=alpha_screen,maxdepth=maxdepth) 
      
      #extract important features
      imp_var[[length(imp_var)+1]] = get_split_names(mytree$tree,data)
    }
    
    # the variables selected from all the modules
    final_var = imp_var[[1]]
    if (length(imp_var)>1){
      # There are at least 2 modules
      for (i in 2:length(imp_var)){
        final_var = c(final_var,imp_var[[i]])
      }
      cat("after screening within modules",final_var,"\n")
      
      # the final selection among all the chosen features 
      regress_var = fixed_regress
      split_var = c(final_var,fixed_split)
      maxdepth = ceiling(maxdepth_factor_select*length(split_var))
      Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                 "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
      mytree = lmertree(Formula, data = data,alpha=alpha_select,maxdepth=maxdepth) 
      final_var = get_split_names(mytree$tree,data)
      cat("final features",final_var)      
    }else{
      # only grey module
      cat("There is only one module, final features",final_var)
    }
    # use the final features as split&regression variables
    split_var = c(final_var,fixed_split)
    maxdepth = length(split_var)
    regress_var = c(fixed_regress,final_var)
    Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                               "|",cluster,"|",
                               paste(split_var,collapse = "+")))
    minsize = round(minsize_multiplier*length(regress_var))
    mytree = lmertree(Formula, data = data,alpha=alpha_predict,maxdepth=maxdepth,
                      minsize = minsize)
    mytree$final_selection = final_var
    return(mytree)           
  }
  if(Fuzzy==FALSE){
    # first do the screening and selecting in non-grey modules
    # Then use those non-grey estimated true features as regressors
    # and grey features as split_var, choose grey features and keep them
    for (name in module_names){
      if(name=="grey"){
        next # next loop item
      }
      split_var = c(module_dic[[name]],fixed_split)
      maxdepth = ceiling(maxdepth_factor_screen*length(split_var))
      regress_var = fixed_regress
      
      # Formula for lmtree
      Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                 "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
      
      # fit the tree
      mytree = lmertree(Formula, data = data,alpha=alpha_screen,maxdepth=maxdepth) 
      
      #extract important features
      imp_var[[length(imp_var)+1]] = get_split_names(mytree$tree,data)
    }# Now imp_var contains important features from modules that are not grey
    if(length(imp_var)==0){
      # only grey module, no other modules
      split_var = c(module_dic[["grey"]],fixed_split)
      maxdepth = ceiling(maxdepth_factor_screen*length(split_var))
      regress_var = fixed_regress
      Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                 "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
      mytree = lmertree(Formula, data = data,alpha = alpha_screen,maxdepth=maxdepth) 
      final_var = get_split_names(mytree$tree,data)
      cat("There is only one module which is grey, final features",final_var)
    }else{
      # at least one non-grey module
      final_var = imp_var[[1]]
      # if only one non-grey module: final_var is the chosen non-grey features
      # if at least two non-grey modules:
      if (length(imp_var)>1){
        # There are at least 2 modules
        for (i in 2:length(imp_var)){
          final_var = c(final_var,imp_var[[i]])
        }
        
        cat("After screening within non-grey modules",final_var,"\n")
        # select from selected non-grey features
        regress_var = fixed_regress
        split_var = c(final_var,fixed_split)
        maxdepth = ceiling(maxdepth_factor_select*length(split_var))
        Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                   "|",cluster,"|",
                                   paste(split_var,collapse = "+")))
        mytree = lmertree(Formula, data = data,alpha = alpha_select,maxdepth=maxdepth) 
        final_var = get_split_names(mytree$tree,data)
      }
      cat("The chosen non-grey features are",final_var,"\n")
      
      # use final_var (chosen non-grey features) to select grey features
      regress_var = c(fixed_regress,final_var)
      split_var = c(module_dic[["grey"]],fixed_split)
      maxdepth = ceiling(maxdepth_factor_screen*length(split_var))
      Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                 "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
      mytree = lmertree(Formula, data = data,alpha = alpha_screen,maxdepth=maxdepth) 
      grey_var = get_split_names(mytree$tree,data)
      cat("The chosen grey features are",grey_var,"\n")
      # use final_var and grey_var do to the final model tree
      final_var = c(final_var,grey_var)    
      cat("final features",final_var)  
    }
    # use the final features as split&regression variables
    split_var = c(final_var,fixed_split)
    maxdepth = length(split_var)
    regress_var = c(fixed_regress,final_var)
    Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                               "|",cluster,"|",
                               paste(split_var,collapse = "+")))
    minsize = round(minsize_multiplier*length(regress_var))
    mytree = lmertree(Formula, data = data,alpha=alpha_predict,maxdepth=maxdepth,
                      minsize=minsize) 
    mytree$final_selection = final_var
    return(mytree)
  }
  
}