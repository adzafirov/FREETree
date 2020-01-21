library("glmertree")
library("WGCNA")
library("pre")

FREEtree = function(data,fixed_regress=NULL,fixed_split=NULL,var_select=NULL,
                    power=6,cluster,maxdepth_factor_screen=0.04,
                    maxdepth_factor_select=0.5,Fuzzy=TRUE,minsize_multiplier = 5,
                    alpha_screen=0.2, alpha_select=0.2, alpha_predict=0.05){
    ### if there are no features to select, just use fixed_regress and fixed_split
    if(length(var_select)==0){
        if (length(fixed_regress)==0){
            if (length(fixed_split)==0){
                stop("no features to split and regress on")
            }
            fixed_regress = "1"
        }
        maxdepth = length(fixed_split)
        Formula = as.formula(paste("y~",paste(fixed_regress,collapse = "+"),
                                       "|",cluster,"|",
                                     paste(fixed_split,collapse = "+")))
        mytree = lmertree(Formula,data=data,alpha=alpha_predict,maxdepth=maxdepth)
        mytree$final_selection = NULL
        return (mytree)
    } ###
    # Now var_select is not empty
    # If don't specify fixed_regress: use PC as regressors at screening step
    if (length(fixed_regress)==0){
        cat("Use FREEtree_PC\n")
        return(FREEtree_PC(data=data,fixed_split=fixed_split,
                    var_select=var_select,
                    power=power,cluster=cluster,
                    maxdepth_factor_screen=maxdepth_factor_screen, 
                    maxdepth_factor_select=maxdepth_factor_select,Fuzzy=Fuzzy,
                    minsize_multiplier=minsize_multiplier,
                    alpha_screen=alpha_screen, alpha_select=alpha_select, 
                    alpha_predict=alpha_predict))
    }###
    # Now we have non-empty var_select,fixed_regress
    cat("Use FREEtree_time\n")
    return(FREEtree_time(data=data,fixed_regress=fixed_regress,
                        fixed_split=fixed_split, var_select=var_select,
                    power=power,cluster=cluster,
                    maxdepth_factor_screen=maxdepth_factor_screen, 
                    maxdepth_factor_select=maxdepth_factor_select,Fuzzy=Fuzzy,
                    minsize_multiplier=minsize_multiplier,
                    alpha_screen=alpha_screen, alpha_select=alpha_select, 
                    alpha_predict=alpha_predict))
                        
}

# FREEtree_time: used when var_select and fixed_regress are non-empty
# FREEtree is equivalent to this FREEtree_time in this case
FREEtree_time = function(data,fixed_regress,fixed_split,var_select,power,cluster,
                         maxdepth_factor_screen,maxdepth_factor_select,
                         Fuzzy,minsize_multiplier,alpha_screen,alpha_select,
                         alpha_predict){
    # Cluster var_select
    data_WGCNA = data[var_select]
    # Must set numericLabels = FALSE so that it uses actual colors like "grey"
    net = blockwiseModules(data_WGCNA, power = power,TOMType = "unsigned", 
                           minModuleSize = 30,reassignThreshold = 0, 
                           mergeCutHeight = 0.25,numericLabels = FALSE, 
                           pamRespectsDendro = FALSE,verbose = 0)
    # the correspondance betweeen feature names and colors
    colors = net$colors # it is a string vector with names (that is the name is V1)
    module_names = unique(colors) # all color names
    #"dictionary"with keys=name of color,values=names of features of that color
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
            next
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

# FREEtree_PC: used when fixed_regress are NULL, use PC as regressors for non-grey module
# FREEtree is equivalent to this FREEtree_PC in this case
FREEtree_PC = function(data,fixed_split, var_select, power,cluster,
                    maxdepth_factor_screen,maxdepth_factor_select,Fuzzy,
                    minsize_multiplier,alpha_screen,alpha_select,
                    alpha_predict){
    # Cluster var_select
    data_WGCNA = data[var_select]
    # Must set numericLabels = FALSE so that it uses actual colors like "grey"
    net = blockwiseModules(data_WGCNA, power = power,TOMType = "unsigned", 
                           minModuleSize = 30,reassignThreshold = 0, 
                           mergeCutHeight = 0.25,numericLabels = FALSE, 
                           pamRespectsDendro = FALSE,verbose = 0)
    # the correspondance betweeen feature names and colors
    colors = net$colors # it is a string vector with names (that is the name is V1)
    module_names = unique(colors) # all color names
    #"dictionary"with keys=name of color,values=names of features of that color
    module_dic = list() 
    for (i in 1:length(module_names)){
        module_dic[[module_names[i]]] = names(colors[colors==module_names[i]])
    }
    
    # extract eigengenes and rename the column
    # The eigengene(1st pricinpal component) is L2 normalized
    eigengene = net$MEs
    # eigengene
    # add eigengen to training data (for grey group, eigengene is meaningless)
    for (name in module_names){
        if (name == "grey"){
            next
        }
        eigen_name = paste("ME",name,sep="")
        data[[eigen_name]] = eigengene[[eigen_name]]
    }
    imp_var = list() # used to store the names of important features
    
    if (Fuzzy==TRUE){
        # first screen then select, just like Fuzzy Forest
        # For each module that is not grey, use model tree as following:
        # use its eigengene as regression variables and all features as splitting ones
        # For grey module, use regular REEM (set regressor = "1")
        # Then select by using all chosen features as split_var and regress="1"
        # Finally, use all the selected features for splitting and regression variables
       
        for (name in module_names){
            split_var = c(module_dic[[name]],fixed_split)
            maxdepth = ceiling(maxdepth_factor_screen*length(split_var))
            # use eigengene as regressor
            if (name == "grey"){
                regress_var = "1"
            }else{
                eigen_name = paste("ME",name,sep="")
                regress_var = eigen_name
            }
            # Formula for lmtree: use PC as regressors
            Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                   "|",cluster,"|",
                                 paste(split_var,collapse = "+")))

            # fit the tree
            mytree = lmertree(Formula, data=data,alpha=alpha_screen,maxdepth=maxdepth) 
            #extract important features
            imp_var[[length(imp_var)+1]] = get_split_names(mytree$tree,data)
        }        
        # the variables selected from all the modules
        final_var = imp_var[[1]]      
        if (length(imp_var)>1){
            for (i in 2:length(imp_var)){
            final_var = c(final_var,imp_var[[i]])
         }
            cat("After screening within modules ",final_var,"\n")
            # select features again
            # use all selected features as split_var with no regressors
            split_var = c(final_var,fixed_split)
            maxdepth = ceiling(maxdepth_factor_select*length(split_var))
            Formula = as.formula(paste("y~","1",
                                       "|",cluster,"|",
                                     paste(split_var,collapse = "+")))
            mytree = lmertree(Formula, data = data,alpha=alpha_select,maxdepth=maxdepth)
            final_var = get_split_names(mytree$tree,data)
            cat("Final features ",final_var,"\n")
            
        }else{
            # length(imp_var) now is 1, only one module
            cat("There is only one module ",final_var,"\n")
        }
        
        # use the final features as split&regression variables
        split_var = c(final_var,fixed_split)
        maxdepth = length(split_var)
        Formula = as.formula(paste("y~",paste(final_var,collapse = "+"),
                                   "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
        minsize = round(minsize_multiplier*length(final_var))
        mytree = lmertree(Formula, data = data,alpha=alpha_predict,maxdepth=maxdepth,
                         minsize = minsize)
        mytree$final_selection = final_var
        return (mytree)
    }
    
    if (Fuzzy== FALSE){
        # select features from non-grey modules and use them as regressors 
        # to select features from grey module. Then use all of them as split and regressor
        
        # for non-grey groups
        for (name in module_names){
            if (name == "grey"){
                next
            }
            split_var = c(module_dic[[name]],fixed_split)
            maxdepth = ceiling(maxdepth_factor_screen*length(split_var))
            # use eigengene as regressor
            eigen_name = paste("ME",name,sep="")
            regress_var = eigen_name
            # Formula for lmtree: use PC as regressors
            Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                   "|",cluster,"|",
                                 paste(split_var,collapse = "+")))

            # fit the tree
            mytree = lmertree(Formula, data=data,alpha=alpha_screen,maxdepth=maxdepth) 
            #extract important features
            imp_var[[length(imp_var)+1]] = get_split_names(mytree$tree,data)
        }
        # Now imp_var contains all the non-grey screened features
        if(length(imp_var)==0){
            # there is only one module which is grey
            # just select from the grey module
            split_var = c(module_dic[["grey"]],fixed_split)
            maxdepth = ceiling(maxdepth_factor_screen*length(split_var))
            regress_var = "1"
            Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                   "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
            mytree = lmertree(Formula, data=data,alpha=alpha_screen,maxdepth=maxdepth) 
            final_var = get_split_names(mytree$tree,data)
            cat("There is only grey module ",final_var,"\n")
        }
        if(length(imp_var)==1){
            # only one non-grey module
            final_var = imp_var[[1]]
            cat("There is only one non-grey module",final_var,"\n")
            # use final_var as regressors
            split_var = c(module_dic[["grey"]],fixed_split)
            maxdepth = ceiling(maxdepth_factor_screen*length(split_var))
            regress_var = final_var
            Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                   "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
            mytree = lmertree(Formula, data=data,alpha=alpha_screen,maxdepth=maxdepth) 
            grey_var = get_split_names(mytree$tree,data)
            final_var = c(final_var,grey_var)
            cat("The final features ",final_var,"\n")
        }
        if(length(imp_var)>1){
            # at least two non-grey modules, select among non-grey modules
            final_var = imp_var[[1]]
            for (i in 2:length(imp_var)){
                final_var = c(final_var,imp_var[[i]])
                }
            cat("After screening from non-grey modules ",final_var,"\n")
            # now final_var contains all the non-grey screened features
            split_var = c(final_var,fixed_split)  
            maxdepth = ceiling(maxdepth_factor_select*length(split_var))
            regress_var = "1"
            Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                   "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
            mytree = lmertree(Formula, data=data,alpha=alpha_select,maxdepth=maxdepth)
            final_var = get_split_names(mytree$tree,data)
            # Now final_var contains final non-grey features
            cat("Final non-grey features ",final_var,"\n")
            # use final_var as regressors and select features from grey features
            split_var = c(module_dic[["grey"]],fixed_split)
            maxdepth = ceiling(maxdepth_factor_screen*length(split_var))
            regress_var = final_var
            Formula = as.formula(paste("y~",paste(regress_var,collapse = "+"),
                                   "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
            mytree = lmertree(Formula, data=data,alpha=alpha_screen,maxdepth=maxdepth) 
            grey_var = get_split_names(mytree$tree,data)
            cat("Final grey features ",grey_var,"\n")
            final_var = c(final_var,grey_var)
            cat("The final features ",final_var,"\n")
        }
        # use the final features as split&regression variables
        split_var = c(final_var,fixed_split)
        maxdepth = length(split_var)
        Formula = as.formula(paste("y~",paste(final_var,collapse = "+"),
                                   "|",cluster,"|",
                                 paste(split_var,collapse = "+")))
        minsize = round(minsize_multiplier*length(final_var))
        mytree = lmertree(Formula, data = data,alpha=alpha_predict,maxdepth=maxdepth,
                         minsize = minsize)
        mytree$final_selection = final_var
        return (mytree)
    }
    
    
}