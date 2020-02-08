library("glmertree")
library("WGCNA")
library("pre")

FREEtree = function(data,fixed_regress=NULL,fixed_split=NULL,var_select=NULL,
                    power=6, minModuleSize = 1, cluster,maxdepth_factor_screen=0.04,
                    maxdepth_factor_select=0.5,Fuzzy=TRUE,minsize_multiplier = 5,
                    alpha_screen=0.2, alpha_select=0.2, alpha_predict=0.05){
    ### if there are no features to select, just use fixed_regress and fixed_split      ### no need to filter, just use GLMM tree right away
    if(length(var_select)==0){
        if (length(fixed_regress)==0){
            if (length(fixed_split)==0){
                stop("no features to split and regress on")
                }
            fixed_regress = "1"
            }
        maxdepth = length(fixed_split)
        Formula = as.formula(paste("y~",paste(fixed_regress,collapse = "+"),
                                       "|",cluster,"|", paste(fixed_split,collapse = "+")))
        mytree = lmertree(Formula,data=data,alpha=alpha_predict,maxdepth=maxdepth)
        mytree$final_selection = NULL
        return (mytree)
    } ###
    # Now var_select is not empty
    # If fixed_regress (eg: time) is not specified: use PCs as regressors at screening step
    if (length(fixed_regress)==0){
        cat("Use FREEtree_PC\n")
        return(FREEtree_PC(data=data,fixed_split=fixed_split,
                    var_select=var_select,
                    power=power, minModuleSize = minModuleSize, cluster=cluster,
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
                        power=power, minModuleSize = minModuleSize, cluster=cluster,
                        maxdepth_factor_screen=maxdepth_factor_screen,
                        maxdepth_factor_select=maxdepth_factor_select,Fuzzy=Fuzzy,
                        minsize_multiplier=minsize_multiplier,
                        alpha_screen=alpha_screen, alpha_select=alpha_select,
                        alpha_predict=alpha_predict))

}



