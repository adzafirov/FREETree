#' Version of FREEtree called when fixed_regress is NULL, uses principal
#' components (PC) as regressors for non-grey modules.
#'
#' @param data data to train or test FREEtree on.
#' @param fixed_split user specified char vector of features to be used in
#'   splitting with certainty.
#' @param var_select a char vector containing features to be selected. These
#'   features will be clustered by WGCNA and the chosen ones will be used in
#'   regression and splitting.
#' @param power soft thresholding power parameter of WGCNA.
#' @param minModuleSize minimum possible module size parameter of WGCNA.
#' @param cluster: the variable name of each cluster (in terms of random effect)
#'   using glmer's implementation.
#' @param Fuzzy boolean to indicate desire to screen like \href{https://cran.r-project.org/web/packages/fuzzyforest/fuzzyforest.pdf}{Fuzzy Forest} if Fuzzy
#'   = TRUE; if Fuzzy= FALSE, first screen within non-grey modules and then
#'   select the final non-grey features within the selected ones from each
#'   non-grey module; Use this final non-grey features as regressors (plus
#'   fixed_regress) and use grey features as split_var to select grey features.
#'   Then use final non-grey features and selected grey features together in
#'   splitting and regression variables, to do the final prediction. Fuzzy=FALSE
#'   is used if there are so many non-grey features and you want to protect grey
#'   features.
#' @param maxdepth_factor_screen when selecting features from one module, the
#'   maxdepth of the glmertree is set to ceiling function of
#'   maxdepth_factor_screen*(features in that module). Default is 0.04.
#' @param maxdepth_factor_select Given screened features (from each modules, if
#'   Fuzzy=FALSE, that is the selected non-grey features from each non-grey
#'   modules), we want to select again from those screened features. The
#'   maxdepth of that glmertree is set to be ceiling of
#'   maxdepth_factor_select*(#screened features). Default is 0.6. for the
#'   maxdepth of the prediction tree (final tree), maxdepth is set to the length
#'   of the split_var (fixed+chosen ones).
#' @param minsize_multiplier At the final prediction tree, the minsize =
#'   minsize_multiplier times the length of final regressors. The default is 5.
#'   Note that we only set minsize for the final prediction tree instead of
#'   trees at the feature selection step since during feature selection, we
#'   don't have to be so careful. Note that when tuning the parameters, larger
#'   alpha and samller minsize_multiplier will result in deeper tree and
#'   therefore may cause overfitting problem. It is recommended to decrease
#'   alpha and decrease minsize_multiplier at the same time.
#' @param alpha_screen alpha used in screening step.
#' @param alpha_select alpha used in selection step.
#' @param alpha_predict alpha used in prediction step.
#' @return a glmertree object (trained tree).
#' @import(WGCNA)
FREEtree_PC = function(data, fixed_split, var_select, power, minModuleSize,
                       cluster, maxdepth_factor_screen, maxdepth_factor_select, Fuzzy, minsize_multiplier,
                       alpha_screen, alpha_select, alpha_predict) {
  # Cluster var_select
  data_WGCNA = data[var_select]
  # Must set numericLabels = FALSE so that it uses actual colors like
  # 'grey'
  net = WGCNA::blockwiseModules(data_WGCNA, power = power, TOMType = "unsigned",
                         minModuleSize = minModuleSize, reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = FALSE, pamRespectsDendro = FALSE, verbose = 0)
  # the correspondance betweeen feature names and colors
  colors = net$colors  # it is a string vector with names (that is the name is V1)
  module_names = unique(colors)  # all color names
  #'dictionary'with keys=name of color,values=names of features of that color
  module_dic = list()
  for (i in 1:length(module_names)) {
    module_dic[[module_names[i]]] = names(colors[colors == module_names[i]])
  }

  # extract eigengenes and rename the column The eigengene(1st pricinpal
  # component) is L2 normalized
  eigengene = net$MEs
  # eigengene add eigengen to training data (for grey group, eigengene is
  # meaningless)
  for (name in module_names) {
    if (name == "grey") {
      next
    }
    eigen_name = paste("ME", name, sep = "")
    data[[eigen_name]] = eigengene[[eigen_name]]
  }
  imp_var = list()  # used to store the names of important features

  if (Fuzzy == TRUE) {
    # first screen then select, just like Fuzzy Forest For each module that
    # is not grey, use model tree as following: use its eigengene as
    # regression variables and all features as splitting ones For grey
    # module, use regular REEM (set regressor = '1') Then select by using
    # all chosen features as split_var and regress='1' Finally, use all the
    # selected features for splitting and regression variables

    for (name in module_names) {
      split_var = c(module_dic[[name]], fixed_split)
      maxdepth = ceiling(maxdepth_factor_screen * length(split_var))
      # use eigengene as regressor
      if (name == "grey") {
        regress_var = "1"
      } else {
        eigen_name = paste("ME", name, sep = "")
        regress_var = eigen_name
      }
      # Formula for lmtree: use PC as regressors
      Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                 "|", cluster, "|", paste(split_var, collapse = "+")))

      # fit the tree
      mytree = lmertree(Formula, data = data, alpha = alpha_screen,
                        maxdepth = maxdepth)
      # extract important features
      imp_var[[length(imp_var) + 1]] = get_split_names(mytree$tree,
                                                       data)
    }
    # the variables selected from all the modules
    final_var = imp_var[[1]]
    if (length(imp_var) > 1) {
      for (i in 2:length(imp_var)) {
        final_var = c(final_var, imp_var[[i]])
      }
      cat("After screening within modules ", final_var, "\n")
      # select features again use all selected features as split_var with no
      # regressors
      split_var = c(final_var, fixed_split)
      maxdepth = ceiling(maxdepth_factor_select * length(split_var))
      Formula = as.formula(paste("y~", "1", "|", cluster, "|", paste(split_var,
                                                                     collapse = "+")))
      mytree = lmertree(Formula, data = data, alpha = alpha_select,
                        maxdepth = maxdepth)
      final_var = get_split_names(mytree$tree, data)
      cat("Final features ", final_var, "\n")

    } else {
      # length(imp_var) now is 1, only one module
      cat("There is only one module ", final_var, "\n")
    }

    # use the final features as split&regression variables
    split_var = c(final_var, fixed_split)
    maxdepth = length(split_var)
    Formula = as.formula(paste("y~", paste(final_var, collapse = "+"),
                               "|", cluster, "|", paste(split_var, collapse = "+")))
    minsize = round(minsize_multiplier * length(final_var))
    mytree = lmertree(Formula, data = data, alpha = alpha_predict,
                      maxdepth = maxdepth, minsize = minsize)
    mytree$final_selection = final_var
    return(mytree)
  }

  if (Fuzzy == FALSE) {
    # select features from non-grey modules and use them as regressors to
    # select features from grey module. Then use all of them as split and
    # regressor

    # for non-grey groups
    for (name in module_names) {
      if (name == "grey") {
        next
      }
      split_var = c(module_dic[[name]], fixed_split)
      maxdepth = ceiling(maxdepth_factor_screen * length(split_var))
      # use eigengene as regressor
      eigen_name = paste("ME", name, sep = "")
      regress_var = eigen_name
      # Formula for lmtree: use PC as regressors
      Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                 "|", cluster, "|", paste(split_var, collapse = "+")))

      # fit the tree
      mytree = lmertree(Formula, data = data, alpha = alpha_screen,
                        maxdepth = maxdepth)
      # extract important features
      imp_var[[length(imp_var) + 1]] = get_split_names(mytree$tree,
                                                       data)
    }
    # Now imp_var contains all the non-grey screened features
    if (length(imp_var) == 0) {
      # there is only one module which is grey just select from the grey
      # module
      split_var = c(module_dic[["grey"]], fixed_split)
      maxdepth = ceiling(maxdepth_factor_screen * length(split_var))
      regress_var = "1"
      Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                 "|", cluster, "|", paste(split_var, collapse = "+")))
      mytree = lmertree(Formula, data = data, alpha = alpha_screen,
                        maxdepth = maxdepth)
      final_var = get_split_names(mytree$tree, data)
      cat("There is only grey module ", final_var, "\n")
    }
    if (length(imp_var) == 1) {
      # only one non-grey module
      final_var = imp_var[[1]]
      cat("There is only one non-grey module", final_var, "\n")
      # use final_var as regressors
      split_var = c(module_dic[["grey"]], fixed_split)
      maxdepth = ceiling(maxdepth_factor_screen * length(split_var))
      regress_var = final_var
      Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                 "|", cluster, "|", paste(split_var, collapse = "+")))
      mytree = lmertree(Formula, data = data, alpha = alpha_screen,
                        maxdepth = maxdepth)
      grey_var = get_split_names(mytree$tree, data)
      final_var = c(final_var, grey_var)
      cat("The final features ", final_var, "\n")
    }
    if (length(imp_var) > 1) {
      # at least two non-grey modules, select among non-grey modules
      final_var = imp_var[[1]]
      for (i in 2:length(imp_var)) {
        final_var = c(final_var, imp_var[[i]])
      }
      cat("After screening from non-grey modules ", final_var, "\n")
      # now final_var contains all the non-grey screened features
      split_var = c(final_var, fixed_split)
      maxdepth = ceiling(maxdepth_factor_select * length(split_var))
      regress_var = "1"
      Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                 "|", cluster, "|", paste(split_var, collapse = "+")))
      mytree = lmertree(Formula, data = data, alpha = alpha_select,
                        maxdepth = maxdepth)
      final_var = get_split_names(mytree$tree, data)
      # Now final_var contains final non-grey features
      cat("Final non-grey features ", final_var, "\n")
      # use final_var as regressors and select features from grey features
      split_var = c(module_dic[["grey"]], fixed_split)
      maxdepth = ceiling(maxdepth_factor_screen * length(split_var))
      regress_var = final_var
      Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                 "|", cluster, "|", paste(split_var, collapse = "+")))
      mytree = lmertree(Formula, data = data, alpha = alpha_screen,
                        maxdepth = maxdepth)
      grey_var = get_split_names(mytree$tree, data)
      cat("Final grey features ", grey_var, "\n")
      final_var = c(final_var, grey_var)
      cat("The final features ", final_var, "\n")
    }
    # use the final features as split&regression variables
    split_var = c(final_var, fixed_split)
    maxdepth = length(split_var)
    Formula = as.formula(paste("y~", paste(final_var, collapse = "+"),
                               "|", cluster, "|", paste(split_var, collapse = "+")))
    minsize = round(minsize_multiplier * length(final_var))
    mytree = lmertree(Formula, data = data, alpha = alpha_predict,
                      maxdepth = maxdepth, minsize = minsize)
    mytree$final_selection = final_var
    return(mytree)
  }


}
