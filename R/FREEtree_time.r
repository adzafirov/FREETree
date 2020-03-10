#' Version of FREEtree called when var_select and fixed_regress are specified,
#'
#' @param data data to train or test FREEtree on.
#' @param fixed_regress user specified char vector of regressors that will never
#'   be screened out; if fixed_regress = NULL, method uses PC as regressor at
#'   screening step.
#' @param consider_split user specified char vector of features that will be considered for
#'   splitting in the final tree and will not be filtered out in the screening step.
#' @param var_select a char vector containing features to be selected. These
#'   features will be clustered by WGCNA and the chosen ones will be used in
#'   regression and splitting.
#' @param power soft thresholding power parameter of WGCNA.
#' @param minModuleSize minimum possible module size parameter of WGCNA.
#' @param cluster the variable name of each cluster (in terms of random effect)
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
#' @import WGCNA
#' @import glmertree
#' @importFrom stats as.formula
#' @return a glmertree object (trained tree).

FREEtree_time <- function(data, fixed_regress, consider_split, var_select,
                         power, minModuleSize, cluster, maxdepth_factor_screen, maxdepth_factor_select,
                         Fuzzy, minsize_multiplier, alpha_screen, alpha_select, alpha_predict) {
  # Cluster var_select
  data_WGCNA = data[var_select]
  # Must set numericLabels = FALSE so that it uses actual colors like
  # 'grey' 'This function performs automatic network construction and
  # module detection on large expression datasets in a block-wise
  # manner.'  !!! attention needs to be paid to minModuleSize
  net = blockwiseModules(data_WGCNA, power = power, TOMType = "unsigned",
                         minModuleSize = minModuleSize, reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = FALSE, pamRespectsDendro = FALSE, verbose = 0)


  # the correspondance between feature names and colors
  colors = net$colors  # it is a string vector with names (that is the name is V1)
  print(colors)
  module_names = unique(colors)  # all color names
  #'dictionary' with keys=name of color,values=names of features of that color
  module_dic = list()
  for (i in 1:length(module_names)) {
    module_dic[[module_names[i]]] = names(colors[colors == module_names[i]])
  }

  imp_var = list()  # used to store the names of important features

  if (Fuzzy == TRUE) {
    # Do the selection step just like Fuzzy Forest: for each module
    # (including grey), use them as split_var and select finally, use all
    # selected ones as split_var and select

    for (name in module_names) {
      # in the formula, add consider_split as split_var, also include the module
      # features
      split_var = c(module_dic[[name]], consider_split)
      maxdepth = ceiling(maxdepth_factor_screen * length(split_var))
      # use fixed_regress as regressor
      regress_var = fixed_regress

      # Formula for lmtree
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
      # There are at least 2 modules
      for (i in 2:length(imp_var)) {
        final_var = c(final_var, imp_var[[i]])
      }
      cat("after screening within modules", final_var, "\n")

      # the final selection among all the chosen features
      regress_var = fixed_regress
      split_var = c(final_var, consider_split)
      maxdepth = ceiling(maxdepth_factor_select * length(split_var))
      Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                 "|", cluster, "|", paste(split_var, collapse = "+")))
      mytree = lmertree(Formula, data = data, alpha = alpha_select,
                        maxdepth = maxdepth)
      final_var = get_split_names(mytree$tree, data)
      cat("final features", final_var)
    } else {
      # only grey module
      cat("There is only one module, final features", final_var)
    }
    # use the final features as split&regression variables
    split_var = c(final_var, consider_split)
    maxdepth = length(split_var)
    regress_var = c(fixed_regress, final_var)
    Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                               "|", cluster, "|", paste(split_var, collapse = "+")))
    minsize = round(minsize_multiplier * length(regress_var))
    mytree = lmertree(Formula, data = data, alpha = alpha_predict,
                      maxdepth = maxdepth, minsize = minsize)
    mytree$final_selection = final_var
    return(mytree)
  }
  if (Fuzzy == FALSE) {
    # first do the screening and selecting in non-grey modules Then use
    # those non-grey estimated true features as regressors and grey
    # features as split_var, choose grey features and keep them
    for (name in module_names) {
      if (name == "grey") {
        next  # next loop item
      }
      split_var = c(module_dic[[name]], consider_split)
      maxdepth = ceiling(maxdepth_factor_screen * length(split_var))
      regress_var = fixed_regress

      # Formula for lmtree
      Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                 "|", cluster, "|", paste(split_var, collapse = "+")))

      # fit the tree
      mytree = lmertree(Formula, data = data, alpha = alpha_screen,
                        maxdepth = maxdepth)

      # extract important features
      imp_var[[length(imp_var) + 1]] = get_split_names(mytree$tree,
                                                       data)
    }  # Now imp_var contains important features from modules that are not grey
    if (length(imp_var) == 0) {
      # only grey module, no other modules
      split_var = c(module_dic[["grey"]], consider_split)
      maxdepth = ceiling(maxdepth_factor_screen * length(split_var))
      regress_var = fixed_regress
      Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                 "|", cluster, "|", paste(split_var, collapse = "+")))
      mytree = lmertree(Formula, data = data, alpha = alpha_screen,
                        maxdepth = maxdepth)
      final_var = get_split_names(mytree$tree, data)
      cat("There is only one module which is grey, final features",
          final_var)
    } else {
      # at least one non-grey module
      final_var = imp_var[[1]]
      # if only one non-grey module: final_var is the chosen non-grey
      # features if at least two non-grey modules:
      if (length(imp_var) > 1) {
        # There are at least 2 modules
        for (i in 2:length(imp_var)) {
          final_var = c(final_var, imp_var[[i]])
        }

        cat("After screening within non-grey modules", final_var,
            "\n")
        # select from selected non-grey features
        regress_var = fixed_regress
        split_var = c(final_var, consider_split)
        maxdepth = ceiling(maxdepth_factor_select * length(split_var))
        Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                   "|", cluster, "|", paste(split_var, collapse = "+")))
        mytree = lmertree(Formula, data = data, alpha = alpha_select,
                          maxdepth = maxdepth)
        final_var = get_split_names(mytree$tree, data)
      }
      cat("The chosen non-grey features are", final_var, "\n")

      # use final_var (chosen non-grey features) to select grey features
      regress_var = c(fixed_regress, final_var)
      split_var = c(module_dic[["grey"]], consider_split)
      maxdepth = ceiling(maxdepth_factor_screen * length(split_var))
      Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                                 "|", cluster, "|", paste(split_var, collapse = "+")))
      mytree = lmertree(Formula, data = data, alpha = alpha_screen,
                        maxdepth = maxdepth)
      grey_var = get_split_names(mytree$tree, data)
      cat("The chosen grey features are", grey_var, "\n")
      # use final_var and grey_var do to the final model tree
      final_var = c(final_var, grey_var)
      cat("final features", final_var)
    }
    # use the final features as split&regression variables
    split_var = c(final_var, consider_split)
    maxdepth = length(split_var)
    regress_var = c(fixed_regress, final_var)
    Formula = as.formula(paste("y~", paste(regress_var, collapse = "+"),
                               "|", cluster, "|", paste(split_var, collapse = "+")))
    minsize = round(minsize_multiplier * length(regress_var))
    mytree = lmertree(Formula, data = data, alpha = alpha_predict,
                      maxdepth = maxdepth, minsize = minsize)
    mytree$final_selection = final_var
    return(mytree)
  }

}
