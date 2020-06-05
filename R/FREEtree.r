#' Initial FREEtree call which then calls actual FREEtree methods depending
#' on parameters being passed through.
#'
#' @param data data to train or test FREEtree on.
#' @param fixed_regress user specified char vector of regressors that will never
#'   be screened out; if fixed_regress = NULL, method uses PC as regressor at
#'   screening step.
#' @param fixed_split user specified char vector of features to be used in
#'   splitting with certainty.
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
#' @param cluster the variable name of each cluster (in terms of random effect) using glmer's implementation.
#' @param minModuleSize WGCNA's minimum module size parameter.
#' @importFrom stats as.formula
#' @return a glmertree object (trained tree).
#' @examples
#' #locate example data file
#' load("data/data.RData")
#' mytree = FREEtree(data,fixed_regress=c("time","time2"), fixed_split=c("treatment"),
#'                   var_select=paste("V",1:400,sep=""), minModuleSize = 5,
#'                   cluster="patient", Fuzzy=TRUE, maxdepth_factor_select = 0.5,
#'                   maxdepth_factor_screen = 0.04, minsize_multiplier = 5,
#'                   alpha_screen = 0.2, alpha_select=0.2,alpha_predict=0.05)
#'
#' @export
FREEtree <- function(data, fixed_regress = NULL, fixed_split = NULL, var_select = NULL,
                    power = 6, minModuleSize = 1, cluster, maxdepth_factor_screen = 0.04,
                    maxdepth_factor_select = 0.5, Fuzzy = TRUE, minsize_multiplier = 5,
                    alpha_screen = 0.2, alpha_select = 0.2, alpha_predict = 0.05) {
  ### if there are no features to select, just use fixed_regress and
  ### fixed_split ### no need to filter, just use GLMM tree right away
  if (length(var_select) == 0)
  {
    if (length(fixed_regress) == 0) {
      if (length(fixed_split) == 0) {
        stop("no features to split and regress on")
      }
      fixed_regress = "1"
    }
    maxdepth = length(fixed_split)
    Formula = as.formula(paste("y~", paste(fixed_regress, collapse = "+"),
                               "|", cluster, "|", paste(fixed_split, collapse = "+")))
    mytree = lmertree(Formula, data = data, alpha = alpha_predict,
                      maxdepth = maxdepth)
    mytree$final_selection = NULL
    return(mytree)
  }  ###
  # Now var_select is not empty If fixed_regress (eg: time) is not
  # specified: use PCs as regressors at screening step
  if (length(fixed_regress) == 0)
  {
    cat("Use FREEtree_PC\n")
    return(FREEtree_PC(data = data, fixed_split = fixed_split,
                       var_select = var_select, power = power, minModuleSize = minModuleSize,
                       cluster = cluster, maxdepth_factor_screen = maxdepth_factor_screen,
                       maxdepth_factor_select = maxdepth_factor_select, Fuzzy = Fuzzy,
                       minsize_multiplier = minsize_multiplier, alpha_screen = alpha_screen,
                       alpha_select = alpha_select, alpha_predict = alpha_predict))
  }  ###
  # Now we have non-empty var_select,fixed_regress
  cat("Use FREEtree_time\n")
  return(FREEtree_time(data = data, fixed_regress = fixed_regress, fixed_split = fixed_split,
                       var_select = var_select, power = power, minModuleSize = minModuleSize,
                       cluster = cluster, maxdepth_factor_screen = maxdepth_factor_screen,
                       maxdepth_factor_select = maxdepth_factor_select, Fuzzy = Fuzzy,
                       minsize_multiplier = minsize_multiplier, alpha_screen = alpha_screen,
                       alpha_select = alpha_select, alpha_predict = alpha_predict))

}


