#' Doubly robust scores for estimating the average conditional local average treatment effect.
#'
#' Given an outcome Y, treatment W and instrument Z, the (conditional) local
#' average treatment effect is tau(x) = Cov[Y, Z | X = x] / Cov[W, Z | X = x].
#' This is the quantity that is estimated with an instrumental forest.
#' It can be intepreted causally in various ways. Given a homogeneity
#' assumption, tau(x) is simply the CATE at x. When W is binary
#' and there are no "defiers", Imbens and Angrist (1994) show that tau(x) can
#' be interpreted as an average treatment effect on compliers. This doubly robust
#' scores provided here are for estimating tau = E[tau(X)].
#'
#' @param forest A trained instrumental forest.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param debiasing.weights A vector of length n (or the subset length) of debiasing weights.
#'               If NULL (default) these are obtained via the appropriate doubly robust score
#'               construction, e.g., in the case of causal_forests with a binary treatment, they
#'               are obtained via inverse-propensity weighting.
#' @param compliance.score An estimate of the causal effect of Z on W, i.e., Delta(X) = E[W | X, Z = 1]
#'               - E[W | X, Z = 0], which can then be used to produce debiasing.weights. If not provided,
#'               this is estimated via an auxiliary causal forest.
#' @param num.trees.for.weights In some cases (e.g., with causal forests with a continuous
#'               treatment), we need to train auxiliary forests to learn debiasing weights.
#'               This is the number of trees used for this task. Note: this argument is only
#'               used when debiasing.weights = NULL.
#' @param ... Additional arguments (currently ignored).
#'
#' @references Aronow, Peter M., and Allison Carnegie. "Beyond LATE: Estimation of the
#'              average treatment effect with an instrumental variable." Political
#'              Analysis 21(4), 2013.
#' @references Chernozhukov, Victor, Juan Carlos Escanciano, Hidehiko Ichimura,
#'             Whitney K. Newey, and James M. Robins. "Locally robust semiparametric
#'             estimation." Econometrica 90(4), 2022.
#' @references Imbens, Guido W., and Joshua D. Angrist. "Identification and Estimation of
#'             Local Average Treatment Effects." Econometrica 62(2), 1994.
#'
#' @return A vector of scores.
#'
#' @exportS3Method grf::get_scores
get_scores.causal_iv_forest <- function(forest,
                                           subset = NULL,
                                           debiasing.weights = NULL,
                                           compliance.score = NULL,
                                           num.trees.for.weights = 500,
                                           ...) {
  if (!all(forest$Z.orig %in% c(0, 1))) {
    stop(paste(
      "Average conditional local average treatment effect estimation",
      "only implemented for binary instruments."
    ))
  }
  subset <- validate_subset(forest, subset)

  W.orig <- forest$W.orig[subset]
  W.hat <- forest$W.hat[subset]
  Y.orig <- forest$Y.orig[subset]
  Y.hat <- forest$Y.hat[subset]
  Z.orig <- forest$Z.orig[subset]
  Z.hat <- forest$Z.hat[subset]

  if (is.null(debiasing.weights)) {
  # The compliance forest estimates the effect of the "treatment" Z on the "outcome" W.
    if (is.null(compliance.score)) {
      clusters <- if (length(forest$clusters) > 0) {
        forest$clusters
      } else {
        1:length(forest$Y.orig)
      }
      compliance.forest <- grf::causal_forest(X = forest$X.orig,
                                         Y = forest$W.orig,
                                         W = forest$Z.orig,
                                         Y.hat = forest$W.hat,
                                         W.hat = forest$Z.hat,
                                         sample.weights = forest$sample.weights,
                                         clusters = clusters,
                                         num.trees = num.trees.for.weights,
                                         seed = forest$seed,
                                         num.threads = forest$num.threads)
      compliance.score <- predict(compliance.forest)$predictions
      compliance.score <- compliance.score[subset]
    } else if (length(compliance.score) == length(forest$Y.orig)) {
      compliance.score <- compliance.score[subset]
    } else if (length(compliance.score) != length(subset))  {
      stop("If specified, compliance.score must have length n or |subset|.")
    }
    if (min(Z.hat) <= 0.01 || max(Z.hat) >= 0.99) {
      rng <- range(Z.hat)
      warning(paste0(
        "Estimated instrument propensities take values between ",
        round(rng[1], 3), " and ", round(rng[2], 3),
        " and in particular get very close to 0 or 1. ",
        "Poor overlap may hurt perfmance for average conditional local average ",
        "treatment effect estimation."
      ))
    }
    if (min(abs(compliance.score)) <= 0.01 * sd(W.orig)) {
      warning(paste0(
        "The instrument appears to be weak, with some compliance scores as ",
        "low as ", round(min(compliance.score), 4)
      ))
    }
    debiasing.weights <- (Z.orig - Z.hat) / (Z.hat * (1 - Z.hat)) / compliance.score
  } else if (length(debiasing.weights) == length(forest$Y.orig)) {
    debiasing.weights <- debiasing.weights[subset]
  } else if (length(debiasing.weights) != length(subset))  {
    stop("If specified, debiasing.weights must have length n or |subset|.")
  }

  tau.hat.pointwise <- predict(forest)$predictions[subset]
  Y.residual <- Y.orig - (Y.hat + tau.hat.pointwise * (W.orig - W.hat))

  tau.hat.pointwise + debiasing.weights * Y.residual
}
