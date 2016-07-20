#' Generalized Boosted Regression Modeling
#' 
#' Fits generalized boosted regression models.
#' 
#' See the gbm vignette for technical details.
#' 
#' This package implements the generalized boosted modeling
#' framework. Boosting is the process of iteratively adding basis
#' functions in a greedy fashion so that each additional basis
#' function further reduces the selected loss function. This
#' implementation closely follows Friedman's Gradient Boosting Machine
#' (Friedman, 2001).
#' 
#' In addition to many of the features documented in the Gradient
#' Boosting Machine, \code{gbm} offers additional features including
#' the out-of-bag estimator for the optimal number of iterations, the
#' ability to store and manipulate the resulting \code{gbm} object,
#' and a variety of other loss functions that had not previously had
#' associated boosting algorithms, including the Cox partial
#' likelihood for censored data, the poisson likelihood for count
#' outcomes, and a gradient boosting implementation to minimize the
#' AdaBoost exponential loss function.
#' 
#' \code{gbm} is a depracated function that now acts as a front-end to \code{gbm2}
#' that uses the familiar R modeling formulas. However, \code{\link[stats]{model.frame}} is
#' very slow if there are many predictor variables. 
#' 
#' @param formula a symbolic description of the model to be fit. The
#' formula may include an offset term (e.g. y~offset(n)+x). If
#' \code{keep.data=FALSE} in the initial call to \code{gbm} then it is
#' the user's responsibility to resupply the offset to
#' \code{\link{gbm_more}}.
#' 
#' @param distribution either a character string specifying the name
#' of the distribution to use or a list with a component \code{name}
#' specifying the distribution and any additional parameters
#' needed. If not specified, \code{gbm} will try to guess: if the
#' response has only two unique values, bernoulli is assumed 
#' (two-factor responses are converted to 0,1);
#' otherwise, if the response has class "Surv", coxph is assumed;
#' otherwise, gaussian is assumed.
#' 
#' Available distributions are "Gaussian" (squared error), "Laplace"
#' (absolute loss), "TDist" (t-distribution loss), "Bernoulli"
#' (logistic regression for 0-1 outcomes), "Huberized" (Huberized
#' hinge loss for 0-1 outcomes), "AdaBoost" (the AdaBoost
#' exponential loss for 0-1 outcomes), "Poisson" (count outcomes),
#' "CoxPH" (right censored observations), "Quantile", or "Pairwise"
#' (ranking measure using the LambdaMART algorithm).
#' 
#' If quantile regression is specified, \code{distribution} must be a
#' list of the form \code{list(name="Quantile",alpha=0.25)} where
#' \code{alpha} is the quantile to estimate. Non-constant weights are
#' unsupported.
#' 
#' If "tdist" is specified, the default degrees of freedom is four and
#' this can be controlled by specifying
#' \code{distribution=list(name="tdist", df=DF)} where \code{DF} is
#' your chosen degrees of freedom.
#' 
#' If "pairwise" regression is specified, \code{distribution} must be a list of
#' the form \code{list(name="Pairwise",group=...,metric=...,max.rank=...)}
#' (\code{metric} and \code{max.rank} are optional, see below). \code{group} is
#' a character vector with the column names of \code{data} that jointly
#' indicate the group an instance belongs to (typically a query in Information
#' Retrieval applications). For training, only pairs of instances from the same
#' group and with different target labels can be considered. \code{metric} is
#' the IR measure to use, one of
#' \describe{
#'   \item{list("conc")}{Fraction of concordant pairs; for binary labels, this is equivalent to the Area under the ROC Curve}
#'   \item{:}{Fraction of concordant pairs; for binary labels, this
#' is equivalent to the Area under the ROC Curve}
#'   \item{list("mrr")}{Mean reciprocal rank of the highest-ranked positive instance}
#'   \item{:}{Mean reciprocal rank of the highest-ranked positive instance}
#'   \item{list("map")}{Mean average precision, a generalization of \code{mrr}
#'   to multiple positive instances}
#'   \item{:}{Mean average precision, a generalization of \code{mrr} to multiple positive instances}
#'   \item{list("ndcg:")}{Normalized discounted cumulative gain. The score is the weighted sum (DCG) of the user-supplied target values, weighted by
#' log(rank+1), and normalized to the maximum achievable value. This is the
#' default if the user did not specify a metric.}
#' }
#' 
#' \code{ndcg} and \code{conc} allow arbitrary target values, while binary
#' targets {0,1} are expected for \code{map} and \code{mrr}. For \code{ndcg}
#' and \code{mrr}, a cut-off can be chosen using a positive integer parameter
#' \code{max.rank}. If left unspecified, all ranks are taken into account.
#' 
#' Note that splitting of instances into training and validation sets follows
#' group boundaries and therefore only approximates the specified
#' \code{train.fraction} ratio (the same applies to cross-validation folds).
#' Internally queries are randomly shuffled before training to avoid bias.
#' 
#' Weights can be used in conjunction with pairwise metrics, however it is
#' assumed that they are constant for instances from the same group.
#' 
#' For details and background on the algorithm, see e.g. Burges (2010).
#'
#' @param data an optional data frame containing the variables in the
#' model. By default the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{gbm} is called. If \code{keep.data=TRUE} in the initial call
#' to \code{gbm} then \code{gbm} stores a copy with the object. If
#' \code{keep.data=FALSE} then subsequent calls to
#' \code{\link{gbm.more}} must resupply the same dataset. It becomes
#' the user's responsibility to resupply the same data at this point.
#'
#' @param weights an optional vector of weights to be used in the
#' fitting process. The weights must be positive but do not need to be
#' normalized. If \code{keep.data=FALSE} in the initial call to
#' \code{gbm}, then it is the user's responsibility to resupply the
#' weights to \code{\link{gbm_more}}.
#'
#' @param subset an optional vector defining a subset of the data to be used
#'
#' @param offset an optional model offset
#'
#' @param var.monotone an optional vector, the same length as the
#' number of predictors, indicating which variables have a monotone
#' increasing (+1), decreasing (-1), or arbitrary (0) relationship
#' with the outcome.
#'
#' @param n.trees the total number of trees to fit. This is equivalent
#' to the number of iterations and the number of basis functions in
#' the additive expansion.
#'
#' @param cv.folds Number of cross-validation folds to perform. If
#' \code{cv.folds}>1 then \code{gbm}, in addition to the usual fit,
#' will perform a cross-validation and calculate an estimate of
#' generalization error returned in \code{cv.error}.
#'
#' @param interaction.depth The maximum depth of variable
#' interactions: 1 builds an additive model, 2 builds a model with up
#' to two-way interactions, etc.
#'
#' @param n.minobsinnode minimum number of observations (not total
#' weights) in the terminal nodes of the trees.
#'
#' @param shrinkage a shrinkage parameter applied to each tree in the
#' expansion. Also known as the learning rate or step-size reduction.
#'
#' @param bag.fraction the fraction of independent training observations (or patients)
#' randomly selected to propose the next tree in the expansion, depending 
#' on the obs.id vector multiple training data rows may belong to
#' a single 'patient'. This introduces randomness into the model fit.
#' If \code{bag.fraction}<1 then running the same model twice will result
#' in similar but different fits. \code{gbm} uses the R random number generator, so
#' \code{set.seed} ensures the same model can be
#' reconstructed. Preferably, the user can save the returned
#' \code{\link{gbm.object}} using \code{\link{save}}.
#'
#' @param train.fraction The first \code{train.fraction * nrows(data)}
#' observations are used to fit the \code{gbm} and the remainder are used for
#' computing out-of-sample estimates of the loss function.
#'
#' @param nTrain An integer representing the number of unique patients,
#' each patient could have multiple rows associated with them, on which
#' to train.  This is the preferred way of specification for
#' \code{gbm_fit}; The option \code{train.fraction} in \code{gbm_fit}
#' is deprecated and only maintained for backward compatibility. These
#' two parameters are mutually exclusive. If both are unspecified, all
#' data is used for training.
#'
#' @param mFeatures Each node will be trained on a random subset of
#' \code{mFeatures} number of features. Each node will consider a new
#' random subset of features, adding variability to tree growth and
#' reducing computation time. \code{mFeatures} will be bounded between
#' 1 and \code{nCols}. Values outside of this bound will be set to the
#' lower or upper limits.
#'
#' @param keep.data a logical variable indicating whether to keep the
#' data and an index of the data stored with the object. Keeping the
#' data and index makes subsequent calls to \code{\link{gbm_more}}
#' faster at the cost of storing an extra copy of the dataset.
#'
#' @param object a \code{gbm} object created from an initial call to
#' \code{\link{gbm}}.
#'
#' @param n.new.trees the number of additional trees to add to
#' \code{object} using \code{gbm.more}.
#'
#' @param verbose If TRUE, gbm will print out progress and performance
#' indicators. If this option is left unspecified for gbm.more then it uses
#' \code{verbose} from \code{object}.
#' 
#' @param class.stratify.cv whether the cross-validation should be
#' stratified by class. Is only implemented for
#' \code{bernoulli}. The purpose of stratifying
#' the cross-validation is to help avoiding situations in which
#' training sets do not contain all classes.
#'
#' @param x,y For \code{gbm_fit}: \code{x} is a data frame or data
#' matrix containing the predictor variables and \code{y} is a matrix of outcomes.
#' Excluding \code{CoxPH} this matrix of outcomes collapses to a vector, in the case of
#' \code{CoxPH} it is a survival object where the event times fill the first 
#' one (or two columns) and the status fills the final column. 
#' The number of rows in \code{x} must be the same as the length of the 1st dimension
#' of \code{y}.
#'
#' @param misc For \code{gbm_fit}: \code{misc} is an R object that is
#' simply passed on to the gbm engine. It can be used for additional
#' data for the specific distribution.  Currently it is only used for
#' with the Cox proportional hazards model to set the tied times method. 
#'
#' @param w For \code{gbm.fit}: \code{w} is a vector of weights of the same
#' length as the 1st dimension of \code{y}.
#'
#' @param var.names For \code{gbm.fit}: A vector of strings of length
#' equal to the number of columns of \code{x} containing the names of
#' the predictor variables.
#'
#' @param response.name For \code{gbm.fit}: A character string label
#' for the response variable.
#'
#' @param group \code{group} used when \code{distribution =
#' 'pairwise'.}
#' 
#' @param tied.times.method For \code{gbm}: This is an optional string used with
#' \code{CoxPH} distribution specifying what method to employ when dealing with tied times. 
#' Currently only "efron" and "breslow" are available; the default value is "efron". 
#' Setting the string to any other value reverts the method to the original CoxPH model
#' implementation where ties are not explicitly dealt with.
#' 
#' @param prior.node.coeff.var Optional double only used with the \code{CoxPH} distribution.
#' It is a prior on the coefficient of variation associated with the hazard 
#' rate assigned to each terminal node when fitting a tree.Increasing its value emphasises 
#' the importance of the training data in the node when assigning a prediction to said node.
#'
#' @param strata Optional vector of integers only used with the \code{CoxPH} distributions.
#' Each integer in this vector represents the stratum the corresponding row in the data belongs to,
#' e. g. if the 10th element is 3 then the 10th data row belongs to the 3rd strata.
#' 
#' @param obs.id Optional vector of integers used to specify which rows of data belong
#' to individual patients.  Data is then bagged by patient id; the default sets each row of the data
#' to belong to an individual patient.
#' 
#' @param n.cores The number of CPU cores to use. The cross-validation
#' loop will attempt to send different CV folds off to different
#' cores. If \code{n.cores} is not specified by the user, it is
#' guessed using the \code{detectCores} function in the
#' \code{parallel} package. Note that the documentation for
#' \code{detectCores} makes clear that it is not reliable and could
#' return a spurious number of available cores.
#'
#' @param fold.id An optional vector of values identifying what fold
#' each observation is in. If supplied, cv.folds can be missing. Note:
#' Multiple observations to the same patient must have the same fold id.
#' 
#' 
#' @usage
#' gbm(formula = formula(data), distribution = "bernoulli",
#' data = list(), weights, subset = NULL, offset = NULL, var.monotone
#' = NULL, n.trees = 100, interaction.depth = 1, n.minobsinnode = 10,
#' shrinkage = 0.001, bag.fraction = 0.5, train.fraction = 1,
#' mFeatures = NULL, cv.folds = 0, keep.data = TRUE, verbose = FALSE,
#' class.stratify.cv = NULL, n.cores = NULL, fold.id=NULL,
#' tied.times.method = "efron", prior.node.coeff.var = 1000, strata = NULL, 
#' obs.id = 1:nrow(data))
#' 
#' gbm.fit(x, y, offset = NULL, misc = NULL, distribution = "bernoulli", 
#' w = NULL, var.monotone = NULL, n.trees = 100, interaction.depth = 1, 
#' n.minobsinnode = 10, shrinkage = 0.001, bag.fraction = 0.5, 
#' nTrain = NULL, train.fraction = NULL, mFeatures = NULL, keep.data = TRUE, 
#' verbose = TRUE, var.names = NULL, response.name = "y", group = NULL,
#' prior.node.coeff.var = 1000, strata = NULL, obs.id = 1:nrow(x))
#'
#'
#' @return \code{gbm} and \code{gbm.fit} return a
#' \code{GBMFit} object.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' 
#' Quantile regression code developed by Brian Kriegler
#' \email{bk@@stat.ucla.edu}
#' 
#' t-distribution code developed by Harry Southworth and
#' Daniel Edwards
#' 
#' Pairwise code developed by Stefan Schroedl \email{schroedl@@a9.com}
#' @seealso \code{\link{gbm2}}, \code{\link{perf_gbm}},
#' \code{\link{plot.GBMFit}}, \code{\link{predict.GBMFit}},
#' \code{\link{summary.GBMFit}}, \code{\link{pretty.gbm.tree}}.
#' @references Y. Freund and R.E. Schapire (1997) \dQuote{A decision-theoretic
#' generalization of on-line learning and an application to boosting,}
#' \emph{Journal of Computer and System Sciences,} 55(1):119-139.
#' 
#' G. Ridgeway (1999). \dQuote{The state of boosting,} \emph{Computing Science
#' and Statistics} 31:172-181.
#' 
#' J.H. Friedman, T. Hastie, R. Tibshirani (2000). \dQuote{Additive Logistic
#' Regression: a Statistical View of Boosting,} \emph{Annals of Statistics}
#' 28(2):337-374.
#' 
#' J.H. Friedman (2001). \dQuote{Greedy Function Approximation: A Gradient
#' Boosting Machine,} \emph{Annals of Statistics} 29(5):1189-1232.
#' 
#' J.H. Friedman (2002). \dQuote{Stochastic Gradient Boosting,}
#' \emph{Computational Statistics and Data Analysis} 38(4):367-378.
#' 
#' B. Kriegler (2007).
#' \href{http://statistics.ucla.edu/theses/uclastat-dissertation-2007:2}{Cost-Sensitive
#' Stochastic Gradient Boosting Within a Quantitative Regression Framework}.
#' PhD dissertation, UCLA Statistics.
#' 
#' C. Burges (2010). \dQuote{From RankNet to LambdaRank to LambdaMART: An
#' Overview,} Microsoft Research Technical Report MSR-TR-2010-82.
#' 
#' \href{http://sites.google.com/site/gregridgeway}{Greg Ridgeway's site}.
#' 
#' The \href{http://www-stat.stanford.edu/~jhf/R-MART.html}{MART} website.
#' @keywords models nonlinear survival nonparametric tree
#' @examples
#'  # A least squares regression example # create some data
#' 
#' N <- 1000
#' X1 <- runif(N)
#' X2 <- 2*runif(N)
#' X3 <- ordered(sample(letters[1:4],N,replace=TRUE),levels=letters[4:1])
#' X4 <- factor(sample(letters[1:6],N,replace=TRUE))
#' X5 <- factor(sample(letters[1:3],N,replace=TRUE))
#' X6 <- 3*runif(N) 
#' mu <- c(-1,0,1,2)[as.numeric(X3)]
#' 
#' SNR <- 10 # signal-to-noise ratio
#' Y <- X1**1.5 + 2 * (X2**.5) + mu
#' sigma <- sqrt(var(Y)/SNR)
#' Y <- Y + rnorm(N,0,sigma)
#' 
#' # introduce some missing values
#' X1[sample(1:N,size=500)] <- NA
#' X4[sample(1:N,size=300)] <- NA
#' 
#' data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
#' 
#' # fit initial model
#' gbm1 <-
#' gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
#'     data=data,                   # dataset
#'     var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease,
#'                                  # +1: monotone increase,
#'                                  #  0: no monotone restrictions
#'     distribution="Gaussian",     # see the help for other choices
#'     n.trees=1000,                # number of trees
#'     shrinkage=0.05,              # shrinkage or learning rate,
#'                                  # 0.001 to 0.1 usually work
#'     interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
#'     bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
#'     train.fraction = 0.5,        # fraction of data for training,
#'                                  # first train.fraction*N used for training
#'     mFeatures = 3,               # half of the features are considered at each node
#'     n.minobsinnode = 10,         # minimum total weight needed in each node
#'     cv.folds = 3,                # do 3-fold cross-validation
#'     keep.data=TRUE,              # keep a copy of the dataset with the object
#'     verbose=FALSE,               # don't print out progress
#'     n.cores=1)                   # use only a single core (detecting #cores is
#'                                  # error-prone, so avoided here)
#' 
#' # check performance using an out-of-bag estimator
#' # OOB underestimates the optimal number of iterations
#' best.iter <- gbm.perf(gbm1,method="OOB")
#' print(best.iter)
#' 
#' # check performance using a 50% heldout test set
#' best.iter <- gbm.perf(gbm1,method="test")
#' print(best.iter)
#' 
#' # check performance using 3-fold cross-validation
#' best.iter <- gbm.perf(gbm1,method="cv")
#' print(best.iter)
#' 
#' # plot the performance # plot variable influence
#' summary(gbm1,n.trees=1)         # based on the first tree
#' summary(gbm1,n.trees=best.iter) # based on the estimated best number of trees
#' 
#' # compactly print the first and last trees for curiosity
#' print(pretty.gbm.tree(gbm1,1))
#' print(pretty.gbm.tree(gbm1,gbm1$n.trees))
#' 
#' # make some new data
#' N <- 1000
#' X1 <- runif(N)
#' X2 <- 2*runif(N)
#' X3 <- ordered(sample(letters[1:4],N,replace=TRUE))
#' X4 <- factor(sample(letters[1:6],N,replace=TRUE))
#' X5 <- factor(sample(letters[1:3],N,replace=TRUE))
#' X6 <- 3*runif(N) 
#' mu <- c(-1,0,1,2)[as.numeric(X3)]
#' 
#' Y <- X1**1.5 + 2 * (X2**.5) + mu + rnorm(N,0,sigma)
#' 
#' data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
#' 
#' # predict on the new data using "best" number of trees
#' # f.predict generally will be on the canonical scale (logit,log,etc.)
#' f.predict <- predict(gbm1,data2,best.iter)
#' 
#' # least squares error
#' print(sum((data2$Y-f.predict)^2))
#' 
#' # create marginal plots
#' # plot variable X1,X2,X3 after "best" iterations
#' par(mfrow=c(1,3))
#' plot(gbm1,1,best.iter)
#' plot(gbm1,2,best.iter)
#' plot(gbm1,3,best.iter)
#' par(mfrow=c(1,1))
#' # contour plot of variables 1 and 2 after "best" iterations
#' plot(gbm1,1:2,best.iter)
#' # lattice plot of variables 2 and 3
#' plot(gbm1,2:3,best.iter)
#' # lattice plot of variables 3 and 4
#' plot(gbm1,3:4,best.iter)
#' 
#' # 3-way plots
#' plot(gbm1,c(1,2,6),best.iter,cont=20)
#' plot(gbm1,1:3,best.iter)
#' plot(gbm1,2:4,best.iter)
#' plot(gbm1,3:5,best.iter)
#' 
#' # do another 100 iterations
#' gbm2 <- gbm.more(gbm1,100,
#'                  verbose=FALSE) # stop printing detailed progress
#' @export
gbm <- function(formula = formula(data),
                distribution = "Bernoulli",
                data = list(),
                weights,
                subset = NULL,
                offset = NULL,
                var.monotone = NULL,
                n.trees = 100,
                interaction.depth = 1,
                n.minobsinnode = 10,
                shrinkage = 0.001,
                bag.fraction = 0.5,
                train.fraction = 1.0,
                mFeatures = NULL,
                cv.folds=0,
                keep.data = TRUE,
                verbose = FALSE,
                class.stratify.cv=NULL,
                n.cores=NULL,
                fold.id = NULL,
                tied.times.method="efron",
                prior.node.coeff.var=1000,
                strata=NULL, obs.id=1:nrow(data)) {
  # Highlight new API
  warning("gbm is depracated - using gbm2...")
  
  # Unravel model 
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "subset", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  m <- mf
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  
  w <- model.weights(mf)
  offset <- model.offset(mf)
  var.names <- attributes(Terms)$term.labels
  
  # Set distribution object - put in all possible additional parameters (this will generate warnings)
  if (is.character(distribution)){ distribution <- list(name=distribution) }
  dist_obj <- gbm_dist(distribution$name, ties=tied.times.method, strata=strata,
                      group=distribution$group, metric=distribution$metric,
                      max.rank=distribution$max.rank, prior_node_coeff_var=prior.node.coeff.var,
                      alpha=distribution$alpha, df=distribution$df, power=distribution$power)
  
  # Set up training parameters
  if(is.null(mFeatures)) mFeatures <- ncol(data)-1
  params <- training_params(num_trees=n.trees, interaction_depth=interaction.depth, min_num_obs_in_node=n.minobsinnode, 
                           shrinkage=shrinkage, bag_fraction=bag.fraction, id=obs.id,
                           num_train=round(train.fraction * length(unqiue(obs.id))),
                           num_features=mFeatures)
  
  # Call gbm2
  gbm_fit_obj <- gbm2(formula, distribution=dist_obj, data, weights=w, offset=offset,
                     train_params= params, var_monotone=var.monotone, var_names=var.names,
                     cv_folds=cv.folds, cv_class_stratify=class.stratify.cv, fold_id=fold.id,
                     keep_gbm_data=keep.data, is_verbose=verbose)

  return(gbm_fit_obj)
}
