#### Preprocessing ####
rescaleY <- function(y) {rr <- range(y); (y - rr[1L]) / diff(rr) - 0.5}

#### Tree representation ####
aTree <- list(
  variable = c("age", NA, "height", NA, NA),
  rule = c(23, NA, 168, NA, NA),
  expectedY = c(NA, 2, NA, 3, 2.5),
  children = list(c(2L, 3L), NA, c(4L, 5L), NA, NA),
  parent = c(0L, 1L, 1L, 3L, 3L)
)
predictTree <- function(tree, x, component) {
  # sanity check
  stopifnot(
    all( c("variable", "rule", "expectedY", "children") %in% names(tree) ),
    identical(is.na(tree$variable), is.na(tree$rule)),
    identical(is.na(tree$variable), !is.na(tree$expectedY)),
    identical(is.na(tree$variable), is.na(tree$children)),
    inherits(x, "data.frame"),
    all( na.omit(tree$variable) %in% names(x) )
  )

  # follow the thread
  followThread <- function(covs, node = 1L) {
    # Encounters end node
    if (is.na(tree$variable[node])) {
      tree$expectedY[[node]][component]
    } else {
      if (covs[tree$variable[node]] <= tree$rule[node]) # Going left
        followThread(covs, tree$children[[node]][1L])
      else { # Going right
        followThread(covs, tree$children[[node]][2L])
      }
    }
  }
  apply(x, 1, followThread)
}

#### One tree at a time ####
b <- function(tree) sum(is.na(tree$variable))
obs_per_node <- function(tree, x) {
  output <- vector(mode = "list", length = length(tree$variable))
  # follow the thread
  followThread <- function(row, node = 1L) {
    # Encounters end node
    if (is.na(tree$variable[node])) {
      list(node = node, obs = row)
    } else {
      df_row <- x[row, ]
      if (df_row[tree$variable[node]] <= tree$rule[node]) # Going left
        followThread(row, tree$children[[node]][1L])
      else { # Going right
        followThread(row, tree$children[[node]][2L])
      }
    }
  }
  for (i in seq(nrow(x))) {
    nn <- followThread(i)
    output[[nn$node]] <- c(output[[nn$node]], nn$obs)
  }
  output
}

M_LN_SQRT_2PI <- 0.5*log(2*pi)

treeMetropolis <- function(tree, x, yb, eta, R1, R0, a_beta, b_beta , sigma21, sigma20, proposal = list(
  grow = 0.25, prune = 0.25, change = 0.4, swap = 0.1
), prior = list(alpha = 0.95, beta = 2)) {
  possibleProposals <- c("grow", "prune", "change", "swap")

  stopifnot(is.list(proposal),
            all(possibleProposals %in% names(proposal)),
            is.list(prior),
            all(c("alpha", "beta") %in% names(prior)),
            all(prior$alpha > 0, prior$alpha < 1, prior$beta >= 0))

  n <- nrow(x)
  # Single node trees can only grow
  possibleProposals <- "grow"
  if (length(tree$variable) > 1L) # Having more than 1 node opens prune and change
    possibleProposals <- c(possibleProposals, "prune", "change")
  # swap is only available if there are children that are not leafs
  if (any(tree$parent[tree$parent] > 0L))
    possibleProposals <- c(possibleProposals, "swap")

  b <- b(tree)

  # likelihood where expected Y has been integrated out
  marginalizedLikelihood <- function(tree) {
    dividedObs <- obs_per_node(tree, x)
    n_ell <- sapply(dividedObs, length)
    # Remove non-terminal nodes
    dividedObs <- dividedObs[n_ell > 0]
    n_ell <- n_ell[n_ell > 0]
    #    y1 <- lapply(dividedObs, function(i) na.omit(y$y1[i]) )
    yb <- lapply(dividedObs, function(i) na.omit(yb[i]) )
    #    y0 <- lapply(dividedObs, function(i) na.omit(y$y0[i]) )

    etab <- lapply(dividedObs, function(i) na.omit(eta[i]) )
    R1d <- lapply(dividedObs, function(i) na.omit(R1[i]) )
    R0d <- lapply(dividedObs, function(i) na.omit(R0[i]) )

    n1 <- lengths(R1d)
    nb <- lengths(yb)
    n0 <- lengths(R0d)

    # 0.5 / sigma2 * sum(sapply(ys, sum) ^ 2 /
    #                          (sigma2 / sigma2_mu + n_ell)) -
    #   0.5 * sum(log(sigma2 + sigma2_mu * n_ell))

    # beta part

    output <- 0
    for(l in 1:length(n_ell)){
      if(n1[l]>0){
        a <- 1/sigma21
        aa <- a / (n1[l] + a)
        R_bar <- mean(R1d[[l]])
        SSE_Z <- sum((R1d[[l]]-R_bar)^2)

        Ltheta1 <-  0.5 * log(aa) - n1[l] * M_LN_SQRT_2PI
        - 0.5 * (SSE_Z + n1[l] * aa * R_bar * R_bar);

        output <- output + Ltheta1
      }
      if(nb[l]>0){
        ss.sum_v_logY <- sum(etab[[l]]*log(yb[[l]]))
        Lu <- lgamma(a_beta + nb[l])
        - (a_beta + nb[l]) * log(b_beta - ss.sum_v_logY);

        output <- output + Lu
      }
      if(n0[l]>0){
        a <- 1/sigma20
        aa <- a / (n0[l] + a)
        R_bar <- mean(R0d[[l]])
        SSE_Z <- sum((R0d[[l]]-R_bar)^2)

        Ltheta0 <-  0.5 * log(aa) - n0[l] * M_LN_SQRT_2PI
        - 0.5 * (SSE_Z + n0[l] * aa * R_bar * R_bar);

        output <- output + Ltheta0
      }

    }
    output
  }

  # For the tree prior
  findDepth <- function(tree, node) {
    depth <- 0
    currentParent <- tree$parent[node]
    while (currentParent > 0) {
      depth <- depth + 1
      currentParent <- tree$parent[currentParent]
    }
    depth
  }
  nonTerminalProb <- function(d) prior$alpha * (1 + d) ^ (-prior$beta)

  # Simplify pruning - find leaves' parents
  prunableNodes <- function(tree) {
    leafs <- is.na(tree$variable)
    unique(tree$parent[leafs])
  }

  # Viability check as per Chipman (1998)
  isViable <- function(tree) {
    allObs <- obs_per_node(tree, x)
    terminals <- which(is.na(tree$variable))
    allObs <- allObs[terminals]
    min(sapply(allObs, length)) >= 5L
  }

  # This is used to calculate prior & proposal ratio for grow and prune
  growprune <- function(newTree, tree, depth, grow = TRUE) {
    (log(proposal$prune) - log(proposal$grow) +
       log(b - !grow) - log(length(prunableNodes(newTree))) +
       log(prior$alpha) + 2 * log1p(-nonTerminalProb(depth + 1)) -
       log((1 + depth) ^ prior$beta - prior$alpha)) * (-1) ^ grow
  }

  foundViableTree <- FALSE
  treeN <- length(tree$variable)

  while (!foundViableTree) {
    # Sample new proposal step every time. If not, there can be infinite loops.
    proposalStep <-
      sample(possibleProposals, 1L,
             prob = as.numeric(proposal[seq_along(possibleProposals)]))
    newTree <- tree # This will be the proposed tree.
    MHprob <- -marginalizedLikelihood(tree)
    if (proposalStep == "grow") {
      ## Select a terminal node and variable.
      node <- which(cumsum(is.na(tree$variable)) == sample.int(b, 1L))[1L]
      variable <- sample(names(x), 1L)

      ## Create proposal tree
      # Turn previous terminal node into non-terminal.
      newTree$variable[node] <- variable
      ux <- unique(x[, variable])
      newTree$rule[node] <- sample(ux, 1L)
      newTree$expectedY[node] <- NA
      newTree$children[[node]] <- treeN + (1L:2L)
      # Add new terminal nodes in the end.
      newTree$variable <- c(newTree$variable, NA, NA)
      newTree$rule <- c(newTree$rule, NA, NA)
      newTree$expectedY <- c(newTree$expectedY, 0, 0) # Will update later.
      newTree$children <- c(newTree$children, NA, NA)
      newTree$parent <- c(newTree$parent, node, node)

      ## There's a bunch of cancellations between priors.
      nodeDepth <- findDepth(tree, node)
      MHprob <- MHprob + growprune(newTree, tree, nodeDepth)
    } else if (proposalStep == "prune") {
      ## Select a pair of terminal nodes to cut off. For this, the parent must
      ## not be a grandparent. Take the parents of all terminal nodes. Whoever
      ## appears twice is eligible.
      possibleParents <- as.integer(
        names(which(table(tree$parent[is.na(tree$variable)]) == 2L))
      )
      parent <- if (length(possibleParents) > 1L)
        sample(possibleParents, 1L) else possibleParents
      # Children are sorted. Reversed are decreasing.
      # It's easier to remove them without index accidents.
      allChildren <- rev(tree$children[[parent]])

      ## Create proposal tree
      # Change parent into leaf
      newTree$variable[parent] <- NA
      newTree$rule[parent] <- NA
      newTree$expectedY[parent] <- 0 # Will update later
      newTree$children[[parent]] <- NA
      # Delete previous children. This is tricky due to indexing.
      for (child in allChildren) {
        newTree$children <- lapply(newTree$children, function(ch) {
          # Subtract 1 or 0 depending on if that index is bigger than child
          ch - (ch > child)
        })
        newTree$parent <- newTree$parent - (newTree$parent > child)
        newTree$variable <- newTree$variable[-child]
        newTree$rule <- newTree$rule[-child]
        newTree$expectedY <- newTree$expectedY[-child]
        newTree$children <- newTree$children[-child]
        newTree$parent <- newTree$parent[-child]
      }
      nodeDepth <- findDepth(tree, parent)
      MHprob <- MHprob + growprune(newTree, tree, nodeDepth, FALSE)
    } else if (proposalStep == "change") {
      ## Select a non-terminal node.
      node <- which(cumsum(is.na(tree$expectedY)) ==
                      sample.int(length(tree$variable) - b, 1L))[1L]

      ## Create proposal tree.
      xNames <- names(x)
      dontRepeatRule <- TRUE
      while (dontRepeatRule) {
        variable <- sample(xNames, 1L)
        ux <- unique(x[, variable])
        rule <- sample(ux, 1L)
        if (any(variable != tree$variable[tree$parent[node]], ## Different rule from parent
                rule != tree$rule[tree$parent[node]])) dontRepeatRule <- FALSE
        else {
          dontRepeatRule <- FALSE
          for (child in tree$children[[node]]) {
            if (is.na(tree$variable[child])) next
            if (all(identical(variable, tree$variable[child]),
                    identical(rule, tree$rule[child]))) dontRepeatRule <- TRUE
          }
        }
      }
      newTree$variable[node] <- variable
      newTree$rule[node] <- rule

      # No need to update MHprob?
    } else if (proposalStep == "swap") {
      ## Choose random child/parent pair, where the former is not terminal
      eligibleChildren <- is.na(tree$expectedY) & (tree$parent != 0L)
      child <- which(cumsum(eligibleChildren) ==
                       sample.int(sum(eligibleChildren), 1L))[1L]
      parent <- tree$parent[child]
      otherChild <- tree$children[[parent]][tree$children[[parent]] != child]

      ## Create proposal tree.
      pair <- c(child, parent); revPair <- rev(pair)
      uxOldParent <- uxOldChild <- 0 # Will take effect in MHprob update if true.
      for (element in c("variable", "rule")) {
        # If the other child has the same
        if (identical(tree$variable[otherChild], tree$variable[child])) {
          newTree[[element]][otherChild] <- tree[[element]][parent]
          uxOldParent <- log(length(unique(x[, tree$variable[parent]])))
          uxOldChild <- log(length(unique(x[, tree$variable[child]])))
        }
        newTree[[element]][pair] <- tree[[element]][revPair]
      }
      MHprob <- MHprob + (uxOldParent - uxOldChild)
    }
    foundViableTree <- isViable(newTree)
  }
  MHprob <- MHprob + marginalizedLikelihood(newTree)

  if (log(runif(1L)) <= MHprob) newTree else tree
}

updateOneMub <- function(y, prior, eta) {
  lambda <- rgamma(1,length(y)+prior$a_beta,prior$b_beta-sum(eta*log(y), na.rm = TRUE))
  if (length(eta) != length(y)) print(list(eta, y))
  log(lambda)
}

updateOneMu01 <- function(r, prior) {
  var <- (1 / prior$variance + length(r) / 1) ^ (-1)
  rnorm(1L, var * (prior$mean / prior$variance + sum(r, na.rm = TRUE) / 1), sqrt(var))
}

updateMus <- function(tree, x, yb, r1, r0, eta, priorMu1, priorMub, priorMu0) {
  decisions <- obs_per_node(tree, x)
  for (i in seq_along(tree$expectedY)) {
    if (is.na(tree$variable[i])){
      tree$expectedY[[i]] <- numeric(3)
      tree$expectedY[[i]][1] <- updateOneMu01(r1, priorMu1)
      tree$expectedY[[i]][2] <- updateOneMub(yb[decisions[[i]]], priorMub, eta[decisions[[i]]])
      tree$expectedY[[i]][3] <- updateOneMu01(r0, priorMu0)
    }
    else
      tree$expectedY[[i]] <- NA
  }
  tree
}

updateSigma2 <- function(trees, x, y, prior) {
  predicted <- rowSums(do.call(cbind, lapply(trees, predictTree, x = x)))
  1 / rgamma(
    1L,
    length(y) / 2 + prior$alpha,
    crossprod(y - predicted) / 2 + prior$beta
  )
}

#' @export
beta_bart <- function(x, y, trees = 200,
                      priorMub = list(a_beta = 2, b_beta = 1), # High prior prob that mub in (-5, 5) and mode = 0
                      priorMu1 = list(mean = 0, variance = 9 / (4 * trees)), # Chipman (2010)
                      priorMu0 = list(mean = 0, variance = 9 / (4 * trees)), # Chipman (2010)
                      mcmc = list(burnin = 1000, sample = 5000)) {
  stopifnot(is.list(mcmc), "sample" %in% names(mcmc))

  # # User can define a list or a pair of numbers
  # if (!is.list(priorS2)) {
  #   topS2 <- var(residuals(lm(y ~ ., data = data.frame(x = x, y = y))))
  #   beta <- uniroot(function(b)
  #     pgamma(1 / topS2, priorS2[1L] / 2, b, lower.tail = FALSE) -
  #       priorS2[2L], c(0, 1000000))$root
  #   priorS2 <- list(alpha = priorS2[1L] / 2, beta = beta)
  # }

  # y1b0 selection
  y.mod <- function(y){
    y1 <- 1*(y==1)
    yb <- ifelse(!y%in%c(0,1),y,NA)
    y0 <- ifelse(y==1,NA,1*(y==0))
    data.frame(y1,yb,y0)
  }

  ## This is for Backfitting : R?
  getR <- function(trees, m, z, onezero) { # m is m-th tree
    output <- z
    if (length(trees) > 1) {
      output <- output -
        rowSums(do.call(cbind, lapply(trees[-m], predictTree, x = x, component = onezero)))
    }
    output
  }

  getEta <- function(trees, m) { # m is m-th tree
    output <- 0
    if (length(trees) > 1) {
      output <- output +
        rowSums(do.call(cbind, lapply(trees[-m], predictTree, x = x, component = 2)))
    }
    exp(output)
  }

  updateOneTree <- function(trees, x, m, z1, z0, yb, priorMu1, priorMub, priorMu0 ) {
    r1 <- getR(trees, m, z1, 1)
    r0 <- getR(trees, m, z0, 3)
    eta <- getEta(trees, m)
    eta[is.na(yb)] <- NA
    tree <- treeMetropolis(trees[[m]], x, yb, eta, r1, r0, priorMub$a_beta, priorMub$b_beta, priorMu1$variance, priorMu0$variance)
    tree <- updateMus(tree, x, yb, r1, r0, eta, priorMu1, priorMub, priorMu0)
    trees[[m]] <- tree

    trees
  }

  sampleZ <- function(y,y.hat){
    abs(rnorm(length(y),y.hat,1))*(-1)^{y==0}
  }
  # Burn
  #sigma2 <- 1
  baseTree <- list(variable = NA, rule = NA, expectedY = list(c(0,0,0)),
                   children = list(NA), parent = 0L)
  currentTrees <- vector(mode = "list", length = trees)
  for (t in seq(trees)) currentTrees[[t]] <- baseTree
  if (mcmc$burnin) {
    cat("Burning in.\n")
    checkpoints <- round(seq(0.1, 1, by = 0.1) * mcmc$burnin)
    for (iter in seq(mcmc$burnin)) {
      Predict1 <- rowSums(do.call(cbind, lapply(currentTrees, predictTree, x = x, component=1)))
      z1 <- sampleZ(y$y1,Predict1)
      Predict0 <- rowSums(do.call(cbind, lapply(currentTrees, predictTree, x = x, component=3)))
      z0 <- sampleZ(y$y0,Predict0)
      for (t in seq(trees)) {
        # cat("\rUpdating tree ", t, "\t", sep = "")
        currentTrees <- updateOneTree(currentTrees, x, t,  z1, z0, y$yb, priorMu1, priorMub, priorMu0)
      }
      if (iter %in% checkpoints)
        cat("\rCompleted ", round(iter / mcmc$burnin * 100, 1), "%.\t",
            sep = "")
    }
    cat("\rBurnin finished.\t\n")
  }


  cat("Sampling now.\n")
  checkpoints <- round(seq(0.1, 1, by = 0.1) * mcmc$sample)

  # For storing output
  posterior <- data.frame(sigma2 = numeric(mcmc$sample))
  for (i in seq(y)) posterior[paste0("y", i)] <- numeric(mcmc$sample)

  # Let's go!
  for (iter in seq(mcmc$sample)) {
    Predict1 <- rowSums(do.call(cbind, lapply(currentTrees, predictTree, x = x, component=1)))
    z1 <- sampleZ(y$y1,Predict1)
    Predict0 <- rowSums(do.call(cbind, lapply(currentTrees, predictTree, x = x, component=3)))
    z0 <- sampleZ(y$y0,Predict0)
    for (t in seq(trees)) {
      # cat("\rUpdating tree ", t, "\t", sep = "")
      currentTrees <- updateOneTree(currentTrees, x, t,  z1, z0, y$yb, priorMu1, priorMub, priorMu0)
    }

    allPredicts1 <- rowSums(do.call(cbind, lapply(currentTrees, predictTree, x = x, component=1)))
    allPredictsb <- rowSums(do.call(cbind, lapply(currentTrees, predictTree, x = x, component=2)))
    allPredicts0 <- rowSums(do.call(cbind, lapply(currentTrees, predictTree, x = x, component=3)))
    for (i in seq(nrow(y))){
      posterior[[paste0("y1", i)]][iter] <- allPredicts1[i]
      posterior[[paste0("yb", i)]][iter] <- allPredictsb[i]
      posterior[[paste0("y0", i)]][iter] <- allPredicts0[i]
    }

    if (iter %in% checkpoints)
      cat("\rCompleted ", round(iter / mcmc$sample * 100, 1), "%.\t",
          sep = "")
  }
  cat("\rDone.\t\t\t\t\n")

  posterior
}
