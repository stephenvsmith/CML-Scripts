############ Data Generation ############

# This code was given by a colleague

generate.data.grid <- function(data.grid = data.frame(network = "child",
                                                      data.type = "continuous",
                                                      n.dat = 10,
                                                      n.obs = 1000,
                                                      c.ratio = 0,
                                                      max.in.degree = Inf,
                                                      lb = 0.5,  # lower bound of coefficients
                                                      ub = 1,  # upper bound of coefficients
                                                      low = 0,  # lower bound of variances if continuous, of number of levels if discrete
                                                      high = 1,  # upper bound of variances if continuous, of number of levels if discrete
                                                      scale = TRUE,
                                                      stringsAsFactors = FALSE),
                               out.dir,
                               array_num,
                               path.start='~/',
                               verbose = TRUE,
                               var_list=NULL){

  ## initialize output directory
  dir.check(out.dir)
  if (verbose) cat('generate.data.grid: Writing data.grid to output directory \n')
  write.table(data.grid, paste0(out.dir, '/data.grid.txt'))  # write data grid in out.dir

  for (i in 1:nrow(data.grid)){

    data.row <- data.grid[i,, drop = FALSE]
    data.dir <- paste0(out.dir, '/', data.row$network,
                       '_',array_num)
    # data.dir <- paste0(out.dir, '/', data.row$network,
    #                    '; n = ', data.row$n.obs,
    #                    '; c = ', data.row$c.ratio)
    dir.check(data.dir)
    write.txt <- function(obj, txt, col.names = FALSE, row.names = FALSE){  # function to write obj to data.dir
      write.table(obj, paste0(data.dir, sprintf('/%s.txt', txt)),
                  col.names = col.names, row.names = row.names)
      # if (verbose) cat(sprintf('generate.data.grid: Writing %s to data directory \n', txt))
    }
    write.txt(data.row, 'data.row', col.names = TRUE, row.names = TRUE)  # write data row in data.dir
    if (verbose) cat(sprintf('generate.data.grid: Generating data for network %s \n', data.row$networks))

    networks <- trimws(strsplit(data.row$network, split = ',')[[1]])
    true <- bnlearn2edgeList(bn.net = networks,
                             c.ratio = data.row$c.ratio,
                             max.parents = data.row$max.in.degree,
                             seed = i,
                             verbose = verbose,
                             net.dir = rds_dir,
                             path.start=path.start,var_list=var_list)
    true_edgeList <- true$edgeList
    trueDAG <- as.matrix(sparsebnUtils::get.adjacency.matrix(true_edgeList))
    trueCPDAG <- pcalg::dag2cpdag(trueDAG)
    trueGroup <- true$group
    n.edge <- sparsebnUtils::num.edges(true_edgeList)
    n.node <- sparsebnUtils::num.nodes(true_edgeList)

    ## generate order matrix
    set.seed(i)
    order.mat <- replicate(data.row$n.dat, sample(1:n.node))

    ## write
    write.txt(trueDAG, 'trueDAG')
    write.txt(trueCPDAG, 'trueCPDAG')
    write.txt(trueGroup, 'trueGroup')
    write.txt(order.mat, 'order.mat')

    ## generate data
    set.seed(i)
    if (data.row$data.type == 'continuous'){

      params.mat <- replicate(data.row$n.dat, {
        coefs <- runif(n.edge, min = data.row$lb, max = data.row$ub) * sample(c(-1, 1), n.edge, replace = TRUE)  # sample correlation from Unif([-ub, -lb] U [lb, ub])
        vars <- runif(n.node, min = data.row$low, max = data.row$high)  # sample noise variance
        params <- c(coefs, vars)  # parameter format for generating data
        return(params)
      })
      write.txt(params.mat, 'params.mat')

      for (j in 1:data.row$n.dat){
        if(verbose) cat(sprintf('generate.data.grid: Generating dataset %s \n', j))
        set.seed(j)

        vars <- utils::tail(params.mat[,j], n.node)  # extract variances
        noise.matrix <- sapply(vars, function(x){rnorm(data.row$n.obs, 0, sqrt(x))})
        dat <- generate_dat(true_edgeList, params.mat[,j], noise.matrix, data.row$n.obs, if_scale = data.row$scale)
        write.txt(dat, paste0('data', j))
      }  # end for
    } else if (data.row$data.type == 'discrete'){

      if (any(is.na(c(data.row$low, data.row$high))) || any(c(data.row$low, data.row$high) < 0)){  # use true levels
        cat(sprintf('generate.data.grid: Using true number of levels ranging from %s to %s \n',
                    min(true$n.levels), max(true$n.levels)))
        N_LEVELS <- lapply(1:data.row$n.dat, function(x) true$n.levels)
      } else{
        low <- data.grid[i,]$low <- data.row$low <- max(floor(data.row$low), 2)
        high <- data.grid[i,]$high <- data.row$high <- max(floor(data.row$high), floor(data.row$low))
        cat(sprintf('generate.data.grid: Generating number of levels ranging from %s to %s \n',
                    low, high))
        N_LEVELS <- lapply(1:data.row$n.dat, function(x){
          sample(rep(low:high, 2), size = n.node, replace = TRUE)
        })
      }

      params.list <- lapply(1:data.row$n.dat, function(x){
        coefs <- discretecdAlgorithm::coef_gen(true_edgeList, n_levels = N_LEVELS[[x]],
                                               FUN = function(n) runif(n, min=data.row$lb, max=data.row$ub))
        return(coefs)
      })
      saveRDS(params.list, paste0(data.dir, '/params.list.rds'))

      for (j in 1:data.row$n.dat){

        set.seed(j)
        while (T){
          dat <- discretecdAlgorithm::generate_discrete_data(graph = true_edgeList, params = params.list[[j]], n = data.row$n.obs, n_levels = N_LEVELS[[j]])
          if (verbose) cat(sprintf('generate.data.grid: Dataset %s has %s/%s variables with only one level \n',
                                   j, sum(sapply(data.frame(dat), sd) == 0), n.node))
          if (sum(sapply(data.frame(dat), sd) == 0) == 0){
            break
          }
        }
        write.txt(dat, paste0('data', j))
      }  # end for
    }  # end discrete
  }
}


#	Function: bn.repository
#
#	Description:
#
#		Loads a bnlearn Bayesian network repository.
#
#	Arguments:
#
#		bn.net - Name of network from bnlear Bayesian network repository.
#   net.dir - Bayesian network directory containing folder rds containing bnlearn networks stored as network.rds files.
#
#	Returns:
#
#		out - Parent edge weight list. Each named element contains parents, indices, and weights.
#
#	Examples:
#
#

bn.repository <- function(bn.net,  # character scalar of
                          path.start='~/',
                          net.dir = NULL){  # directory of bn.net .rda files
  if (is.null(net.dir)){
    net.dir <- paste0(path.start,"Desktop/Research/projects/bn_data_generation/networks")
  }
  # load(sprintf(paste0(net.dir, "/rda/%s.rda"), tolower(bn.net)))  # load bn.fit object from rda file
  bn <- readRDS(sprintf(paste0(net.dir, "/%s.rds"), tolower(bn.net)))

  return(bn)
}


#	Function: bnlearn2edgeList
#
#	Description:
#
#		Creates edgeList (sparsebnUtils) object, by connecting networks if desired, from the bnlearn Bayesian network repository.
#
#	Arguments:
#
#		bn.net - Character object with name(s) of networks in bnlearn repository. Repeat network name k times for k identically structured subgraphs.
#   c.ratio - Connectivity ratio constant c. Subgraphs are connected with (c * total edges within subgraphs) edges.
#   max.parents - Maximum number of parents per node (maximum in-degree). If a node has max.parents or more parents, an edge to the node will not be considered.
#   net.dir - Bayesian network directory containing folder rds containing bnlearn networks stored as network.rds files.
#   seed - Random seed.
#
#	Returns:
#
#		List containing:
#
#       edgeList - Fused sparsebn edgeList object.
#       group - True group of each node.
#       n.edges - Numeric vector containing number of edges within subgraphs (s_w) and number of edges between subgraphs (s_b).
#       n.parents - Number of parents for each node.
#       seed - Random seed argument.
#
#	Examples:
#
#   > pigs <- bnlearn2edgeList(bn.net = "pigs")
#   > plot(pigs$edgeList)
#		> mix <- bnlearn2edgeList(bn.et = c("pigs", "pathfinder", "hailfinder"),
#                            c.ratio = 0.1)
#   > plot(mix$edgeList)
#

bnlearn2edgeList <- function(bn.net,
                             c.ratio = 0,
                             max.parents = Inf,
                             net.dir = NULL,
                             seed = NULL,
                             verbose = FALSE,
                             path.start='~/',
                             var_list=NULL){

  if (is.null(seed)){
    seed <- 1
  }
  set.seed(seed)

  if (!is.null(var_list)) net.dir <- var_list$rds_dir

  bn.net.list <- lapply(bn.net, function(x){
    bn <- readRDS(sprintf(paste0(net.dir, "/%s.rds"), tolower(x)))
    # load(sprintf(paste0(net.dir, "/rda/%s.rda"), tolower(x)))  # load bn.fit object from rda file
    sparsebnUtils::edgeList(lapply(bn, function(x) match(x$parents, names(bn))))  # get edge list
  })

  n.net <- length(bn.net)
  net.len <- sapply(bn.net.list, length)
  net.edgeList <- do.call("c",
                          lapply(1:n.net, function(x)
                            lapply(bn.net.list[[x]], function(y)
                              y + sum(net.len[1:x]) - net.len[x]))
  )  # end do.call
  names(net.edgeList) <- paste0("V", 1:length(net.edgeList))
  # names(net.edgeList) <- 1:length(net.edgeList)

  ## add random edges between
  if(n.net > 1 &&  # no c if only 1 network
     c.ratio){
    n.btwn <- ceiling(c.ratio * length(unlist(net.edgeList)))  # pick c * k * s_sub edges
    cum.net.len <- c(0, cumsum(net.len))

    if (verbose) cat(sprintf('bnlearn2edgeList: Connecting subgraphs with %s edges \n', n.btwn))
    continue <- TRUE

    ## repeat until success
    while (continue){

      ## regenerate edge list
      net.edgeList <- do.call(c,
                              lapply(1:n.net, function(x)
                                lapply(bn.net.list[[x]], function(y)
                                  y + sum(net.len[1:x]) - net.len[x]))
      )  # end do.call
      names(net.edgeList) <- paste0("V", 1:length(net.edgeList))
      # names(net.edgeList) <- 1:length(net.edgeList)
      net.adjMat <- as.matrix(sparsebnUtils::get.adjacency.matrix(sparsebnUtils::edgeList(net.edgeList)))  # store as adjacency matrix
      valid.child <- which(colSums(net.adjMat) < max.parents)  # valid children to search have fewer than maximum parents

      ## generate all possible edges between
      btwn.edges <- do.call(rbind, lapply(1:n.net, function(x){  # all possible combinations
        expand.grid(cum.net.len[x] + 1:net.len[x],  # from group x
                    setdiff(1:length(net.edgeList), cum.net.len[x] + 1:net.len[x]))  # to a different group
      }))  # includes both directions
      btwn.edges <- btwn.edges[btwn.edges[,2] %in% valid.child,]  # choose all rows that have child that is valid
      btwn.edges <- btwn.edges[sample(1:nrow(btwn.edges)),]  # reorder btwn edges to be chosen in random order

      ## add edges one at a time
      for (i in 1:nrow(btwn.edges)){  # for each candidate
        net.adjMat[btwn.edges[i,1], btwn.edges[i,2]] <- 1  # parents in rows, child in column

        if (! gRbase::is.DAG(net.adjMat) |  # if not a DAG
            sum(net.adjMat[, btwn.edges[i,2]]) > max.parents){  # or if child has too many parents
          net.adjMat[btwn.edges[i,1], btwn.edges[i,2]] <- 0  # remove edge
          if (i == nrow(btwn.edges)){  # if not a DAG and exhausted candidates
            cat(sprintf('bnlearn2edgeList: Failed attempt with %s edges added \n', sum(net.adjMat) - length(unlist(net.edgeList))))
          }
        } else if (sum(net.adjMat) == (length(unlist(net.edgeList)) + n.btwn)){  # if is DAG and all edges acquired
          if (verbose) cat('bnlearn2edgeList: Successfully connected subgraphs \n')
          continue <- FALSE
          break  # no longer search through edges
        }

        if (verbose > 1 && (i %% 100) == 0){
          cat(sprintf('bnlearn2edgeList: Successfully added %s edges \n', i))
        }
      }
      net.edgeList <- sparsebnUtils::to_edgeList(igraph::graph_from_adjacency_matrix(net.adjMat))
    }  # end while
  } else n.btwn <- 0

  # groups
  group <- rep(1:n.net, net.len)

  # get number of original levels for each node
  bn.list <- lapply(bn.net, function(x){
    bn <- readRDS(sprintf(paste0(net.dir, "/%s.rds"), tolower(x)))
    # load(sprintf(paste0(net.dir, "/rda/%s.rda"), tolower(x)))  # load bn.fit object from rda file
    return(bn)
  })
  n.levels <- unlist(lapply(bn.list, function(x){
    sapply(x, function(y){
      dim(y$prob)[1]
    })
  }))
  #names(n.levels) <- names(net.edgeList)

  # get number of parents per node
  n.parents <- sapply(net.edgeList, length)

  out <- list(edgeList = sparsebnUtils::edgeList(net.edgeList),
              group = group,
              n.edges = c(s_w = length(unlist(net.edgeList)) - n.btwn, s_b = n.btwn),
              n.levels = n.levels,
              n.parents = n.parents,
              seed = seed)

  return(out)
}


#	Function: sparse_to_edgeWeightList
#
#	Description:
#
#		Converts a sparse object to an parent edge weight list.
#
#	Arguments:
#
#		x - sparsebnUtils sparse object.
#   nodes - Names of the nodes.
#
#	Returns:
#
#		out - Parent edge weight list. Each named element contains parents, indices, and weights.
#
#	Examples:
#
#

sparse_to_edgeWeightList <- function(x, nodes) {
  stopifnot(sparsebnUtils::is.sparse((x)))
  # sp <- sparsebnUtils::as.sparse(x) # NOTE: no longer a bottleneck under sparsebnUtils v0.0.4

  # nodes <- colnames(x)
  stopifnot(x$dim[1] == x$dim[2])

  out <- lapply(vector("list", length = x$dim[1]), function(z) list(parents = character(0), index = integer(0), weights = numeric(0)))
  names(out) <- nodes
  for(j in seq_along(x$cols)){
    child <- x$cols[[j]]
    parent <- x$rows[[j]]
    weight <- x$vals[[j]]
    parents <- c(out[[child]]$parents, nodes[parent]) # !!! THIS IS SLOW
    index <- c(out[[child]]$index, parent) # !!! THIS IS SLOW
    weights <- c(out[[child]]$weights, weight) # !!! THIS IS SLOW
    out[[nodes[child]]] <- list(parents = parents, index = index, weights = weights)
  }

  return(out)
}


#	Function: generate_dat
#
#	Description:
#
#		Generates data.
#
#	Arguments:
#
#		true_edgeList - True parent edge list.
#   params - Parameters containing coefficients for each edge and noise variance of each variable.
#   noise_matrix - Pre-generated noise matrix using the noise variance of each variable in params.
#   n_obs - Number of observations to generate.
#   if_scale - Logical value indicating whether or not to scale each column to have variance = 1 in the generating process.
#
#	Returns:
#
#		data_matrix - Generated data matrix.
#
#	Examples:
#
#

generate_dat <- function(true_edgeList, params, noise_matrix, n_obs, if_scale = FALSE) {
  node_name <- names(true_edgeList)
  true_graph <- sparsebnUtils::to_igraph(true_edgeList)
  if (is.null(node_name)) {
    ts <- igraph::topo_sort(true_graph)
  } else {
    ts <- match(names(igraph::topo_sort(true_graph)), node_name)
  }
  n_node <- length(true_edgeList)

  vars <- utils::tail(params, n_node) # parameters associated with variances
  names(vars) <- node_name
  coefs <- params[1:(length(params) - n_node)] # parameters associated with edge weights
  sp <- sparsebnUtils::as.sparse(true_edgeList)
  sp$vals <- coefs # previous line leaves NAs for values in sparse object; need to fill these in
  edgelist <- sparse_to_edgeWeightList(sp, node_name)

  my_adj <- matrix(0,nrow = n_node,ncol = n_node)
  for (i in 1:length(edgelist)){
    ll <- edgelist[[i]]
    ind <- ll$index
    ws <- ll$weights
    my_adj[ind,i] <- ws
  }
  
  # noise_matrix <- sapply(vars, function(x){rnorm(n_obs, 0, sqrt(x))})

  data_matrix <- matrix(0, nrow = n_obs, ncol = n_node)
  for(i in 1:n_node) {
    child <- ts[i]
    parent <- true_edgeList[[child]]
    n_parent <- length(parent)
    if (n_parent) {
      beta <- matrix(edgelist[[child]]$weights, ncol = 1)
      new_data <- data_matrix[, parent, drop = FALSE] %*% beta + noise_matrix[, child, drop = FALSE]
    } else {
      new_data <- noise_matrix[, child, drop = FALSE]
    }
    if (if_scale == TRUE) {
      new_data <- scale(new_data)
    }
    data_matrix[, child] <- new_data
  }
  return(data_matrix)
}



#	Function: dir.check
#
#	Description:
#
#		Checks if a directory exists, creating the directory if necessary.
#
#	Arguments:
#
#		dir - Character working directory.
#   min.depth - Minimum depth. If folders earlier than min.depth do not exist (likely because of typos), stop.
#
#	Returns:
#
#		None. Createsmissing working directory.
#
#	Examples:
#
#

dir.check <- function(dir, min.depth = 0){

  folders <- strsplit(dir, split = '/')[[1]]
  folders <- folders[folders != '']  # end(s)

  stopifnot(length(folders) >= min.depth)  # stop if not enough folders

  for (i in 1:length(folders)){

    temp.dir <- paste0(c('', folders[1:i]), collapse = '/')
    if (i < min.depth){  # below minimum depth
      stopifnot(dir.exists(temp.dir))  # if directory does not exist, stop
    } else if (!dir.exists(temp.dir)) {  # if at least minimum depth and doesn't exist
      cat(sprintf('dir.check: Creating folder `%s`\n', folders[i]) )
      dir.create(file.path(temp.dir))  # create directory
    }
  }
}
