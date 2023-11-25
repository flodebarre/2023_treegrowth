# Generate genealogies and plot tree for samples

# Generic stuff ####

mai <- 0.05 + c(0, 0, 0.05, 0) # margins

# Functions ####

#
# Generate population and its history, discrete non-overlapping generations
#
genealogy <- function(mOff = 1.15, vOff = NA, distribOff = "Poisson", nGen = 10, n0 = 10, nmax = 1000){
  # Input parameters
  # mOff  Average number of offspring
  # vOff  Variance of the number of offspring (not needed for Poisson)
  # distribOff  "Poisson"
  # nGen  Number of generations
  # n0  Initial population size
  # nmax  Max population size, has to be multiple of 10
  
  if(nmax%%10 != 0) stop("nmax has to be a multiple of 10")
  nPad <- log(nmax, 10) # Number of zeros for padding the result
  # Function to format the output
  formatID <- function(v){
    formatC(v, width = nPad, format = "d", flag = "0")
  }

    # Function to draw the numbers of offspring
  drawOffNb <- function(n, mOff, vOff){
    if(distribOff == "Poisson"){
      out <- rpois(n, mOff)
    }
    
    if(distribOff == "Constant"){
      # Draw parents
      p <- sample(1:n, size = n, replace = TRUE)
      # Count values
      off1 <- table(p)
      tb1 <- data.frame(id = formatID(as.numeric(names(off1))), noff = c(off1))
      tb2 <- data.frame(id = formatID(1:n))
      tb <- merge(tb1, tb2, all = TRUE)
      tb[is.na(tb$noff), "noff"] <- 0
      out <- c(tb$noff)
    }
    out
  } 
  
  
  # How to store the results? 
  # Let's try a list of vectors
  pop <- list() # Initialize output
  pop[[1]] <- formatID(1:n0)
  
  for(igen in 2:nGen){
    # Draw numbers of offspring for each individual of the igen-1 generation
    nAdu <- length(pop[[igen-1]]) # Number of adults
    noff <- drawOffNb(n = nAdu, mOff = mOff, vOff = vOff)
    ntot <- sum(noff) # total number of offspring
    
    # If beyond max pop size, stop and return population vector
    if(ntot > nmax){
      return(pop)
    }
    # Build output vector
    v <- c() # Initialize with empty vector
    for(i in (1:nAdu)[noff != 0]){ # Do it only for the adults with offspring
      # Define offspring IDs
      ids <- formatID(1:noff[i])
      # Aggregate offspring of each adult
      v <- c(v, paste(pop[[igen-1]][i], ids, sep = "."))
    }
    # Save it
    pop[[igen]] <- v
  }
  
  pop
}

#
# Plot the whole population history
#
plotGenealogy <- function(pop, maincol = gray(0.85), dx = 0.2, dy = 0.1, maxSize = NA){
  # pop  population list, as output of the genealogy function
  # maincol  main color
  
  nGen <- length(pop) # Number of generations
  if(is.na(maxSize)){
    maxSize <- max(unlist(lapply(pop, length))) # Max population size
  }

  # Initialize plot (empty)
  par(mai = mai)
  plot(0, xlim = dx * c(1, nGen), ylim = dy * c(-(maxSize/2), (maxSize/2)), 
       axes = FALSE, type = "n", 
       xlab = "", ylab = "")
  
  # Draw individuals and save their positions
  pos <- list() 
  for(iGen in 1:nGen){
    nn <- length(pop[[iGen]])
    pos[[iGen]] <- seq(-(nn-1)/2, (nn-1)/2)
    points(dx * rep(iGen, nn), dy * seq(-(nn-1)/2, (nn-1)/2), 
           pch = 1, col = maincol)
  }
  
  # Get padding length
  nPad <- nchar(pop[[1]][1])
  
  # Draw parent-offspring links
  for(iGen in 2:nGen){
    adults <- pop[[iGen-1]] # IDs of the previous generation individuals
    offspring <- pop[[iGen]] # IDs of the current generation individuals
    parents <- substr(offspring, 1, nchar(adults[1])) # IDs of the parents
    positionOffspring <- 1:length(offspring) # Positions of the offspring
    positionParents <- match(parents, adults) # Positions of their parents
    # Draw parent-offspring links
    segments(x0 = dx * (iGen - 1), x1 = dx * iGen, 
             y0 = dy * pos[[iGen - 1]][positionParents], 
             y1 = dy * pos[[iGen]][positionOffspring], col = maincol)
  }
  pos
}

#
# Draw a sample of the final population
#
getSample <- function(size, pop){
  # size  number of leaves to draw
  # pop  population list
  sort(sample(pop[[length(pop)]], size = 8, replace = FALSE))
}

#
# Plot the tree corresponding to the sample, on the whole pop's tree
#
plotTree <- function(sampleLeaves, pop, pos, colCoal = "red", colLink = "black", colLeaf = "orange", dx = 0.2, dy = 0.1){
  # sample  sample leaves
  # pop  population list
  # pos  positions list
  
  # Initialize
  iGen <- length(pop)
  # Current number of leaves
  nLeaves <- length(sampleLeaves)
  # Show the end leaves
  points(dx * rep(iGen, nLeaves), dy * pos[[iGen]][match(sampleLeaves, pop[[iGen]])], pch = 16, col = colLeaf)
  
  currentLeaves <- sampleLeaves
  
  for(iGen in seq(length(pop), 2, -1)){
    adults <- pop[[iGen - 1]]
    youngs <- pop[[iGen]]
    parentsLeaves <- substr(currentLeaves, 1, nchar(adults[1]))
    # Get positions
    positionLeaves <- match(currentLeaves, youngs) # Positions of the offspring
    positionParents <- match(parentsLeaves, adults) # Positions of their parents
    # Draw parent-offspring links
    segments(x0 = dx * (iGen - 1), x1 = dx * iGen, 
             y0 = dy * pos[[iGen - 1]][positionParents], 
             y1 = dy * pos[[iGen]][positionLeaves], col = colLink)
    # Identify coalescence,  if any
    if(any(duplicated(parentsLeaves))){
      # Plot points
      coalNodes <- parentsLeaves[which(duplicated(parentsLeaves))]
      points(rep(dx * (iGen -1), length(coalNodes)), 
             dy * pos[[iGen - 1]][match(coalNodes, adults)], pch = 16, col = colCoal)
    }
    # Update leaves
    currentLeaves <- unique(parentsLeaves)
    # If we have reached the MRCA, stop
    if(length(currentLeaves) == 1){
      return(0)
    }
  }
}

#
# Plot tree only for the sample
#
plotTree2 <- function(sampleLeaves, pop, colLink = "black", dx = 0.2, dy = 0.1){
  # sample  sample leaves
  # pop  population list
  
  # Here we do a simple, regular tree
  
  # Initialize plot
  nGen <- length(pop)
  par(mai = mai)
  maxSize <- length(sampleLeaves) 
  plot(0, xlim = dx * c(1, nGen), ylim = dy * c(1, length(sampleLeaves)), 
       axes = FALSE, type = "n", 
       xlab = "", ylab = "")
  
  currentLeaves <- sampleLeaves
  timeLeaves <- rep(nGen, length(currentLeaves))
  posLeaves <- seq_along(sampleLeaves)
    
  for(iGen in seq(length(pop), 2, -1)){
    adults <- pop[[iGen - 1]]
    youngs <- pop[[iGen]]
    parentsLeaves <- substr(currentLeaves, 1, nchar(adults[1]))
    
    # Identify coalescence, if any
    if(any(duplicated(parentsLeaves))){
      # ID of the parents of coalescing nodes
      coalNodes <- unique(parentsLeaves[which(duplicated(parentsLeaves))])

      newPos <- posLeaves
      newTimes <- timeLeaves
      
      for(coal in coalNodes){ # For each coalescence (usually just 1)
        # Identify the indices of the coalescing nodes
        icoal <- which(parentsLeaves == coal)
        # Their position is the average of the offspring
        newPos[icoal] <- mean(posLeaves[icoal])
        # Change their time
        newTimes[icoal] <- iGen - 1
        
        # Draw the tree
        #   lines to leaves
        segments(x0 = dx * newTimes[icoal], x1 = dx * timeLeaves[icoal], 
                 y0 = dy * posLeaves[icoal], y1 = dy * posLeaves[icoal])
        #   jonction
        segments(x0 = dx * min(newTimes[icoal]), x1 = dx * max(newTimes[icoal]), 
                 y0 = dy * min(posLeaves[icoal]), y1 = dy * max(posLeaves[icoal]))
      }
      # Update vectors
      iNoDupl <- which(!duplicated(parentsLeaves))
      posLeaves <- newPos[iNoDupl]
      timeLeaves <- newTimes[iNoDupl]
      currentLeaves <- parentsLeaves[iNoDupl]
    }else{
      currentLeaves <- parentsLeaves
    }
    
    # If we have reached the MRCA, stop
    if(length(currentLeaves) == 1){
      return(0)
    }
  }
  # If not MRCA yet, draw segments
  if(length(currentLeaves > 1)){
    segments(x0 = dx * rep(1, length(posLeaves)), x1 = dx * timeLeaves, 
             y0 = dy * posLeaves, y1 = dy * posLeaves)
  }
}

#
# Draw populations and plot them, 
# comparing exponential growth and constant population size
#
compareGrowth <- function(nGens = 200, nSamp = 8, n0 = 10, mOff = 1.2, export = TRUE, suffix = ""){
  # nGens  number of generations to simulate
  # nSamp  number of samples to draw
  # n0  initial population size (exp growth)
  # mOff  average number of offspring
  # export  whether to export the figures as png
  # suffic  suffix to add to the figure names

  # Export parameters
  wpng <- 15 # pic width
  hpng1 <- 15 # pic height (1; whole pop)
  hpng2 <- 4 # pic height (2; tree)
  res <- 200 # resolution
  units <- "cm" # unit of pic size
  
  # Exponential growth
  #  Generate population
  pop.exp <- genealogy(mOff = 1.2, nGen = nGens, n0 = 10, distribOff = "Poisson")
  #  Draw sample
  sampleLeaves.exp <- getSample(nSamp, pop.exp)
  #  Plot population history and sample
  if(export) png(paste0("../figs/fullpop_exponential", suffix, ".png"), width = wpng, height = hpng1, res = res, units = units)
  pos.exp <- plotGenealogy(pop.exp)
  plotTree(sampleLeaves.exp, pop.exp, pos.exp)
  if(export) dev.off()
  #  Plot regular tree for the sample
  if(export) png(paste0("../figs/tree_exponential", suffix, ".png"), width = wpng, height = hpng2, res = res, units = units)
  plotTree2(sampleLeaves.exp, pop.exp)
  if(export) dev.off()
  
  # Constant population size
  
  # We use the average number of individuals per generation from the previous simulation
  # (Same code as previously; not repeating comments)
  pop.cst <- genealogy(mOff = NA, nGen = nGens, n0 = mean(unlist(lapply(pop.exp, length))), distribOff = "Constant")
  sampleLeaves.cst <- getSample(nSamp, pop.cst)
  #
  if(export) png(paste0("../figs/fullpop_constant", suffix, ".png"), width = wpng, height = hpng1, res = res, units = units)
  pos.cst <- plotGenealogy(pop.cst, maxSize = max(unlist(lapply(pop.exp, length))))
  plotTree(sampleLeaves.cst, pop.cst, pos.cst)
  if(export) dev.off()
  #
  if(export) png(paste0("../figs/tree_constant", suffix, ".png"), width = wpng, height = hpng2, res = res, units = units)
  plotTree2(sampleLeaves.cst, pop.cst)
  if(export) dev.off()
  
}

# Run ####

compareGrowth(nGens = 200, nSamp = 8, n0 = 10, mOff = 1.2, export = TRUE, suffix = "")

