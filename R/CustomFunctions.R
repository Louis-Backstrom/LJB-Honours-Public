hoover <- function (x, distribution = NULL) {
  if (is.null(distribution)) {
    return (0.5 * sum(abs(x - mean(x))) / sum(x))
  }
  
  if (length(x) != length(distribution)) {
    stop ("Vector lengths are not equal")
  }
  
  if (isTRUE(all.equal(sum(x), sum(distribution))) == FALSE) {
    new_x <- x * sum(distribution) / sum(x)
    warning ("Vector sums are not equal. Applying transformation to make them equal: this may have unexpected results!")
    return (c(0.5 * sum(abs(new_x - distribution)) / sum(new_x), sum(distribution) / sum(x)))
  }
  
  return (0.5 * sum(abs(x - distribution)) / sum(x))
}

sign_chr <- function (x) {
  if (x < 0) {
    return ("-")
  }
  
  if (x == 0) {
    return ("")
  }
  
  if (x > 0) {
    return("+")
  }
}

flexi_round <- function (x, to) {
  return (round(x / to) * to)
}

simulate_hoover <- function (groups, observations) {
  observations <- sample(groups, observations, replace = TRUE) %>% 
    tibble("values" = .) %>% 
    group_by(values) %>% 
    summarise(observations = n()) %>% 
    left_join(tibble("values" = 1:groups), ., by = "values") %>% 
    mutate(observations = replace_na(observations, 0)) %>% 
    pull(observations)
  
  return (hoover(observations))
}

occ_change <- function (start, change, length.out) {
  if (change < 0) {
    x <- start * (1 + change) ^ c(0:(length.out - 1))
  }
  
  if (change == 0) {
    x <- rep(start, length.out)
  }
  
  if (change > 0) {
    x <- 1 - ((1 - start) * (1 - change) ^ c(0:(length.out - 1)))
  }
  
  return (x)
}

rmse <- function (errors) {
  x <- sqrt(mean(errors ^ 2))
  
  return (x)
}

rss <- function (errors) {
  x <- sum(errors ^ 2)
  
  return (x)
}

overlaps <- function (x, y) {
  overlaps <- FALSE
  
  if (min(x) >= min(y) && min(x) <= max(y)) {
    overlaps <- TRUE
  }
  
  if (max(x) >= min(y) && max(x) <= max(y)) {
    overlaps <- TRUE
  }
  
  if (min(y) >= min(x) && min(y) <= max(x)) {
    overlaps <- TRUE
  }
  
  if (max(y) >= min(x) && max(y) <= max(x)) {
    overlaps <- TRUE
  }
  
  return(overlaps)
}

