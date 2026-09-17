# Scheme 3: schedule generation and observation mapping are separate operations.
# No event time is an input to the schedule generator.
visit.plan0914 <- function(administrative.end, mean.gap=7) {
  stopifnot(length(administrative.end)==1L, is.finite(administrative.end),
            administrative.end>0, length(mean.gap)==1L,
            is.finite(mean.gap), mean.gap>0)
  visits <- numeric()
  last <- 0
  repeat {
    # Generate strictly positive iid exponential increments, including the first
    # visit beyond C. Stopping at C does not use either generated event time.
    size <- max(64L, ceiling((administrative.end-last)/mean.gap*1.1)+10L)
    gaps <- rexp(size, rate=1/mean.gap)
    stopifnot(all(gaps>0))
    batch <- last+cumsum(gaps)
    over <- which(batch>administrative.end)
    if(length(over)) return(c(visits,batch[seq_len(over[1])]))
    visits <- c(visits,batch)
    last <- tail(visits,1)
  }
}

observe.visits0914 <- function(plan, onset, followup.end) {
  stopifnot(length(plan)>0, all(is.finite(plan)), all(diff(c(0,plan))>0),
            length(onset)==1L, !is.na(onset), onset>0,
            length(followup.end)==1L, is.finite(followup.end), followup.end>0)
  completed <- plan[plan<=followup.end]
  positive <- which(completed>=onset)
  if(length(positive)) {
    j <- positive[1]
    L <- if(j==1) 0 else completed[j-1]
    return(c(L1=L,R1=completed[j],d1=1,last.negative.visit=L,
             completed.visits=length(completed)))
  }
  # Agreed working observation convention: no detected onset -> censor at Y2.
  # This is not a new examination at Y2, and does not change latent event truth.
  c(L1=followup.end,R1=Inf,d1=0,
    last.negative.visit=if(length(completed)) tail(completed,1) else 0,
    completed.visits=length(completed))
}
