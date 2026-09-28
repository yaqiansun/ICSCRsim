# Visit schedule and observed interval for the non-terminal event.
#
# visit.plan(): visit times 0 < V1 < V2 < ... with gaps V(j) - V(j-1) ~ Exponential(mean 7 days),
#   continued up to the first visit after the censoring time C.
#   Only C is an input, so the schedule does not depend on T1 or T2.
visit.plan <- function(administrative.end, mean.gap=7) {
  stopifnot(length(administrative.end)==1L, is.finite(administrative.end),
            administrative.end>0, length(mean.gap)==1L,
            is.finite(mean.gap), mean.gap>0)
  visits <- numeric()
  last <- 0
  repeat {
    # Draw gaps in batches of about (C - last)/7 * 1.1 + 10 (at least 64), e.g. 2839 for C = 18000.
    # Unused gaps are discarded; another batch is drawn only if the first does not reach C.
    size <- max(64L, ceiling((administrative.end-last)/mean.gap*1.1)+10L)
    gaps <- rexp(size, rate=1/mean.gap)
    stopifnot(all(gaps>0))
    batch <- last+cumsum(gaps)
    over <- which(batch>administrative.end)
    # Return all visits up to and including the first visit after C.
    if(length(over)) return(c(visits,batch[seq_len(over[1])]))
    visits <- c(visits,batch)
    last <- tail(visits,1)
  }
}

# observe.visits(): observed (L1, R1, d1) from one schedule, onset T1 and Y2 = min(T2, C).
# Example: visits 5, 12, 20, 31, 40; T1 = 15; Y2 = 35
#   visits <= Y2 (completed): 5, 12, 20, 31
#   completed visits >= T1:   20, 31
#   R1 = 20 (first of these), L1 = 12 (the visit before) -> 12 < T1 = 15 <= 20,
#   i.e. onset lies in (L1, R1] and is first recorded at R1.
observe.visits <- function(plan, onset, followup.end) {
  stopifnot(length(plan)>0, all(is.finite(plan)), all(diff(c(0,plan))>0),
            length(onset)==1L, !is.na(onset), onset>0,
            length(followup.end)==1L, is.finite(followup.end), followup.end>0)
  completed <- plan[plan<=followup.end]
  positive <- which(completed>=onset)
  if(length(positive)) {
    # d1 = 1: R1 = first completed visit >= T1; L1 = previous visit (0 if none), so L1 < T1 <= R1 <= Y2.
    j <- positive[1]
    L <- if(j==1) 0 else completed[j-1]
    return(c(L1=L,R1=completed[j],d1=1,last.negative.visit=L,
             completed.visits=length(completed)))
  }
  # d1 = 0 (no completed visit >= T1; e.g. T1 = 33, Y2 = 35 above): L1 = Y2, R1 = Inf.
  # last.negative.visit (last completed visit, 0 if none) and completed.visits are diagnostics only.
  c(L1=followup.end,R1=Inf,d1=0,
    last.negative.visit=if(length(completed)) tail(completed,1) else 0,
    completed.visits=length(completed))
}
