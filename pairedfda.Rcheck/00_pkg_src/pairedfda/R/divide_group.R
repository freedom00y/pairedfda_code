divide_group<-function(n,K)
{
  num_per_group = floor(n/K)
  no_group1 = n-K*num_per_group
  #no_group0 = K-no_group1
  set.seed(123*n)
  group_ind = sample(1:K,no_group1)
  no_group_member = rep(num_per_group,K)
  no_group_member[group_ind] = num_per_group+1
  
  out = rep(0,n)
  ind = cumsum(c(0,no_group_member))
  ord = sample(n)
  for(i in 1:K)
  {
    ind1 = ind[i]+1
    ind2 = ind[i+1]
    if(ind1>ind2)
      next
    out[ord[ind1:ind2]] = i
  }
  return(list(group_label = out,
              num_in_group = no_group_member))
}