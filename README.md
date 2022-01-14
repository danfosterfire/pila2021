
# Todos

## Reorganize to better track workflow plan

## Implement posterior retrodictive checks

## Implement posterior predictive checks

## Divy dataset into training and validation piles

Done

## Switch from midpoint to class median tree sizes for IPM (done)

I think this is more biologically realistic for small size classes - 
treating all the new recruits as if they're 2.5" dbh is bugging me. Many of 
the observed new recruits are seedlings much smaller than 2.5" dbh. This is 
in line with the findings/recommendations from the Doak paper on IPMs/MPMs, 
which suggests using class median sizes rather than midpoints. I believe 
my case is an exception to their caution against using the median rule when 
different size classes are sampled at different intensity: Because my class 
boundaries have been carefully constructed such that sampling intensity 
is consisting *within a class*, the within-class sample medians are still 
realistic pictures of the population size distribution within a class. This 
does get tricky for larger size classes though - the macro plot DBH cutoff 
varies from plot to plot, so some large trees were sampled at higher intensity 
than others. I don't think this is likely to affect the results much - theres 
only a couple of size classes where there is ambiguity about the sampling 
intensity (24-36 in), and I don't think lumping trees of different intensities 
together is going to throw off the class median all that much. Will probably 
need to do median sizes of all species to make sure we have a decent sample 
in each size class. 

## Consider log size 

Maybe binning sizes based on a log scale, or using log-size in the survival
and or fecundity models. 