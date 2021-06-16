#------------------------
# evaluating the distributions
#------------------------
# generating sample output from a model ------------------
var.names = c ('x','y','z')
ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
ms['x','y'] = .5; ms['y','z'] = .75
bs = c (.75, .5, .5)
a <- param_variability(ms = ms,
                       bs = bs,
                       ms_conc = .1,
                       bs_conc = .1,
                       nSamples = 6)

# I'm not sure how to do this systematically - we should discuss
cpj_eval <- function(ms, bs, model) {}

guess_fun <- function(dist, guess_prob) {
  dist2 = (1-p) * dist$p
  dist2[dist2$x == .5] = dist2[dist2$x == .5] + guess_prob
  return(dist2)
}
