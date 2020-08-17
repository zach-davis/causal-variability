#----------------------------------------------------------
# Unit testing
#     Test conditional sampling
#----------------------------------------------------------
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source ("../utilities/utilities.R")
source ("../utilities/joint.utilities.R")
source ("../utilities/cgm.R")
source ("../models/mutsampler.R")

# CONDITIONAL SAMPLING TESTS

# Define a five variable "student" network (from Pearl and Figure 2 of Davis & Rehder, 2020).
var.names = c ('a','b','c','d','e')
ms = matrix (0, nrow=5, ncol=5, dimnames = list (var.names, var.names))
ms['a','c'] = ms['b','c'] = ms['b','d'] = ms['c','e'] = .75
bs = c (.5, 1/3, 1/3, 1/3, 1/3)
student.joint = joint.cgm.generic (ms, bs)
student.neighbors = neighbors.of.joint (student.joint)

cat ('*** Test of conditional sampling.\n')
cat ('*** Conditionalized sampling yields an answer that is closer to normative.\n')
condition = c (a=0,b=1,c=1,d=1)
cat ('Normative: p(e=1|a=0,b=1,c=1,d=1) =', v.p (student.joint, 'e', condition), '\n')

exact.student.joint = mutsampler.exact (student.joint, chainLen, student.neighbors)
cat ('Normal mutatation sampling: p(e=1|a=0,b=1,c=1,d=1) =', v.p (exact.student.joint, 'e', condition), '\n')

# Now create alternative proposal distribution that biases sampling toward states in the condition.
# This can mean that the 
q.ps = conditional.proposal.distribution (student.joint, condition, student.neighbors)
conditional.joint = mutsampler.exact (student.joint, chainLen, student.neighbors, q.ps=q.ps)
cat ('Conditionalized mutatation sampling: p(e=1|a=0,b=1,c=1,d=1) =', v.p (conditional.joint, 'e', condition), '\n')


