function y=GEN_howmanyterms(decay_exp,tol)

y=round(exp( -log(tol)/decay_exp ));