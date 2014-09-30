function alpha=GEN_sqrt(alpha_sq)

alpha=sqrt(alpha_sq); alpha=sign((imag(alpha)>=0)-.5).*alpha;
