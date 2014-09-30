function dAn=OP_diff_legendre(An,m,Lc,do_test)

want_coeffs=~iscell(An);
if ~want_coeffs
  Ngg=An{1};
  An=ones(Ngg+1,1);
else
  Ngg=length(An)-1;
end
Nout=Ngg-m;




if nargin>=3%%allow for intervals other than [-1 1]
  dAn=dAn/Lc^m;
end

if do_test

end