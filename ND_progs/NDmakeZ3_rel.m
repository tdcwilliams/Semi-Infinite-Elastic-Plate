function [Rts,rts,H,el,del0,L,char]=NDmakeZ3_rel(Z,hhh,N)
%% MAIN NONDIMENSIONALISATION PROGRAM WHEN 3 SHEET THICKNESSES ARE INVOLVED.
%% Differs from NDmakeZ3.m in that it nondimensionalises with respect to
%% the thickest ice sheet, whereas NDmakeZ3.m always ND's wrt the lh ice sheet.
%%
%% CALL: [Rts,rts,H,el,del0,L,char]=NDmakeZ3_rel(Z,[h0 h1 h2],N)
%% INPUTS: Z=T or [T,th] or {h,T,H_dim} or {h,T,H_dim,th},
%% N is the number of imaginary roots wanted.
%%
%% OUTPUTS: Rts={gam1,K,gam2}, where Rts{j} is the (N+3)- or (N+1)-vector
%% of roots to the dispersion relation for the jth ice sheet;
%% el=K(1)*sin(theta);
%% rts{j}=(Rts{j}.^2-el^2) is the wavenumber in the x direction for the jth
%% region (square roots are taken from the upper complex half-plane);
%% H is the nondimensional depth, del0=[lam,-mu], L is the natural length
%% and char=[L_ice, T_ice], where L_ice is the characteristic length
%% for the lh ice and T_ice is its characteristic time.

if length(Z)==1
  T      = Z;
  Z      = NDphyspram([]);
  Z{2}   = T;
elseif length(Z)==2
  T      =Z(1);
  theta  = Z(2);
  Z      = NDphyspram([]);
  Z{2}   = T;
  Z{4}   = theta;
elseif length(Z)==3
  T      = Z{1};
  th     = Z{2};
  H_dim  = Z{3};
  Z      = {1,T,H_dim,th};
end

[hhh2,jj]   = sort(hhh);
hmax        = hhh(jj(3));
Z{1}        = hmax;
if hhh(jj(2))==hhh(jj(1))
  hr     = hhh(jj(1))/hmax;
  y      = ND2(Z,hr,N);
  H      = y{3};
  th     = pi*y{4}/180;
  del0   = y{5};
  char   = y{6};
  L      = y{7};
  %%
  Rts{jj(3)}   = y{1};
  Rts{jj(1)}   = y{2};
  Rts{jj(2)}   = y{2};
elseif hhh(jj(2))==hhh(jj(3))
  hr     = hhh(jj(1))/hmax;
  y      = ND2(Z,hr,N);
  H      = y{3};
  th     = pi*y{4}/180;
  del0   = y{5};
  char   = y{6};
  L      = y{7};
  %%
  Rts{jj(3)}   = y{1};
  Rts{jj(1)}   = y{2};
  Rts{jj(2)}   = y{1};
else
  hr     = [hhh(jj(1)) hhh(jj(2))]/hmax;
  y      = ND3(Z,hr,N);
  Gam    = y{2};%={K,{gam1,gam2},H,th,del0,char,L}
  H      = y{3};
  th     = pi*y{4}/180;
  del0   = y{5};
  char   = y{6};
  L      = y{7};
  %%
  Rts{jj(3)}   = y{1};
  Rts{jj(1)}   = Gam{1};
  Rts{jj(2)}   = Gam{2};
end

gam   = Rts{1};
el    = gam(1)*sin(th);
for j=1:3
  gam    = Rts{j};
  alp    = sqrt(gam.^2-el^2);
  rts{j} = alp.*sign( (imag(alp)>=0)-.5 );
end
