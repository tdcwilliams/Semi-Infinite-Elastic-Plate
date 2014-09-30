function [K,del] = RTS_ice_doubleroots(H,nc)
%% CALL: K = RTS_ice_doubleroots(H,nc);
%% calculates the double roots in the nc-th interval of
%% the dispersion relation 1/(K^4 + del) = K*tanh(K*H).
%%
%% DEPENDENCIES: RTS_imag_roots_ice.c (.mexglx),
%%  RTS_imag_root_ice.c (.mexglx).
do_test=0;
TOL=1e-12;

[H3,del3,k3]=RTS_ice_tripleroot(nc);
if abs(H-H3)<TOL
  K=k3*[1;1];
  del=del3*[1;1];
else
  w3=k3*H3/1i;
  [w20,del20,w21,del21]=double_roots(H,nc,{H3,del3,w3});
  K=1i/H*[w20;w21];
  del=[del20;del21];
end

function [w20,del20,w21,del21]=double_roots(H,nc,Z)

H3=Z{1};
del3=Z{2};
w3=Z{3}; %[H3,H]
%%
I1=nc*pi;
I0=I1-pi;
tol=1e-8;
I0=max([I0 tol]);
jc=2*nc-1;
Hc=( (jc*pi)^4/4 )^.2;
wc=jc*pi/2;
if H>H3%%no roots.
  w20=Inf;
  w21=Inf;
elseif H==H3%%get a triple root.
  w20=w3;
  w21=w3;
else
  %%know w21 is always in (w3,I1)
  i0=w3;
  i1=I1;
  w21=double_root(H,i0,i1);
  if H==Hc %%know w20=wc.
    w20=wc;
  elseif H<Hc %%know w20 is to the left of wc
    i0=I0;
    i1=wc;
    w20=double_root(H,i0,i1);
  else %%w20 is in (I0,w3)
    i0=wc;
    i1=w3;
    w20=double_root(H,i0,i1);
  end
end
del20=-((w20/H)^4+H*cot(w20)/w20);
del21=-((w21/H)^4+H*cot(w21)/w21);

function w=double_root(H,i0,i1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=1e-8;
w=fzero(@q_fxn,[i0, i1],optimset('TolX',tol),H);

function q=q_fxn(w,H)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c2=cos(2*w);
s2=sin(2*w);
q=0*w-2*H^5;
j=find(w~=0);
q(j)=2*w(j).^4.*(1-c2(j))-H^5*(1+s2(j)/2./w(j));