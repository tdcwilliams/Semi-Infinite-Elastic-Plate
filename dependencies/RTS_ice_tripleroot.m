function [H,del,K]=RTS_ice_tripleroot(n)
%% triple_roots: returns the triple roots of the dispersion
%%  relation in the nth interval. If
%%   q(w)=2*w^4*(1-c2)-H^5*(1+s2/2/w), where c2=cos(2*w)
%%    and s2=sin(2*w)), the roots are found by setting
%%     q=q'=0 and eliminating H^5.

if 1
  i0=(2*n-1)*pi/2;
  I=i0+[0 .8*pi/2];
%    I=[(2*n-1)*pi/2 n*pi]%r_fxn(I)
  %  I=[(n-1) n]*pi;
  w=fzero(@r_fxn,I,optimset('TolX',1e-8));%[I w]/pi
else
  w0=(2*n-1)*pi/2;
  tol=1e-12;
  err=2*tol;
  while err>tol
    [y,dy]=r_fxn(w0);
    dw=y/dy;
    w=w0-dw;%{n w/pi}
    w0=w;
    err=abs(dw);
  end
end

num=2-2*cos(2*w)+w*sin(2*w);
den=2*w*cos(2*w)-sin(2*w);
H5=8*w^5*num/den;
H=H5^.2;
del=-( (w/H)^4+H*cot(w)/w );
K=1i*w/H;

function [y,dy]=r_fxn(w)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to be "zeroed" to get triple root
c2=cos(2*w);
s2=sin(2*w);
lhs=(4*w+2*s2).*(2-2*c2+w.*s2);
rhs=(1-c2).*(2*w.*c2-s2);
y=rhs-lhs;
%%
d_lhs=(4+4*c2).*(2-2*c2+w.*s2) + (4*w+2*s2).*(5*s2+2*w.*c2);
d_rhs=2*s2.*(2*w.*c2-s2) + (1-c2).*(-4*w.*s2);
dy=d_rhs-d_lhs;