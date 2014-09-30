function x=RTS_imag_roots_ice_matlab(del,H,N0,N1)

EPS=1e-9;
MAXIT=50;
x=zeros(N1+1-N0,1);

for j=N0:N1
  %% initial guess
  w=j*pi;
  i1=w;
  i0=i1-pi/2;
  if j==N0
    i0=i1-pi;
  end
  for its=1:MAXIT
    w4=w*w*w*w;
    H4=H*H*H*H;
    p=(w4+del*H4)*w*sin(w)+H4*H*cos(w);
    dp=( 5*w4+H4*(del-H) )*sin(w)+(w4+del*H4)*w*cos(w);
    w1=w;
    w=w1-p/dp;
    if abs(w-w1) <= EPS
      break;
    end
  end
  %%check root has converged & is in correct interval:
  if ( (its < MAXIT) & ((w-i0+EPS)*(i1+EPS-w)>=0) )
    x(j+1-N0)=w;
  else
    x(j+1-N0)=0;
  end
end