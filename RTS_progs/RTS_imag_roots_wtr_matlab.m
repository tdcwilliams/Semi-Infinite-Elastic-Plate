function x=RTS_imag_roots_wtr_matlab(lam,H,N1)

EPS=1e-9;
MAXIT=50;
N0=1;
x=zeros(N1+1-N0,1);

for j=N0:N1
  %% initial guess
  w=j*pi;
  i1=w;
  i0=i1-pi/2;
  for its=1:MAXIT
    w4=w*w*w*w;
    H4=H*H*H*H;
    p=lam*w*sin(w)+H*cos(w);
    dp=(lam-H)*sin(w)+lam*w*cos(w);
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