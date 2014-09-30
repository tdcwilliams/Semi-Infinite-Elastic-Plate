function K = RTS_wtr_roots(lam,H,N)
%% CALL: K=RTS_wtr_roots(lam,H,N)
%% calculates the real root and N imaginary solutions of 
%% the dispersion relation for open water, 1/lam = K*tanh(K*H).
%% DEPENDENCIES: RTS_imag_roots_wtr.c (.mexglx), RTS_imag_root_wtr.c (.mexglx).

USE_MEX=0;
if USE_MEX
  FXN1   = @RTS_imag_roots_wtr;
  FXN2   = @RTS_imag_root_wtr;
else
  FXN1   = @RTS_imag_roots_wtr_matlab;
  FXN2   = @RTS_imag_root_wtr_matlab;
end

do_test  = 0;
%% GET REAL ROOT:
tol   = 1e-8;
ans1  = 1/lam;
k     =  ans1 - f_on_df(ans1,lam,H);
test  = abs(k-ans1)/abs(k);
while (test > tol)
  ans1   = k;
  k      = ans1 - f_on_df(ans1,lam,H);
  test   = abs(k-ans1)/abs(k);
end
K  = k;

if N>0%% GET IMAGINARY ROOTS:
%    w=RTS_imag_roots_wtr(lam,H,N);
  w   = feval(FXN1,lam,H,N);
    %%w=-i*K*H (real for imaginary K);
  j_notconverged  = find(w==0);
  for r=1:length(j_notconverged)
    j    = j_notconverged(r);
    i1   = j*pi;
    i0   = i1-pi/2;
    w(j) = feval(FXN2,lam,H,i0,i1);
%     w(j)=RTS_imag_root_wtr(lam,H,i0,i1);
  end
  K   = [ k; 1i*w/H ];

  if do_test%% do some tests
    Del  = pi/H;
    subplot(1,2,1), plot(real(K),imag(K)/Del,'x'), hold on;
    plot(-real(K),-imag(K)/Del,'x');
    Xlim = (K(1)+1)*[-1 1];
    Ylim = 5*[-1 1];
    plot(Xlim,0*Xlim,':r')
    plot(0*Xlim,Ylim,':r'), hold off;
    xlim(Xlim), ylim(Ylim);
    %%
    n       = 1000;
    NP      = 15;
    wf      = (0:n)'*NP*Del/n;
    pf      = lam*wf.*sin(wf*H)+cos(wf*H);
    scale   = abs(wf)+1;
    subplot(1,2,2), plot(wf/Del,[pf./scale,0*wf]);
    hold on;
    %%
    jj   = find(abs(real(K)) < tol);
    plot(real(-i*K(jj))/Del,0*K(jj),'.r'); hold off;
    xlim([0 NP]);
  end
end%%%

function y = f_on_df(K,lam,H)
%y=f/f', f being a multiple of the dispersion function;
%this is the correction term in Newton's method for finding zeros in f

if abs(real(K))<=7.5
  f   = lam*K*sinh(K*H)-cosh(K*H);
  df  = lam*K*H*cosh(K*H)+(lam-H)*sinh(K*H);
else
  f   = lam*K*tanh(K*H)-1;
  df  = lam*K*H+(lam-H)*tanh(K*H);
end
y  = f/df;
