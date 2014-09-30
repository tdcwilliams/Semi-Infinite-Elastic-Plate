function [tn,wn]=GEN_numint_exp(N)
%% CALL: [tn,wn]=GEN_numint_exp(N);
%% does integral of exp(pi*i*n*t) over t=-1..1
%%  using 2N quadrature points;
%%  exact for -(2*N-1)<=n<=(2*N-1)
%%  => gives Fourier coefficients f_m
%%       for m=-(N-1)..(N-1);
%% NB this means if you want Fourier coefficients
%%  n=-N..N, you need 2*(N+1) points
%%   ie ask for GEN_numint_exp(N+1).

USE_CENTREPOINTS=1;%% CHANGE THIS TO '0'
                   %% IF WANT TO USE THE END POINTS:
jj=(-N:N-1)';%[length(jj),2*N]
Del=1/N;
wn=Del+0*jj;
tn=jj*Del+USE_CENTREPOINTS*Del/2;

if 0%%do test
  for jt=-2*N:2*N
    ff=exp(i*pi*jt*tn);
    int=sum(wn.*ff);
    disp(jt),disp(int)
  end
end