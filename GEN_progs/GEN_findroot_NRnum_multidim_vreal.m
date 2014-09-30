function x=GEN_findroot_NRnum_multidim_vreal(...
             FXN,guess,varargin)

N=length(guess);
I_N=eye(N);
DF=0*I_N;
critter=1;
x0=guess;
eps=1e-12;
TOL=1e-14;

MAXITS=100;
its=0;

while critter
  its=its+1;
  FF=feval(FXN,x0,varargin{:});

  %% GET JACOBIAN:
  for j=1:N
    xtilde=x0+1i*eps*I_N(:,j);
    Ftilde=feval(FXN,xtilde,varargin{:});
    DF(:,j) = imag( Ftilde/eps);
  end

  dx=-DF\FF;
  x=x0+dx;
  ERR=max(abs(dx));
  x0=x;
  critter=(ERR>TOL & its<MAXITS);
end

if its==MAXITS
  disp('warning (GEN_findroot_NRnum_multidim_vreal.m):');
  disp('root not converged');
end