function x=GEN_findroot_NRnum_multidim(...
             FXN,guess,varargin)

if nargin==0
  FXN=@tester;
  guess=[3;4.6];
  x_target=[pi;1.5*pi]
end


N=length(guess);
I_N=eye(N);
DF=0*I_N;
critter=1;
x0=guess;
eps=1e-5;
TOL=1e-8;

MAXITS=500;
its=0;

while critter
  its=its+1;
  FF=feval(FXN,x0,varargin{:});

  %% GET JACOBIAN:
  for j=1:N
    xtilde2=x0+eps*I_N(:,j);
    xtilde1=x0-eps*I_N(:,j);
    Ftilde2=feval(FXN,xtilde2,varargin{:});
    Ftilde1=feval(FXN,xtilde1,varargin{:});
    DF(:,j) = (Ftilde2-Ftilde1)/(2*eps);
  end

  dx=-DF\FF;
  x=x0+dx;
  x0=x;
%  abs(dx),pause
  ERR=max(abs(dx));
  critter=(ERR>TOL & its<MAXITS);
end

if its==MAXITS
  disp('warning (GEN_findroot_NRnum_multidim_vreal.m):');
  disp('root not converged');
end

function f=tester(xy,prams)

x=xy(1);
y=xy(2);

f=[sin(x)*cos(y);(x-pi)*(y-3*pi/2)];