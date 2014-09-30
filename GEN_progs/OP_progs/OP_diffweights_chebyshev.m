function [ip_1,ip_t,ip_acos,ip_t_acos]=OP_diffweights_chebyshev(ipT)
%% integrates \int_{-1}^1w(t)f(t)dt=ip_w*[f(t_j)]
%% using the chebyshev quadrature points (t_j), when
%% w(t)=1, t, acos(t)/pi, t*acos(t)/pi;

if ~iscell(ipT)
  N=size(ipT,1)-1;
else
  N=ipT{1};
  [tt,ww]=OP_numint_chebyshev(N+1);
  ipT=OP_inprod_chebyshev(tt,ww,N);
end
nn=(0:N)';
nn_ev=(0:2:N)';
nn_odd=(1:2:N)';

%% ip_1*ff=\int_{1}^{1}f(t)dt=\sum_{n=0}^Nf_n\int_{-1}^{1}T_n(t)dt
int_n=zeros(N+1,1);
int_n(nn_ev+1)=2./(1-nn_ev.^2);
ip_1=int_n'*ipT;

%% ip_t*ff=\int_{1}^{1}[f(t)*t]dt=\sum_{n=0}^Nf_n\int_{-1}^{1}[T_n(t)*t]dt
int_n=zeros(N+1,1);
int_n(nn_odd+1)=2./(4-nn_odd.^2);
ip_t=int_n'*ipT;

%% ip_acos*ff=\int_{1}^{1}[f(t)*acos(t)/pi]dt
%% = \sum_{n=0}^Nf_n\int_{-1}^{1}[T_n(t)*acos(t)]dt
int_n=zeros(N+1,1)-1/4;
nn_=nn;
nn_(2)=[];
int_n(nn_+1)=(-1).^nn_./(1-nn_.^2);
ip_acos=int_n'*ipT;

%% ip4*ff=\int_{1}^{1}[f(t)*t*acos(t)/pi]dt
%% = \sum_{n=0}^Nf_n\int_{-1}^{1}[T_n(t)*t*acos(t)]dt
int_n=zeros(N+1,1)-1/16;
nn_=nn;
nn_(3)=[];
int_n(nn_+1)=-(-1).^nn_./(4-nn_.^2);
ip_t_acos=int_n'*ipT;