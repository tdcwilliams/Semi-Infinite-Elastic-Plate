function Nroots=...
GEN_Nroots_complex_polygon(FXN,zvec,varargin)

Nroots=0;
%  if 0%%test complex differentiation routine
%    z=2*exp(i*pi/7);
%    x=real(z);
%    y=imag(z);
%    f=cos(z);
%    theta=angle(f);%[f abs(f)*exp(1i*theta)]
%    df=-sin(z)
%    %%
%    eps=1e-15;
%    phi=theta+[0 -pi/2];
%    ztilde=z+i*eps*exp(i*phi);
%    df_ap=imag(cos(ztilde)*exp(-i*theta)/eps)*[1;1i]
%  end
%  
%  if 0%%test complex differentiation routine
%    z=2*exp(i*pi/7);
%    x=real(z);
%    y=imag(z);
%    f=cos(z);
%    u=real(f);
%    v=imag(f);
%    df=-sin(z)
%    %%
%    eps=1e-3;
%    ztilde=z+i*eps;
%    utilde=real(cos(ztilde));imag(utilde/eps)
%    vtilde=imag(cos(ztilde));
%    cos(ztilde)/eps
%  end

DO_TEST=0;
if nargin==0
  zvec=[1+1i;-1+1i;-1-1i;1-1i];
  DO_TEST=1;
end
Ratio=15;%% = ratio of points to winding no;
Narms=length(zvec);
zvec=[zvec;zvec(1)];
critter=1;
ratio=Ratio;
Npts=70;

while critter
  Npts=ceil(Ratio/ratio*Npts);
  tt=(0:Npts-1)'/Npts;
  Ntot=Npts*Narms;
  zz=zeros(Ntot+1,1);

if DO_TEST
  char={'-k','-r','-b','-g'};
  fxn=...
    inline('zz.^2.*sin(zz-.2).*(zz-.2).^17.*(zz+.6).*cos(14*zz+.4i)');
end

  for j=1:Narms
    zz0=zvec(j) + (zvec(j+1)-zvec(j))*tt;
    zz( (1:Npts)+(j-1)*Npts )=zz0;
%  ff0=fxn(zz0);
%  GEN_complex_plot( ff0,char{j} );
%  hold on, GEN_complex_plot(ff0(1),'*r'), pause;
  end

  zz(end)=zz(1);
  if ~DO_TEST
%  varargin{:}
    FF=feval(FXN,zz,varargin{:});
  else
    FF=fxn(zz);
  end
  arg=angle(FF);
  Arg=arg;
%  plot(arg),pause

  Jswp=find(abs(arg(2:end)-arg(1:end-1))>pi);
  chg=0;
  for j=1:length(Jswp)
    j0=Jswp(j);
    arg0=arg(j0);
    arg1=arg(j0+1);
    if arg0>arg1
      Arg(j0+1:end)=Arg(j0+1:end)+2*pi;
    elseif arg0<arg1
      Arg(j0+1:end)=Arg(j0+1:end)-2*pi;
    end
  end
%  plot(arg/2/pi,'g'), hold on;
%  plot(Arg/2/pi)
  Nroots=abs(round((Arg(end)-Arg(1))/2/pi));
  ratio=Ntot/Nroots;
  critter=(ratio<Ratio);
end