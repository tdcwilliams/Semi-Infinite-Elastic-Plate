function y=GEN_findroots_complex_rectangle(FXN,verts,varargin);

if nargin==0
  FXN=@tester;
  verts={0,[1 1]};
end

centre=verts{1};
ab=verts{2};
c=max(ab);
jmax=find(ab==c);
jmin=find(ab~=c);
z0=a+1i*b;
ERR=abs(z0);
TOL=1e-12;
%%
zvec0=[z0;-z0';-z0;z0'];
zvec=centre+zvec0;
Nroots=GEN_Nroots_complex_polygon(FXN,zvec,varargin)
return

if Nroots==0
  y=[];
  return;
elseif ERR<TOL
  y(1:Nroots)=centre;
  return;
end

if a>=b%% split vertically
  verts1={centre+a/2,[a/2,b]};
  y1=GEN_findroots_complex_rectangle(FXN,verts1,varargin);
  verts2={centre-a/2,[a/2,b]};
  y2=GEN_findroots_complex_rectangle(FXN,verts2,varargin);
  y=[y1;y2];
elseif
  verts1={centre+.5i*b,[a,b/2]};
  y1=GEN_findroots_complex_rectangle(FXN,verts1,varargin);
  verts2={centre-.5i*b,[a,b/2]};
  y2=GEN_findroots_complex_rectangle(FXN,verts2,varargin);
  y=[y1;y2];
end

function y=tester(zz)
y=(zz-3).^2.*sin(zz-.2).*(zz-.2).^1.*(zz+.6).*cos(14*zz+.4i)
