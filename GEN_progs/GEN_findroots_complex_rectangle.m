function [y,Nroots]=...
  GEN_findroots_complex_rectangle(FXN,verts,...
    rts_found,varargin);

DO_TEST=0;
if nargin==0
  DO_TEST=1;
  FXN=@tester;
  verts={0.0001,1,1};
  prams=[];
end
prams=[];

centre=verts{1};
a=verts{2};
b=verts{3};
z0=a+1i*b;
ERR=abs(z0);
TOL=1e-12;
%%
zvec0=[z0;-z0';-z0;z0'];
zvec=centre+zvec0;
Nroots=GEN_Nroots_complex_polygon(...
         FXN,zvec,varargin{:});
%  {centre,rts_found,Nroots},pause
%  zvec,{Nroots,centre,ERR,ERR<TOL},pause
Nfound=length(rts_found);%,pause

if Nroots<=Nfound
  %% if we've already found all the roots in
  %% the box no need to iterate further:
  y=[];
  return;
elseif Nroots==1
  %% use Muller's method to see if we can find
  %% the root faster:
  x3=centre+[z0,0,-z0];
  y_muller=...
    GEN_findroots_muller(FXN,x3,varargin{:});
  if is_inside(y_muller,{centre,a,b})
    y=y_muller;
    return;
  end
elseif ERR<TOL
  %% box is small enough to just
  %% take the centre as the root:
  y(1:Nroots-Nfound,1)=centre;
  return;
end

if a>=b%% split vertically:
  verts1={centre+a/2,a/2,b};
  j1=find(is_inside(rts_found,verts1));
  rf1=rts_found(j1);
  %%
  verts2={centre-a/2,a/2,b};
  j2=find(is_inside(rts_found,verts2));
  rf2=rts_found(j2);
  %%
  y1=GEN_findroots_complex_rectangle(...
       FXN,verts1,rf1,varargin{:});
  y2=GEN_findroots_complex_rectangle(...
       FXN,verts2,rf2,varargin{:});
  y=sort([y1;y2]);
else%%split horizontally:
  Verts1={centre+.5i*b,a,b/2};
  j1=find(is_inside(rts_found,Verts1));
  rf1=rts_found(j1);
  %%
  Verts2={centre-.5i*b,a,b/2};
  j2=find(is_inside(rts_found,Verts2));
  rf2=rts_found(j2);
  %%
  Y1=GEN_findroots_complex_rectangle(...
       FXN,Verts1,rf1,varargin{:});
  Y2=GEN_findroots_complex_rectangle(...
       FXN,Verts2,rf2,varargin{:});
  y=sort([Y1;Y2]);
end

if length(y)~=Nroots
%    disp('warning: not all roots found')
%    [length(y),Nroots]
%    verts,pause
end

function y=tester(zz,prams)
y=(zz-.3).^2.*sin(zz-.2).*(zz-.2).^9.*...
    (zz+.6).*cos(14*zz+.4i).*(zz-.2-1e-9i).^4;

function critter=is_inside(y,verts);

centre=verts{1};
a=verts{2};
b=verts{3};
%%
cre=real(centre);
cim=imag(centre);
yre=real(y);
yim=imag(y);
critter= (yre>cre-a & yre<cre+a) & ...
         (yim>cim-b & yim<cim+b);