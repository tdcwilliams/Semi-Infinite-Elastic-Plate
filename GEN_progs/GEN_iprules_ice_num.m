function M=GEN_iprules_ice_num(gam1,gam2,zlims,zfloors,HH)

do_test=0;
if do_test
  gam1=[1 1.3]';
  gam2=[1 1.6]';
  zlims=[-1 0];
  zfloors=[-1 -1.5];
  HH=[1 2];
end

z0=zlims(1);
z1=zlims(2);
dz=z1-z0;
%%
zf1=zfloors(1);
zf2=zfloors(2);
%%
H1=HH(1);
H2=HH(2);
[Gam2,Gam1]=meshgrid(gam2,gam1);
[sonc10,conc10]=get_sonc_conc(gam1,[z0-zf1,H1]);
[sonc11,conc11]=get_sonc_conc(gam1,[z1-zf1,H1]);
[sonc20,conc20]=get_sonc_conc(gam2,[z0-zf2,H2]);
[sonc21,conc21]=get_sonc_conc(gam2,[z1-zf2,H2]);
%%
Sum=Gam1+Gam2;
Mp=( sonc11*conc21.' + conc11*sonc21.' )./(2*Sum);
Mp=Mp-( sonc10*conc20.' + conc10*sonc20.' )./(2*Sum);
%%
Mm=sonc11*conc21.' - conc11*sonc21.';
Mm=Mm-( sonc10*conc20.' - conc10*sonc20.' );
Diff=Gam1-Gam2;
jp=find(Diff~=0);
Mm(jp)=Mm(jp)./(2*Diff(jp));
%%
j0=find(Diff==0);
if ~isempty(j0)
  [sonc1f,conc1f]=get_sonc_conc(gam1,[zf1,H1]);
  [sonc2f,conc2f]=get_sonc_conc(gam2,[zf2,H2]);
  Mm0=conc1f*conc2f.'-sonc1f*sonc2f.';
  Mm(j0)=dz/2*Mm0(j0);
end
%%
M=Mp+Mm;

if do_test
  m=length(gam1);
  n=length(gam2);
  zz=(0:-.1:-1)';
%    ig=test_ig(zz,gam1(m),gam2(n),zfloors,HH);
%    ig=[ig, cosh(gam1(m)*(zz-zf1)).*cosh(gam2(n)*(zz-zf2))/...
%              ( cosh(gam1(m)*H1)*cosh(gam2(n)*H2) ) ];
%    ig=[ig,test_ig1(zz,gam1(m),gam2(n),zfloors,HH,[1 1])]
  for j=1:m
    for r=1:n
      Mtst=[quad(@(z)test_ig(z,gam1(j),gam2(r),zfloors,HH),z0,z1),M(j,r)]
%    Mtst=[Mtst(1),antidiff(z1,gam1(m),gam2(n),zfloors,HH,[1 1])-...
%  	antidiff(z0,gam1(m),gam2(n),zfloors,HH,[1 1]),Mtst(2)]

%    Mtst=[quad(@(z)test_ig1(z,gam1(m),gam2(n),zfloors,HH,[1 0]),z0,z1),Mp(m,n)];
%    Mtst=[Mtst(1),antidiff(z1,gam1(m),gam2(n),zfloors,HH,[1 0])-...
%    	antidiff(z0,gam1(m),gam2(n),zfloors,HH,[1 0]),Mtst(2)]
%    antidiff(z1,gam1(m),gam2(n),zfloors,HH,[1 0])
  %Mtst=[quad(@(z)test_ig1(z,gam1(m),gam2(n),zfloors,HH,[0 1]),z0,z1),Mm(m,n)]
    end
  end
end

function [sonc,conc]=get_sonc_conc(gam,ZH);
%%sonc=sinh(gam*Z)/cosh(gam*H), conc=cosh(gam*Z)/cosh(gam*H);
Z=ZH(1);
H=ZH(2);
denom=1+exp(-2*gam*H);
sonc=( exp(-gam*(H-Z))-exp(-gam*(H+Z)) )./denom;
conc=( exp(-gam*(H-Z))+exp(-gam*(H+Z)) )./denom;

function y=test_ig(zz,gam1,gam2,zfloors,HH)
%%y=cosh(gam1*(z-zf1))/cosh(gam1*H1)*cosh(gam2*(z-zf2))/cosh(gam2*H2)
for j=1:length(zz)
  z=zz(j);
  [sonc,conc1(j,1)]=get_sonc_conc(gam1,[z-zfloors(1),HH(1)]);
  [sonc,conc2(j,1)]=get_sonc_conc(gam2,[z-zfloors(2),HH(2)]);
end
y=conc1.*conc2;

function y=test_ig1(zz,gam1,gam2,zfloors,HH,coeffs)

zf1=zfloors(1);
zf2=zfloors(2);
H1=HH(1);
H2=HH(2);
c1=coeffs(1);
c2=coeffs(2);

for j=1:length(zz)
  z=zz(j);
  sg=gam1+gam2;
  dg=gam1-gam2;
  y(j,1)=c1*cosh( sg*z-(gam1*zf1+gam2*zf2) )+...
			c2*cosh( dg*z-(gam1*zf1-gam2*zf2) );
end
y=y/2/cosh(gam1*H1)/cosh(gam2*H2);

function y=antidiff(z,gam1,gam2,zfloors,HH,coeffs)

zf1=zfloors(1);
zf2=zfloors(2);
H1=HH(1);
H2=HH(2);
c1=coeffs(1);
c2=coeffs(2);

sg=gam1+gam2;
dg=gam1-gam2;
y=c1/sg*sinh( sg*z-(gam1*zf1+gam2*zf2) )+...
			c2/dg*sinh( dg*z-(gam1*zf1-gam2*zf2) );
y=y/2/cosh(gam1*H1)/cosh(gam2*H2);