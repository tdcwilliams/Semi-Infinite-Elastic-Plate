function out=SF_gamma_cx_v1(z,TAKE_LOG)
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ==========================================================
%     Purpose: This program computes the gamma function â(z)
%     or ln[â(z)]for a complex argument using
%     subroutine CGAMA
%     CALL: out=SF_gamma_cx_v1(z,TAKE_LOG)
%     Input :
%       z, a complex number;
%       TAKE_LOG==0 -> out = â(z)
%       TAKE_LOG==1 -> out = ln[â(z)]
%       TAKE_LOG==2 -> out = (e./z).^z.*â(z)
%
%     Examples:
%     x         y           Re[â(z)]Im[â(z)]
%     --------------------------------------------------------
%     2.50      5.00     .2267360319D-01    -.1172284404D-01
%     5.00     10.00     .1327696517D-01     .3639011746D-02
%     2.50     -5.00     .2267360319D-01     .1172284404D-01
%     5.00    -10.00     .1327696517D-01    -.3639011746D-02
%     x         y          Re[lnâ(z)]Im[lnâ(z)]
%     ---------------------------------------------------------
%     2.50      5.00    -.3668103262D+01     .5806009801D+01
%     5.00     10.00    -.4285507444D+01     .1911707090D+02
%     2.50     -5.00    -.3668103262D+01    -.5806009801D+01
%     5.00    -10.00    -.4285507444D+01    -.1911707090D+02
%     ==========================================================

x=real(z);
y=imag(z);
gr=0;
gi=0;

if nargin==1
  TAKE_LOG=0;
end

kf=~TAKE_LOG;
[gr,gi]=cgama(x,y,kf);
out=gr+1i*gi;
if TAKE_LOG==2
  log_out=out+z.*(1-log(z))+log(z);
  out=exp(log_out);
  out(find(z==0))=1;
end
return;


function [gr,gi]=cgama(x,y,kf);
%     =========================================================
%     Purpose: Compute the gamma function â(z)or ln[â(z)]
%     for a complex argument
%     Input :  x  --- Real part of z
%     y  --- Imaginary part of z
%     KF --- Function code
%     KF=0 for ln[â(z)]
%     KF=1 for â(z)
%     Output:  GR --- Real part of ln[â(z)]or â(z)
%     GI --- Imaginary part of ln[â(z)]or â(z)
%     ========================================================
a=zeros(1,10);
x1=0.0;
pi=3.141592653589793d0;
gr=0*x;
gi=0*x;

%% coefficients of expansion of \log(\G):
a(:)=[8.333333333333333d-02,...
      -2.777777777777778d-03,...
      7.936507936507937d-04,...
      -5.952380952380952d-04,...
      8.417508417508418d-04,...
      -1.917526917526918d-03,...
      6.410256410256410d-03,...
      -2.955065359477124d-02,...
      1.796443723688307d-01,...
      -1.39243221690590d+00];

%% if a negative integer:
jneg_int=find(y == 0.0d0&x == fix(x)&x <= 0.0d0);
jOK=find(~(y == 0.0d0&x == fix(x)&x <= 0.0d0));
gr(jneg_int)=Inf;
x(jneg_int)=[];
y(jneg_int)=[];
%%
if isempty(jOK)
  return;
end

%%case 1: x<0
j_neg_hp=find(x<0);
%  x1=x(j_neg_hp);%%x1 now original x;
%  y1=y(j_neg_hp);
x(j_neg_hp)=-x(j_neg_hp);%% now x>0;
y(j_neg_hp)=-y(j_neg_hp);
x0=x;

%% make x bigger if x too small:
j_small=find(x <= 7.0);
%  x(j_small)
na=fix(7-x(j_small));
x0(j_small)=x(j_small)+na;
%%
z1=sqrt(x0.*x0+y.*y);
th=atan(y./x0);
gr(jOK)=(x0-.5d0).*log(z1)-th.*y-x0+0.5d0.*log(2.0d0.*pi);
gi(jOK)=th.*(x0-0.5d0)+y.*log(z1)-y;
for  k=1:10;
  t=z1.^(1-2.*k);
  gr(jOK)=gr(jOK)+a(k).*t.*cos((2.0d0.*k-1.0d0).*th);
  gi(jOK)=gi(jOK)-a(k).*t.*sin((2.0d0.*k-1.0d0).*th);
end
k=10+1;
%  tstg0=[gr(jOK),gi(jOK)]

%% reverse changes made to small x:
gr1=0*x(j_small);
gi1=0*x(j_small);
NA=max(na);
for  j=0:NA-1;
%    j
  FAC=(j<na);
  gr1=gr1+.5d0.*FAC.*log((x(j_small)+j).^2+y(j_small).*y(j_small));
  gi1=gi1+FAC.*atan(y(j_small)./(x(j_small)+j));
end
%  gr1
%  j=na-1+1;
%  [gr(jOK(j_small)),gr1]
gr(jOK(j_small))=gr(jOK(j_small))-gr1;
gi(j_small)=gi(j_small)-gi1;
%  tstg1=[gr(jOK),gi(jOK)]

%% reverse changes made to x<0:
z1=sqrt(x(j_neg_hp).*x(j_neg_hp)+y(j_neg_hp).*y(j_neg_hp));
th1=atan(y(j_neg_hp)./x(j_neg_hp));
sr=-sin(pi.*x(j_neg_hp)).*cosh(pi.*y(j_neg_hp));
si=-cos(pi.*x(j_neg_hp)).*sinh(pi.*y(j_neg_hp));
z2=sqrt(sr.*sr+si.*si);
th2=atan(si./sr);
th2(find(sr < 0.0d0))=pi+th2(find(sr < 0.0d0));
%%
JJ=jOK(j_neg_hp);
gr(JJ)=log(pi./(z1.*z2))-gr(JJ);
gi(JJ)=-th1-th2-gi(JJ);

%% if kf==1, exponentiate to get
%% \G(z) itself
if(kf == 1);
  g0=exp(gr(jOK));
  gr(jOK)=g0.*cos(gi(jOK));
  gi(jOK)=g0.*sin(gi(jOK));
end;
return;
%  end