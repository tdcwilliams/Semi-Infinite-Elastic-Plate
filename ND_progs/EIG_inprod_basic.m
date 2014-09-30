function Y=EIG_inprod_basic(gam1,gam2,HH)
%% CALL: Y=EIG_inprod_basic(gam1,gam2,HH)
%% HH=[H_1,H_2];
%% HM=min(HH);
%% Y(j,r)=\int_{-HM}^0.varf_1(z,gam1(j)).varf_2(z,gam2(r)).dz,
%% where
%% varf_n(z,gam)=cosh(gam*(z+HM))/cosh(gam*H_n).

H1=HH(1);
H2=HH(2);
HM=min(HH);
%%
[Gam2,Gam1]=meshgrid(gam2,gam1);
Sum=Gam2+Gam1;
Diff=Gam2-Gam1;
jp=find(abs(Diff)>0);
jz=find(abs(Diff)==0);
%%
Y=0*Sum;
Y(jp)=.5*sinh(Sum(jp)*HM)./Sum(jp) +...
             .5*sinh(Diff(jp)*HM)./Diff(jp);
Y(jz)=.5*sinh(Sum(jz)*HM)./Sum(jz) + HM/2;
Y=Y./cosh(Gam1*H1)./cosh(Gam2*H2);