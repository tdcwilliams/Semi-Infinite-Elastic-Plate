function Y=EIG_inprod(gam1,gam2,Z1,Z2)
%% CALL: Y=EIG_inprod(gam1,gam2,Z1,Z2)
%% Z1={Dr1,del1,H};
%% Z2={Dr2,del2,H};
%% Y(j,r)=-B10(gam1(j)*\int_{-H}^0.varf_1(z,gam1(j))*...
%%     varf_2(z,gam2(r)).dz - ...
%%       B10(gam1(j)*Dr1*( gam1(j)^2+gam2(r)^2 ),
%% where
%% varf_n(z,gam)=cosh(gam*(z+H))/cosh(gam*H).

Dr1=Z1{1};
del1=Z1{2};
H=Z1{3};
%%
Dr2=Z2{1};
del2=Z2{2};
%%
[Gam2,Gam1]=meshgrid(gam2,gam1);
Lam11=Dr1*Gam1.^4+del1;
B11=2*Gam1.*Lam11.*calc_res(Z1,Gam1);
%%
Lam12=Dr1*Gam2.^4+del1;
Lam22=Dr2*Gam2.^4+del2;
%%
Diff=Gam2-Gam1;
jp=find(abs(Diff)>0);
jz=find(abs(Diff)==0);
%%
Y=0*Sum;
Y(jz)=1;
Y(jp)=B11(jp).*( 1-Lam12(jp)./Lam22(jp) )./(Gam2.^2-Gam1.^2);

function y=calc_res(Z2,gam)
%% y=calc_res(Z2,hr,gamma)=Res(1/f(K),gamma_n),
%% where gamma_n is a root of the dispersion relation
%% f=1/K/tanh(KH)-(Dr*K^4+del);
%% Z2={Dr,del,H}.

Dr=Z2{1};
del=Z2{2};
H=Z2{3};
%%
Lam=Dr*gam.^4+del;
Lampr=Lam+4*Dr*gam.^4;
denom=H*(Lam.^2.*gam.^2-1)+Lampr;
y=-gam./denom;