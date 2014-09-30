function [y,sing_coeffs]...
=BES_diff_hankel_series(r,A,gam,n,keepsing)

if nargin==4
  keepsing=0;
end
sing_coeffs=zeros(n+1,1);

if n==0
  [y,sing_coeffs]=BES_hankel_series(r,A,gam,0,keepsing);
	%%~c1*log(r)
elseif n==1
  [y,sc]=BES_hankel_series(r,A,gam,1,keepsing);
  	%%~c1*log(r)+c2/r^2
  sing_coeffs=[0;sc(2)];
  y=r.*y;
  if keepsing==0
    jp=find(r);
    y(jp)=y(jp)+sc(1)*r(jp).*log(r(jp));
  end
elseif n==2
  [Del2W,sc1]=BES_hankel_series(r,-A.*gam.^2,gam,0,keepsing);
	%%~c1*log(r)
  [negDW,sc]=BES_hankel_series(r,-A,gam,1,keepsing);
	%%~c1*logr+c2/r^2
  sing_coeffs(1)=sc(1)+sc1;
  sing_coeffs(3)=sc(2);
  y=Del2W+negDW;
elseif n==3
  [D1Del2W,sc1]=BES_hankel_series(r,-A.*gam.^2,gam,1,keepsing);
	%%~c1*logr+c2/r^2
  [negD2W,sc]=BES_hankel_series(r,-A,gam,2,keepsing);
	%%~c1*log(r)+c2/r^2+c3/r^4
  sing_coeffs(2)=sc(2)+sc1(2);
  sing_coeffs(4)=sc(3);
  y=r.*(D1Del2W+negD2W);
  if keepsing==0
    jp=find(r);
    C=sc(1)+sc1(1);
    y(jp)=y(jp)+C*r(jp).*log(r(jp));
  end
elseif n==4
  Cpoly=zeros(2,1);%%poly is deg 1 in r^2
  Cpoly_log=zeros(2,1);%%logterm is deg 1 in r^2 times r^2*log(r)
  jj=1+[2 4];%%1/r^m coeffs; m=2,4
  r2=r.^2;
  %%
  [facD2W,sc2]=BES_hankel_series(r,3*A,gam,2,keepsing);
	%%~c1*lg+c2/r^2+c3/r^4
  y=facD2W;
  sing_coeffs([1,jj])=sc2;%%log coeff as well
  %%
  [facD3W,sc3]=BES_hankel_series(r,6*A,gam,3,keepsing);
	%%~[c1*lg+c2/r^2+c3/r^4+c4/r^6]; r^m->r^(m+2)
  y=y+r2.*facD3W;
  Cpoly(2)=sc3(2);%%r^0
  Cpoly_log(2)=sc3(1);%%r^2log(r)
  sing_coeffs(jj)=sing_coeffs(jj)+sc3(3:4);
  %%
  [D4W,sc4]=BES_hankel_series(r,A,gam,4,keepsing);
	%%~c1*lg+c2/r^2+c3/r^4+c4/r^6+c5/r^8;r^m->r^(m+4)
  y=y+r2.^2.*D4W;
  Cpoly(1:2)=Cpoly(1:2)+sc4(2:3);%%r^2,r^0
  Cpoly_log(1)=sc4(1);%%r^4log(r)
  sing_coeffs(jj)=sing_coeffs(jj)+sc4(4:5);
  %%
  %y=facD2W+r.^2.*facD3W+r.^4.*D4W;
  if keepsing==0
    poly=Cpoly(1)*r2+Cpoly(2);
    y=y+poly;
    %%
    jp=find(r);
    logr=log(r(jp));
    r2=r2(jp);
    logterm=r2.*logr.*(Cpoly_log(1)*r2+Cpoly_log(2));
    y(jp)=y(jp) + logterm;
  end
elseif n==5
  Cpoly=zeros(2,1);%%poly is deg 1 in r^2 times r
  Cpoly_log=zeros(3,1);%%poly is deg 2 in r^2 times r.log(r)
  jj=1+[1 3 5];%%1/r^m sing's, m=1,3,5
  r2=r.^2;
  %%
  [facD3W,sc3]=BES_hankel_series(r,15*A,gam,3,keepsing);
	%%~c1.log(r)+c2/r^2+c3/r^4+c4/r^6
  y=r.*facD3W;
  sing_coeffs(jj)=sc3(2:4);
  Cpoly_log(3)=sc3(1);
  %%
  [facD4W,sc4]=BES_hankel_series(r,10*A,gam,4,keepsing);
	%%~c1.log(r)+c2/r^2+c3/r^4+c4/r^6+c5/r^8
  y=y+r.^3.*facD4W;
  sing_coeffs(jj)=sing_coeffs(jj)+sc4(3:5);
  Cpoly_log(2)=sc4(1);
  Cpoly(2)=sc4(2);
  %%
  [D5W,sc5]=BES_hankel_series(r,A,gam,5,keepsing);
	%%~c1.log(r)+c2/r^2+c3/r^4+c4/r^6+c5/r^8+c6/r^10
  y=y +r.^5.*D5W;
  sing_coeffs(jj)=sing_coeffs(jj)+sc5(4:6);
  Cpoly_log(1)=sc5(1);
  Cpoly(1:2)=Cpoly(1:2)+sc5(2:3);
  if keepsing==0
    poly=( Cpoly(1)*r2+Cpoly(2) ).*r;
    y=y+poly;
    %%
    jp=find(r);
    logr=log(r(jp));
    r2=r2(jp);
    poly_log=( Cpoly_log(1)*r2+Cpoly_log(2) ).*r2 + Cpoly_log(3);
    y(jp)=y(jp) + r(jp).*logr.*poly_log;
  end
elseif n==6
  Cpoly=zeros(3,1);%%poly is deg 2 in r^2
  Cpoly_log=zeros(3,1);%%poly is deg 2 in r^2 times r^2log(r)
  jj=1+[2 4 6];
  r2=r.^2;
  %%
  [facD3W,sc3]=BES_hankel_series(r,15*A,gam,3,keepsing);
	%%~c1.log(r)+c2/r^2+c3/r^4+c4/r^6
  y=facD3W;
  sing_coeffs([1,jj])=sc3;
  %%
  [facD4W,sc4]=BES_hankel_series(r,45*A,gam,4,keepsing);
	%%~c1.log(r)+c2/r^2+c3/r^4+c4/r^6+c5/r^8
  y=y+r2.*facD4W;
  sing_coeffs(jj)=sing_coeffs(jj)+sc4(3:5);
  Cpoly_log(3)=sc4(1);
  Cpoly(3)=sc4(2);
  %%
  [facD5W,sc5]=BES_hankel_series(r,15*A,gam,5,keepsing);
	%%~c1.log(r)+c2/r^2+c3/r^4+c4/r^6+c5/r^8+c6/r^10
  y=y+r2.^2.*facD5W;
  sing_coeffs(jj)=sing_coeffs(jj)+sc5(4:6);
  Cpoly_log(2)=sc5(1);
  Cpoly(2:3)=Cpoly(2:3)+sc5(2:3);
  %%
  [D6W,sc6]=BES_hankel_series(r,A,gam,6,keepsing);
	%%~c1.log(r)+c2/r^2+c3/r^4+c4/r^6+c5/r^8+c6/r^10+c7/r^12
  y=y+r2.^3.*D6W;
  sing_coeffs(jj)=sing_coeffs(jj)+sc6(5:7);
  Cpoly_log(1)=sc6(1);
  Cpoly(1:3)=Cpoly(1:3)+sc6(2:4);
  %%
  %y=facD3W + r.^2.*facD4W + r.^4.*facD5W + r.^6.*D6W;
  if keepsing==0
    poly=( Cpoly(1)*r2+Cpoly(2) ).*r2 + Cpoly(3);
    y=y+poly;
    %%
    jp=find(r);
    logr=log(r(jp));
    r2=r2(jp);
    poly_log=( Cpoly_log(1)*r2+Cpoly_log(2) ).*r2 + Cpoly_log(3);
    y(jp)=y(jp) + r2.*logr.*poly_log;
  end
end

