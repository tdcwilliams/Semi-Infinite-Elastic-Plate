function [Matrix_left,Matrix_right] = ...
  OPcheby_ODE_matrix_constcoeff(A,B,C,BCS,NT)
%% function to compute matrix corresponding to
%% operator L=A\pa_y^2+B\pa_y+C
%% where A, B and C are constant coefficient matrices.
%% NB no bc's put in yet.

DO_TEST=0;
if nargin==0
  NT=10;
  DO_TEST=1;
  if 0
    C=[3  1; 0 2];
    A=0*C;
    B=0*C;
    BCS={ [],[],[],[] };
  elseif 0
    A=zeros(2,2);
    B=eye(2);
    C=eye(2);
    BCS={ zeros(2,2),zeros(2,2),eye(2),zeros(2,2) };
    BCS=BCS([1 2 4 3]);
  elseif 0
    C=[3  1; 0 2];
    A=0*C;
    B=C;
    BCS={ zeros(2,2),zeros(2,2),eye(2),zeros(2,2) };
  elseif 0
    A=2*eye(2);
    B=3*eye(2);
    C=1*eye(2);
    z22=zeros(2,2);
    I2=eye(2);
    BCS={ zeros(4,2),zeros(4,2),[1 0;0 0;0 1;0 0],[0 0;1 0; 0 0;0 1] };
  elseif 1
    A=[4 1; 6 8];
    B=[2 .5;1/3 6];
    C=[3  1; 0 2];
    z42=zeros(4,2);
    BCS={ z42,z42,eye(4,2),fliplr(flipud(eye(4,2))) };
    %%
    B2=[A,0*A;0*A,eye(2)];
    A2=0*B2;
    C2=[B,C;-eye(2),0*C];
    z44=zeros(4,4);
    BCS2={ z44,z44,[BCS{1},BCS{3}],[BCS{2},BCS{4}] };
  end
end

N=size(A,1);
%% Minprods={tt,uu,ut,cc,cu,ct};
%% T_n-> n*U_{n-1} -> 2*n*C^{(2)}_{n-2}
J0=(1:NT);
ddtt=[1,1./J0];

if 0%% SLOW BUT SURE WAY:
  Nint=NT+5;
  [tT,wT]=OP_numint_chebyshev(Nint);
  %%
  [ipT,hnT,Tvals]=...
    OP_inprod_chebyshev(tT,wT,NT);
  alpU=1;
  [ipU,hnU,Uvals]=...
    OP_inprod_gegenbauer(tT,wT.*(1-tT.^2),alpU,NT-1);
  alpC=2;
  [ipC,hnC,Cvals]=...
    OP_inprod_gegenbauer(tT,wT.*(1-tT.^2).^2,alpC,NT-2);
  Minprods={diag(ddtt),eye(NT),ipU*Tvals*diag(ddtt),...
             2*eye(NT-1),ipC*Uvals,ipC*Tvals*diag(ddtt)};
else%% GET INNER PRODUCT METHODS ANALYTICALLY:
  Minprods={diag(ddtt),[],[],[],[],[]};
  %%
  chqA=sum(abs(A),2);
  if ~isempty(find(chqA))
    %% if A has any non-zero rows will need some
    %% C_m inner products:
    Minprods{4}=2*eye(NT-1);
    Mcu=diag(1./J0(1:end-1),-1);
    Mcu(1,:)=[];
    Mcu(1:NT-2,3:end)=Mcu(1:NT-2,3:end)-diag(1./(3:NT));
    Minprods{5}=Mcu;
    %%
    ddct=.5*ddtt(1:end-2)./J0(1:end-1);
    ddct(1)=1;
    ddct2=.5*ddtt(5:end)./J0(3:end-1);
    J1=J0(2:end);
    ddct1=-ddtt(3:end).*J1./(J1.^2-1);
    Mct=diag(ddct,-2);
    Mct(1:2,:)=[];
    Mct(1:NT-3,5:end)=Mct(1:NT-3,5:end)+diag(ddct2);
    Mct(:,3:end)=Mct(:,3:end)+diag(ddct1);
    Minprods{6}=Mct;
  end

  chqB=sum(abs(B),2);
  for j=1:N
    if chqA(j)==0 & chqB(j)~=0
      %% if A has any zero rows may need some
      %% U_m inner products:
      Minprods{2}=eye(NT);
      ddut=ddtt(1:end-1)/2;
      ddut(1)=1;
      ddut2=ddtt(3:end)/2;
      Mut=diag(ddut,-1);
      Mut(1,:)=[];
      Mut(1:NT-1,3:end)=Mut(1:NT-1,3:end)-diag(ddut2);
      Minprods{3}=Mut;
      break;
    end
  end
end

%  for j=1:6
%    Minprods{j},Minprods0{j},pause
%  end,return

Matrix_left=zeros(N*(NT+1),N*(NT+1));
Matrix_right=Matrix_left;
J_discard=[];

%% SOLVE GOVERNING EQUATIONS:
for j=1:N
  if norm(A(j,:))~=0%% 2nd order ode:

    JL=(1:NT-1)+(j-1)*(NT+1);
    JL_all=(1:NT+1)+(j-1)*(NT+1);
    J_discard=[J_discard, (NT:NT+1)+(j-1)*(NT+1)];
    Matrix_right(JL,JL_all)=Minprods{6};

    for r=1:N
      JR=(3:NT+1)+(r-1)*(NT+1);
      Matrix_left(JL,JR)=A(j,r)*Minprods{4};
      %%
      JR=(2:NT+1)+(r-1)*(NT+1);
      Matrix_left(JL,JR)=Matrix_left(JL,JR) +...
                           + B(j,r)*Minprods{5};
      %%
      JR=(1:NT+1)+(r-1)*(NT+1);
      Matrix_left(JL,JR)=Matrix_left(JL,JR) + ...
                           + C(j,r)*Minprods{6};
    end

  elseif norm(B(j,:))~=0%% 1st order ode:

    JL=(1:NT)+(j-1)*(NT+1);
    JL_all=(1:NT+1)+(j-1)*(NT+1);
    J_discard=[J_discard, (NT+1)+(j-1)*(NT+1)];
    Matrix_right(JL,JL_all)=Minprods{3};

    for r=1:N
      JR=(2:NT+1)+(r-1)*(NT+1);
      Matrix_left(JL,JR)=B(j,r)*Minprods{2};
      %%
      JR=(1:NT+1)+(r-1)*(NT+1);
      Matrix_left(JL,JR)=Matrix_left(JL,JR) +...
                           + C(j,r)*Minprods{3};
    end

  else%% algebraic eqn:

    JL=(1:NT+1)+(j-1)*(NT+1);
    Matrix_right(JL,JL)=Minprods{1};

    for r=1:N
      JR=(1:NT+1)+(r-1)*(NT+1);
      Matrix_left(JL,JR)=C(j,r)*Minprods{1};
    end

  end
end

Matrix_left_bc=Matrix_left(J_discard,:);

%% APPLY BOUNDARY CONDITIONS:
T1=ddtt;
Tm1=T1;
Tm1(2:2:end)=-Tm1(2:2:end);
%%
alpU=1;
dT1=[0,OP_interp_gegenbauer(1,alpU,{NT-1})];
dTm1=dT1;
dTm1(3:2:end)=-dTm1(3:2:end);
%%
bc1=BCS{1};%% coeffiencts for \bfy'(-1)
bc2=BCS{2};%% coeffiencts for \bfy'(1)
bc3=BCS{3};%% coeffiencts for \bfy(-1)
bc4=BCS{4};%% coeffiencts for \bfy(1)
Nbc=size(bc1,1);%% no of boundary conditions to apply
for j=1:Nbc
  for r=1:N
    JR=(1:NT+1)+(r-1)*(NT+1);
    Matrix_left_bc(j,JR) = bc1(j,r)*dTm1 +...
                            + bc2(j,r)*dT1 +...
                            + bc3(j,r)*Tm1 +...
                            + bc4(j,r)*T1;
  end
end
Matrix_left(J_discard,:)=Matrix_left_bc;

if DO_TEST
  inv00=inv(D00);
  if 0
    fT=cos(3*tT);
    gT=2*cos(1.6*tT);
    fn=inv00*ipT*fT;
    gn=inv00*ipT*fT;
%      plot(tT,fT), hold on;
%      plot(tT,Tvals*(D00*fn),'--r'), hold off;
    yn=gn/2;
    xn=fn/3-gn/6;
    xy=Matrix_left\Matrix_right*[fn;gn];
    [xy,[xn;yn]]
  elseif 0
    xn=inv00*ipT*( tT.*sin(tT) );%% x(-1)=x(1)=sin(1)
    fn=inv00*ipT*( tT.*(cos(tT)+sin(tT))+sin(tT) );%% x'(t)=f(t)
    %%
%      Matrix_left
%      Matrix_right
    RHS=Matrix_right*[fn;fn];
    RHS(J_discard)=sin(1);
    xx=Matrix_left\RHS;
    [xx,[xn;xn]]
  elseif 0
    xn=inv00*ipT*( tT.*sin(tT) );%% x(-1)=x(1)=sin(1)
    ftilde=tT.*(cos(tT)+sin(tT))+sin(tT);
    xm1=sin(1);
    %%
%      yn=inv00*ipT*( tT.^2./(tT.^2+5) );%% y(-1)=y(1)=1/6;
%      gtilde=(-tT.^4+2*tT.^3-5*tT.^2+10*tT)./(5+tT.^2).^2;
%      ym1=1/6;
    yn=inv00*ipT*( tT.*cos(tT) );%% y(-1)=-y(1)=-cos(1)
    ym1=-cos(1);
    gtilde=tT.*(cos(tT)-sin(tT))+cos(tT);
%      yn=xn;
%      gtilde=ftilde;
%      ym1=xm1;
    %%
    fn=inv00*ipT*(3*ftilde+gtilde);
    gn=inv00*ipT*(2*gtilde);
    %%
    RHS=Matrix_right*[fn;gn];
    RHS(J_discard)=[xm1;ym1];
    xy=Matrix_left\RHS;
    [xy,[xn;yn]]
  elseif 0
    xn=inv00*ipT*( tT.*cos(tT) );
    xm1=-cos(1);
    xp1=cos(1);
    ff=(-2*sin(tT)-tT.*cos(tT))*2+tT.*cos(tT)+...
        3*( cos(tT)-tT.*sin(tT) );
    fn=inv00*ipT*(ff);
    RHS=Matrix_right*[fn;fn]
    RHS(J_discard)=[xm1;xp1;xm1;xp1];
    xy=Matrix_left\RHS;
    [xy,[xn;xn]]
    for j=1:2
      for r=1:2
        [j r]
        Matrix_left((1:NT+1)+(j-1)*(NT+1),(1:NT+1)+(r-1)*(NT+1)),pause
      end
    end

  elseif 1
    ff=tT.*(cos(tT)+sin(tT))+sin(tT);
    gg=tT.*(cos(tT)-sin(tT))+cos(tT);
    fn=inv00*ipT*ff;
    gn=inv00*ipT*gg;
    %%
    RHS=Matrix_right*[fn;gn];
    xy=Matrix_left\RHS;
    %%
    [ML2,MR2]=OPcheby_ODE_matrix_constcoeff(...
       A2,B2,C2,BCS2,NT);
    RHS2=MR2*[fn;gn;0*fn;0*fn];
%      for j=1:4
%        for r=1:4
%          [j r]
%          ML2((1:NT+1)+(j-1)*(NT+1),(1:NT+1)+(r-1)*(NT+1)),pause
%        end
%      end

    xy2=ML2\RHS2;
    [xy,xy2(2*NT+3:end)]
  end
end