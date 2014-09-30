%%testOP_logint1_gegenbauer.m
clear;

alf=.6;
alf2=alf-.5;
tt=1/9;
tau=2*tt-1;
%%
F=inline('2*(1-s.^2).^(alpha-.5).*log(abs(t-s))')
quad(@(s)F(alf,s,tt(1)),0,1,1e-8)
%%
G=inline('(1-sig).^(alpha-.5).*( (3+sig)/4 ).^(alpha-.5).*log(abs(tau-sig)/2)')
int=quad(@(sig)G(alf,sig,tau(1)),-1,1,1e-8)
%%
G2A=inline('(1-sig).^(alpha-.5).*log(abs(tau-sig)).*( ((3+sig)/4).^(alpha-.5) - ((3+tau)/4).^(alpha-.5))')
int2A=quad(@(sig)G2A(alf,sig,tau(1)),-1,1)
%  sig=(-.9:.001:.9)';
%  plot(sig,[G0(alf,sig,tau(1)),0*sig])

G2=inline('(1-sig).^(alpha-.5).*log(abs(tau-sig)).*( ((3+sig)/4).^(alpha-.5) - ((3+tau)/4).^(alpha-.5) -.25*(alpha-.5)*(sig-tau).* ((3+tau)/4).^(alpha-1.5))')
int2=quad(@(sig)G2(alf,sig,tau(1)),-1,1)

%%
Nint=1500;
[sig,wJ]=OP_numint_jacobi(alf2,0,Nint);
ig=log(abs(tau-sig)).*( ((3+sig)/4).^alf2 -...
     ((3+tau)/4).^alf2 -.25*alf2*(sig-tau).* ((3+tau)/4).^(alf2-1));
%[size(wJ'),size(ig)]
[wJ'*ig]
%%


G1=inline('(1-sig).^(alpha-.5).*( (3+sig)/4 ).^(alpha-.5).*log(.5)')
int1=quad(@(sig)G1(alf,sig),-1,1,1e-8)

F0=inline('(1-s).^(alpha-.5).*log(abs(t-s))')
I0=quad(@(s)F0(alf,s,tau(1)),-1,1,1e-8);
int0=((3+tau(1))/4)^alf2*I0
%%
F1=inline('(1-s).^(alpha-.5).*(t-s).*log(abs(t-s))')
I1=quad(@(s)F1(alf,s,tau(1)),-1,1,1e-8);
int1B=-.25*alf2*((3+tau(1))/4)^(alf2-1)*I1
%%
[int,int0+int1+int2A,int0+int1+int2+int1B]