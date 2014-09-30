function [tt,ww]=OP_numint_logjacobi(alf,bet,Nint)
%% CALL: [tt,ww]=OP_numint_logjacobi(alf,bet,Nint)
%% approximates integral
%%  \int_{-1}^1.w(t).f(t)~\sum_{n=1}^Nint.[ww]_n*f([tt]_n),
%%    where w(t)=(1+t)^beta*(1-t)^alpha*log(1+t)
%% makes a substitution u=log(2)-log(1+t)
%%  to change it into a Gauss-Laguerre type integral.

[vv,wg_lag]=OP_numint_laguerre(alf,Nint);
uu=vv/(1+bet);
tt=2*exp(-uu)-1;
%%
w0=(log(2)-uu).*( (1-exp(-uu))./uu ).^alf;
ww=2^(alf+bet+1)/(1+bet)^(1+alf)*wg_lag.*w0;