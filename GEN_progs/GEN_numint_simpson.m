function [w,x]=GEN_numint_simpson(a,numpts)
%% CALL: [w,x]=GEN_numint_simpson(a,numpts)
%% \int_0^a.f(\xi).d\xi ~ w'*f(x),
%% where x=(0:numpts-1)'*a/(numpts-1);
%% NB numpts must be odd.
%% Result is exact if f is a cubic.

w=zeros(numpts,1);
Del=a/(numpts-1);
w(1)=Del/3;
w(numpts)=Del/3;
w(2:2:end-1)=4*Del/3;
w(3:2:end-1)=2*Del/3;
%%
if nargout==2
	x=(0:numpts-1)'*Del;
end
