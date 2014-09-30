function D=GEN_numdiff_from_right(n,Del,varargin)
%% CALL D=GEN_numdiff_from_right(n,Del);
%% if nn=(0:n)' (where n=0,1,2 or 3), x=nn*Del & df=(d/dx).^nn.f(x0),
%% then this function calculates a matrix D so that
%% df = D*f(x0+x) + O(Del) --- ie it estimates n derivatives by using
%% values of the function to the right of the point in question.
%% if a 3rd argument is entered, D is calculated so that
%% df = D*f(x0-flipud(x)) + O(Del) --- ie it estimates n derivatives by using
%% values of the function to the left of the point in question.

D=eye(n+1);
diff_from_left=(nargin==3);
m=(-1)^diff_from_left; Del=m*Del;
%%
if n>=1
	D(2,1:2)=[-1 1]/Del;
end
if n>=2
	D(3,1:3)=[1 -2 1]/Del^2;
end
if n==3
	D(4,:)=[-1 3 -3 1]/Del^3;
end
%%
if diff_from_left
	D=fliplr(D);
end
