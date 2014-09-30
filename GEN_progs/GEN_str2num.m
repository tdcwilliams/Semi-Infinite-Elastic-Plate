function Y  = GEN_str2num(X)
%% CALL: Y  = GEN_str2num(X)
%% Y is a number corresponding to the string X,
%% which is a visual representation of a number;
%% ie. X=num2str(Y);

N  = size(X,1);
Y  = zeros(N,1);
%%
for j=1:N
   x  = X(j,:);
   %%
   x(x==32) = [];%%get rid of spaces;
   f        = 1;
   if x(1)==45
      f     = -1;
      x(1)  = [];
   end
   z  = x-48;
   nn = 10.^(length(z)-1:-1:0)';
   %%
   Y(j)  = f*z*nn;
end
