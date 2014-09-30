function y=GEN_run_sum_matlab(x)

N=size(x,1);
S=size(x,2);
y=zeros(N+1,S);

for j=1:N
  y(j+1,:)=y(j,:)+x(j,:);
end