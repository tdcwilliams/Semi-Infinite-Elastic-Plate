function cc=GEN_mult_poly(aa,bb)

aa=fliplr(aa);
bb=fliplr(bb);
M=length(aa)-1;
N=length(bb)-1;
cc=zeros(1,M+N+1);

%% pad aa & bb with extra zeros to avoid errors:
cc(1:M+1)=aa;
aa=cc;
cc=zeros(1,M+N+1);
%%
cc(1:N+1)=bb;
bb=cc;
cc=zeros(1,M+N+1);

for j=0:M+N
  for k=max(0,j-N):j
    cc(j+1)=cc(j+1)+aa(k+1)*bb(1+j-k);
  end
end
cc=fliplr(cc);