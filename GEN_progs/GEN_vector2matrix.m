function M=GEN_vector2matrix(v)

jp=2:length(v);
M=WT_rearr(v(jp),v(jp-1),1);
M=M+i*WT_rearr(v(jp)/i,v(jp-1)/i,1);
M(:,1)=v;
