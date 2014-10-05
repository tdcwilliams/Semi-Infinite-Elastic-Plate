function [kinf,jout] = GEN_sort_iceroots_inf(rts);

jneg  = find(rts<0);%%rh half-plane
kneg  = rts(jneg);
if imag(kneg(1))<0
   kneg  = flipud(kneg);
end
%%
jpos  = find(rts>0);%%rh half-plane
kpos  = rts(jpos);
%%
jin   = {[2 3],[1 3],[1 2]};
for j=1:3
   k0    = kpos(j);
   kk    = kpos(jin{j});
   d0(j) = min(abs( k0-conj(kk) ));
end

jmax  = find(d0==max(d0));
k0    = kpos(jmax);
kcx   = kpos(jin{jmax});

%%k0 is real root; sort kcx
if imag(kcx(1))<0
   kcx   = flipud(kcx);
end
kinf  = [k0;kcx;kneg];

if nargout==2
   for j=1:5
      %jout(j)  = find(rts(j)==kinf);
      jout(j)  = find(rts==kinf(j));
   end
end

kinf(1)  = real(kinf(1));
