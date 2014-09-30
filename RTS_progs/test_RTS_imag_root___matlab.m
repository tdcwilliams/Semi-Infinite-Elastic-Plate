function test_RTS_imag_root___matlab(pram,H,ice_or_water)

if ice_or_water==1%%
  %% del=pram;
  fxn=@RTS_imag_root_ice_matlab;
  fxn0=@RTS_imag_root_ice;
else
  %% lam=pram;
  fxn=@RTS_imag_root_wtr_matlab;
  fxn0=@RTS_imag_root_wtr;
end

j=1;
i1=j*pi;
i0=i1-pi;
w=feval(fxn,pram,H,i0,i1);
w0=feval(fxn0,pram,H,i0,i1);
tst=[w w0]