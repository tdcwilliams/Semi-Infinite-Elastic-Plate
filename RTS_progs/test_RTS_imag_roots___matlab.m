function test_RTS_imag_roots_matlab(pram,H,ice_or_water)

if ice_or_water==1%%
  %% del=pram;
  N0=1;
  N1=10;
  fxn=@RTS_imag_roots_ice_matlab;
  fxn0=@RTS_imag_roots_ice;
  w=feval(fxn,pram,H,N0,N1);
  w0=feval(fxn0,pram,H,N0,N1);
else
  %% lam=pram;
  N1=10;
  fxn=@RTS_imag_roots_wtr_matlab;
  fxn0=@RTS_imag_roots_wtr;
  w=feval(fxn,pram,H,N1);
  w0=feval(fxn0,pram,H,N1);
end

tst=[w w0]