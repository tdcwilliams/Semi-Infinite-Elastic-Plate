function gam=RTS_ice_roots_shallow(del,H)

x=roots([H 0 del*H -1]);



%gam(1)=sqrt(x(find(
gam=sqrt(x);
gam=gam-2*gam.*( angle(gam)<0 );
