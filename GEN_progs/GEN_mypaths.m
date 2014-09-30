%% PRINT USER-DEFINED PATHS:
D=path; j=[0,find(D==':')];
r=1; crit=( D(j(r)+2:j(r)+4)~='usr' );

if crit==0
	disp('There are no user-defined path directories.');
else
	disp('User-defined path directories:');
	while crit
		jj=j(r)+1:j(r+1)-1;
		disp(D(jj));
		r=r+1; crit=( D(j(r)+2:j(r)+4)~='usr' );
	end
end
                                                                          
clear;
