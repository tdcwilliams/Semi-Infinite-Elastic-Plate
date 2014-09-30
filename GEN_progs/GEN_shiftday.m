function day   = GEN_shiftday(day,year,refyear)

if refyear>year
   disp('ERROR in GEN_shiftday.m: refyear>year');
   return;
end

day0  = 0;
for yr=refyear:1:year-1
    day0 = day0+cat_year(yr);
end
day   = day0+day;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [days_in_year,LYEAR] = cat_year(year);

if mod(year,100)==0
   LYEAR = (mod(year,400)==0);
else
   LYEAR = (mod(year,4)==0);
end
   days_in_year   = 366;
if LYEAR%%add another day to February in leap years;
   days_in_year   = 366;
else
   days_in_year   = 365;
end
