function date_str  = GEN_get_datestring(day,year,DAY1_IS_ZERO)
%% CALL: date_str  = GEN_get_datestring(day,year,DAY1_IS_ZERO)
%% if 'day' is the day number in 'year',
%% 'date_str' is the date, in a string of the form yyyymmdd;
%% DAY1_IS_ZERO = 1 means first day is day 0;
%% DAY1_IS_ZERO = 0 means first day is day 1;
%% check against Francois script:
%% - if DAY1_IS_ZERO then should agree with;
%% ~fanf/bin/jultodate day refyear 1 0

%%convert day into a date which is a string;

if ~exist('DAY1_IS_ZERO')
   DAY1_IS_ZERO   = 0;
end
if DAY1_IS_ZERO
%  day
%  day   = day-1
   day   = day+1;
end

[days_in_year,LYEAR] = cat_year(year);
while day>days_in_year
   day                  = day-days_in_year;
   year                 = year+1;
   [days_in_year,LYEAR] = cat_year(year);
end

if mod(year,100)==0
   LYEAR = (mod(year,400)==0);
else
   LYEAR = (mod(year,4)==0);
end

months         = [31 28 31 30 31 30 31 31 30 31 30 31];%%days in each month;
if LYEAR%%add another day to February in leap years;
   months(2)      = 29;
end


cutoffs  = cumsum(months);
%%
jj    = find(day<=cutoffs);
mon   = jj(1);
if day<=31
   day2  = day;
else
   day2  = day-cutoffs(mon-1);
end
%%
if day2<10
   day_str  = ['0',num2str(day2)];
else
   day_str  = num2str(day2);
end
if mon<10
   mon_str  = ['0',num2str(mon)];
else
   mon_str  = num2str(mon);
end
date_str = [num2str(year),mon_str,day_str];

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
