function day_out  = GEN_get_daynumber(cdate_in,DAY1_IS_ZERO,refyear)
%% CALL: day  = GEN_get_daynumber(cdate,DAY1_IS_ZERO,refyear)
%% cdate can be a vector, with rows being
%% dates of the form yyyymmdd, which may be numbers or strings,
%% day is the day number in that year;
%% if DAY1_IS_ZERO=1, 1 Jan = day 0
%% can add refyear as a reference from which to count;

if ~exist('DAY1_IS_ZERO')
   DAY1_IS_ZERO   = 0;
end

N        = size(cdate_in,1);
day_out  = zeros(N,1);

for n=1:N
   cdate = cdate_in(n,:);

   %DAY1_IS_ZERO
   if isstr(cdate)
      string   = cdate;
   else
      string   = num2str(cdate);
   end
   year  = GEN_str2num(string(1:end-4));
   mon   = GEN_str2num(string(end-3:end-2));
   day   = GEN_str2num(string(end-1:end));
   if n==1
      if ~exist('refyear')
         refyear  = year;
      end
   end

   if mon==1
      if DAY1_IS_ZERO
         day   = day-1;
      end
      day   = GEN_shiftday(day,year,refyear);
   else
      %%convert day into a date which is a string;
      months   = [31 28 31 30 31 30 31 31 30 31 30 31];%%days in each month;
      if mod(year,100)==0
         LYEAR = (mod(year,400)==0);
      else
         LYEAR = (mod(year,4)==0);
      end
      if LYEAR%%add another day to February in leap years;
         months(2)   = 29;
      end

      cutoffs  = cumsum(months);
      day      = day+cutoffs(mon-1);

      if DAY1_IS_ZERO
         day   = day-1;
      end
      day   = GEN_shiftday(day,year,refyear);
   end
   day_out(n)  = day;
end
