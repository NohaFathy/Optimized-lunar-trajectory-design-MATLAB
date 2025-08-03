function [jd] = julian_day(year, month, day,UT)

j0 = 367*year- fix(7*(year + fix((month + 9)/12))/4) ...
 + fix(275*month/9) + day + 1721013.5;
jd =   j0 + UT/24;

end




