function res = get_f(t)

global start_time;
global end_time;
global f_start;
global f_end;

res = (f_start + (f_end-f_start)*(t-start_time)/(end_time-start_time));

end

