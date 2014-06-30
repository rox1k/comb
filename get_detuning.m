function detuning = get_detuning(t)

global start_time;
global end_time;
global detuning_start;
global detuning_end;

detuning = (detuning_start + (detuning_end-detuning_start)*(t-start_time)/(end_time-start_time));

end


