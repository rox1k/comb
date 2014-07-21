function da = coupledModeEquations(t,a)

global kappa;
global omega;
global omega0;
global d1;
global start_time;
global end_time;
global done;
global modes_number;
global detuning_profile;
global pump_profile;

da = zeros(modes_number,1);
detuning_indx=round(length(detuning_profile)*(t-start_time)/(end_time-start_time));
if detuning_indx == 0 
    detuning = detuning_profile(1);
else 
    detuning = detuning_profile(detuning_indx);
end
force_indx=round(length(pump_profile)*(t-start_time)/(end_time-start_time));
if force_indx == 0
    force=pump_profile(1);
else
    force=pump_profile(force_indx);
end

% coupled mode equations
% nonlinear term
fa = fft(a);
fNL = fa.*conj(fa).*fa;
NL = ifft(fNL);

for k = 1:modes_number
  da(k)=-(1+1i*2/kappa*double((omega(k)-omega0+(detuning*kappa)-(k-round(modes_number/2))*d1)))*a(k)+1i*NL(k);
end
% adding pump to the central mode (pumped mode)
da(round(modes_number/2)) = da(round(modes_number/2)) + force;

% estimate elapsed progress
done_tmp = 100*(t-start_time)/(end_time-start_time);
global progress;
if done_tmp - done > 1
    done = 100*(t-start_time)/(end_time-start_time);
    set(progress,'String',[num2str(round(done)+1) ' %']);
    pause(.01);
end
end