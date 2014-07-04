function combsim()

tic

global done;
global detuning_profile;
global start_time;
global end_time;
global sweep_speed;
global timesteps_cme;
global initial_conditions;
global modes_number;
global filename;

start_time = 0;
end_time = (detuning_profile(end)-detuning_profile(1))/sweep_speed+start_time;

ode_options = ''; % default options for solver

% if Initial == 0
% %     initial_conditions = zeros(nModes,1);
%     initial_conditions =sqrt(2*g/kappa)*0.5*(randn(1,nModes)+1i*randn(1,nModes));
% else
%     initial_conditions = Initial;
% end

% simulating
done = 0;

inc = initial_conditions;
Y = zeros(timesteps_cme,modes_number);
T = zeros(timesteps_cme,1);

% start calculation with equation from
% coupledModeEquations function. We divide time interval on PNum cutoffs
% and solve equations on each interval consecutive.
% total_steps=0;
for kk=1:timesteps_cme
    t1 = start_time + (kk-1)*(end_time - start_time)/timesteps_cme; %setting time interval
    t2 = start_time + kk*(end_time - start_time)/timesteps_cme; % setting time interval
    [Tx,Yx] = ode23(@coupledModeEquations, [t1  t2], inc, ode_options); % solving
%     total_steps=total_steps+size(Tx,1);
    Y(kk,:) = Yx(end,:); % save solution for current step
    T(kk) = Tx(end); % save time
    inc = Yx(end,:); % set new intial conditions for the next step
end
toc 

% apply log10 for plotting spectrum
spectrum_plot= abs(Y);
% % cut-off multiplier
% multipl=10^2;
% spectrum_plot(spectrum_plot>0)=log10(multipl*spectrum_plot(spectrum_plot>0));
% % filter negative specrtum after log()
% spectrum_plot(spectrum_plot<0)=0;

spectrum_plot_transposed=spectrum_plot.';


% ugly way to save global params in file
global coupling;
global modes_number;
global pump_freq;
global fsr;
global d2;
global d3;
global linewidth;
global refr_index;
global nonlin_index;
global mode_volume; 
global sweep_speed;
global pump_profile;
global reslist;

a_coupling=coupling;
a_detuning_profile=detuning_profile;
a_modes_number=modes_number;
a_pump_freq=pump_freq;
a_fsr=fsr;
a_d2=d2;
a_d3=d3;
a_pump_profile=pump_profile;
a_linewidth=linewidth;
a_refr_index=refr_index;
a_nonlin_index=nonlin_index;
a_mode_volume=mode_volume;
a_sweep_speed=sweep_speed;
a_mode_list = reslist;
a_initial_conditions=initial_conditions;

save([filename '.mat']);

end