function simulate()

global done;
global detuning_profile;
global sweep_speed;
global timesteps_cme;
global initial_conditions;
global modes_number;
global filename;

tic
start_time = 0;
% TODO: not correct for nonlinear detuning
end_time = (detuning_profile(end)-detuning_profile(1))/sweep_speed+start_time;
done = 0;
timepoints = linspace(start_time,end_time,timesteps_cme);
[~,Y] = ode23(@coupledModeEquations,timepoints,initial_conditions,'');
toc 

Y_wo_pump = Y;
Y_wo_pump(:,round(modes_number/2)) = zeros(size(Y_wo_pump,1),1);
% dB = 10*log10(abs(Y_wo_pump).^2);
dB = 10*log10(abs(Y).^2);
dB = dB-max(max(dB))+100;
spectrum_plot=dB;
% spectrum_plot=abs(Y_wo_pump);

% apply log10 for plotting spectrum
% spectrum_plot= abs(Y);
% % % cut-off multiplier
% multipl=10^2;
% spectrum_plot(spectrum_plot>0)=log10(multipl*spectrum_plot(spectrum_plot>0));
% % % filter negative specrtum after log()
% spectrum_plot(spectrum_plot<0)=0;
spectrum_plot_transposed=spectrum_plot.';

% ugly way to save global variables in file
globals = who('global');
for kk = 1:numel(globals)
  eval(sprintf('global %s', globals{kk}));
end
save([filename '.mat']);

end