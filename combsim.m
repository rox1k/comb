function combsim(ModeList, SweepSpeed, ...
    Initial, Detuning, Power, Linewidth, Coupling, ...
    RefIndex, NonlinIndex, Veff, FileName)

% Mode List: List of resonance frequencies
% StartTime: 
% EndTime:
% Initial: 
% Detuning: [startDetuning endDetuning] in units of cavity linewidth
% Power: [startPower endPower] in units of Watt
% Linewidths: cavity linewidths in units of Hz
% PumpFreq: optical frequency of the pump laser in units of Hz
% Coupling: coupling strength, 1/2 for critical coupling
% RefIndex: refractive index
% NonlinIndex: nonlinear refractive index in units of m^2/W
% Veff: nonlinear mode volume
%%note: 
% (1) detuning and power will ramp linearly with time.
% (2) blue detuning corresponds to a negative detuning.

tic

global d1;
global f_start;
global f_end;
global g;
global kappa;
global omega0;
global done;
global nModes;
global detuning_start;
global detuning_end;
global start_time;
global end_time;
global omega;
global fsr;

% constants & parameters
hbar = 1.05457148e-34;
c = 299792458;
pi = 3.14159;

%read resonances

omega = 2*pi*ModeList; % resonance frequencies 
nModes = size(omega,1); % number of considered modes
% d1= (1-1e-6)*0.5*abs(omega(round(nModes/2)+1)-omega(round(nModes/2)-1)); % D1 coefficient that is 2*pi*FSR
d1 = 2*pi*fsr;
kappa = Linewidth*2*pi; % cavity decay rate (full)
omega0 = omega(round(nModes/2)); % central pumped frequency
eta = Coupling; % coupling coefficient
n0 = RefIndex; % refractive index
n2 = NonlinIndex; % nonlinear refractive index
g = hbar*omega0^2*c*n2/n0^2/Veff; % nonlinear coupling coefficient
s_start = sqrt(Power(1)/hbar/omega0); 
s_end = sqrt(Power(2)/hbar/omega0);
f_start = sqrt(8*eta*g/kappa^2)*s_start; % normalized amplitutde of input field in the beggining
f_end = sqrt(8*eta*g/kappa^2)*s_end; % normalized amplitutde of input field in the end
detuning_start = Detuning(1); % starting detuning
detuning_end = Detuning(2); % end detuning
start_time = 0; % start time
end_time = (detuning_end-detuning_start)/SweepSpeed+start_time; % end time

ode_options = ''; % default options for solver

if Initial == 0
%     initial_conditions = zeros(nModes,1);
    initial_conditions =sqrt(2*g/kappa)*0.5*(randn(1,nModes)+1i*randn(1,nModes));
else
    initial_conditions = Initial;
end


% simulating
done = 0;

PNum = 2048; % number of time points that will be use for calculations. 2048 is OK.
inc = initial_conditions;
Y = zeros(PNum,nModes);
T = zeros(PNum,1);

% start calculation with equation from
% coupledModeEquations function. We divide time interval on PNum cutoffs
% and solve equations on each interval consecutive.
% total_steps=0;
for kk=1:PNum
    t1 = start_time + (kk-1)*(end_time - start_time)/PNum; %setting time interval
    t2 = start_time + kk*(end_time - start_time)/PNum; % setting time interval
    [Tx,Yx] = ode23(@coupledModeEquations, [t1  t2], inc, ode_options); % solving
%     total_steps=total_steps+size(Tx,1);
    Y(kk,:) = Yx(end,:); % save solution for current step
    T(kk) = Tx(end); % save time
    inc = Yx(end,:); % set new intial conditions for the next step
end
toc 

% apply log10 for plotting spectrum
spectrum_plot= abs(Y);
% cut-off multiplier
multipl=10^2;
spectrum_plot(spectrum_plot>0)=log10(multipl*spectrum_plot(spectrum_plot>0));
% filter negative specrtum after log()
spectrum_plot(spectrum_plot<0)=0;

spectrum_plot_transposed=spectrum_plot.';


% ugly way to save global params in file
global coupling;
global detune_s;
global detune_e;
global modes_number;
global pump_freq;
global d2;
global d3;
global pump_power;
global linewidth;
global refr_index;
global nonlin_index;
global mode_volume; 
global sweep_speed;

a_coupling=coupling;
a_detune_s=detune_s;
a_detune_e=detune_e;
a_modes_number=modes_number;
a_pump_freq=pump_freq;
a_fsr=fsr;
a_d2=d2;
a_d3=d3;
a_pump_power=pump_power;
a_linewidth=linewidth;
a_refr_index=refr_index;
a_nonlin_index=nonlin_index;
a_mode_volume=mode_volume;
a_sweep_speed=sweep_speed;
a_mode_list = ModeList;

save([FileName '.mat']);

end