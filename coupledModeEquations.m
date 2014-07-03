function da = coupledModeEquations(t,a)

global kappa;
global omega;
global omega0;
global d1;
global start_time;
global end_time;
global done;
global nModes;

da = zeros(nModes,1);
detuning = get_detuning(t);

% coupled mode equations
% nonlinear term
fa = fft(a); 
fNL = fa.*conj(fa).*fa;
NL = ifft(fNL);

for k = 1:nModes
      da(k) =  ...
        -(1+1i/kappa*(omega(k)-omega0+(detuning*kappa)-(k-round(nModes/2))*d1))*a(k) ...                
        + 1i*NL(k);
%         + sqrt(2*g/kappa)*0.5*(randn(1)+1i*randn(1));
end

% TODO: redo
% adding pump to the central mode (pumped mode)
da(round(nModes/2)) = da(round(nModes/2)) + get_f(t);

% estimate elapsed progress
done_tmp = 100*(t-start_time)/(end_time-start_time);
global progress;
if done_tmp - done > 1
    done = 100*(t-start_time)/(end_time-start_time);
    set(progress,'String',[num2str(round(done)+1) ' %']);
    pause(.01);
%     display([num2str(round(done)) ' %'])
end
end