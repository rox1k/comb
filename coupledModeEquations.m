% Copyright (c) 2013, T. Herr, I. Mirgorodskiy, G. Lihachev, M.L. Gorodetsky
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
% OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
% nonlinear terü
fa = fft(a); 
fNL = fa.*conj(fa).*fa;
NL = ifft(fNL);

for k = 1:nModes
      da(k) =  ...
        -(1+1i/kappa*(omega(k)-omega0+(detuning*kappa)-(k-round(nModes/2))*d1))*a(k) ...                
        + 1i*NL(k);
%         + sqrt(2*g/kappa)*0.5*(randn(1)+1i*randn(1));
end

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