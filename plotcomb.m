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

function plotcomb(filename,Snapshot,flag)

c = 299792458;
file = load(filename);
reslist=file.a_mode_list;
detune_start=file.a_detune_s;
detune_end=file.a_detune_e;

if strfind(filename, 'coupledeq')==1
    % get parameters from coupled eq file
    nModes = file.nModes; % number of modes
    Y = file.Y; % field amplitudes for every mode
    T = file.T; % time
    det=linspace(detune_start,detune_end,size(T,1));
    fi=linspace(-pi,pi,size(Y,2));    

    N = (1-round(nModes/2):1:round(nModes/2)-1); % mode numbers
%     w = file.omega; % frequency corresponding to each mode
    spectrum_plot=file.spectrum_plot;
    spectrum_plot_transposed=file.spectrum_plot_transposed;
    %det = get_detuning(T); % detuning for given time
    int_pow_a = sum(Y.*conj(Y),2); % intracavity power by summing squared absolute values of field ampl. over each raw
%     fi = pi*N/(round(nModes/2)-1); % angular coordinate.
    % number of points in the range of [-pi;pi] or [0;2pi]
    max_amp=max(int_pow_a);
    ind = find(det>=Snapshot,1);
    size_fi = max(size(fi));
    pulse = zeros(1,size_fi);    
    for jj=1:size_fi
        pulse(1,jj) = sum(Y(ind,:).*exp(-1i*N*fi(jj)),2);
    end
else
    % get parameters from llessfm file
    E=file.E;
%     detune_step=(detune_end-detune_start)/size(E,1);
%     det=detune_start:detune_step:detune_end;
    det=linspace(detune_start,detune_end,size(E,1));
    ind=find(det>=Snapshot,1);
    power=(abs(E(ind,:)).^2);
    int_field = sum((abs(E).^2),2);
    max_amp=max(int_field);
    fi=linspace(-pi,pi,size(E,2));
end

%plotting
switch flag
    case 'all'
        % plot for coupled equation data
        if strfind(filename, 'coupledeq')==1
            subplot(2,2,1:2)
            set(plot(det,int_pow_a,'m.'),'MarkerSize',4)
            xlabel('\zeta_0')
            ylabel('$\sum |a_i|^2$','interpreter','latex')
            grid on
            line([Snapshot Snapshot],[0 max_amp], 'Color','b','LineStyle','--')
            
            %subplot(3,2,3:4)
            %contourf(det,N,spectrum_plot_transposed);
            %shading flat;
            %xlabel('time')
            %ylabel('mode number')
            %colormap('hsv')
            
            subplot(2,2,3)
            for ii=1:1:length(reslist)
                x_bar(ii)=(10^(9)*c)/reslist(ii);
            end
            bar(x_bar,spectrum_plot(ind,:),'r');
            % ylabel('$ |a_i|$','interpreter','latex')
            xlabel('$\lambda$, nm','interpreter','latex');
            ylabel('amplitude, a.u.');
            title ('Optical spectrum');
            % set(gca, 'FontSize', 12);
            
            subplot(2,2,4)
            plot(fi,abs(fftshift(pulse)));
            xlabel('$\phi$','interpreter','latex')
            ylabel('$\Psi$','interpreter','latex')
            title ('Waveform');
        else
            subplot(2,2,1:2)
            plot(det,int_field);
            xlabel ('Detuning, (halflinewidths)');
            ylabel ('Intracavity power (a.u.)');
            grid on
            line([Snapshot Snapshot],[0 max_amp], 'Color','b','LineStyle','--')
                             
            subplot(2,2,3)            
            for ii=1:1:length(reslist)
                x_bar(ii)=(10^(9)*c)/reslist(ii);
            end
            bar(x_bar,abs(fftshift(fft(power))),'r');
            xlabel('$\lambda$, nm','interpreter','latex');
            ylabel('amplitude, a.u.');
            title ('Optical spectrum');
                        
            subplot(2,2,4)
            plot(fi,power);  % plotting power in time domain in resonator
            xlabel('$\phi$','interpreter','latex')
            ylabel('$\Psi$','interpreter','latex')
            title ('Waveform');
        end
        % plot in separate figures
    case 'total_field'        
        figure
        if strfind(filename, 'coupledeq')==1
            set(plot(det,int_pow_a,'m.'),'MarkerSize',4)
            xlabel('\zeta_0')
            ylabel('$\sum |a_i|^2$','interpreter','latex')
            grid on
        else
            plot(det,int_field);
            xlabel ('detuning, un. of linewidth');
            ylabel ('intracavity power, a.u.');
            grid on
        end
    case 'amps'
        figure
        if strfind(filename, 'coupledeq')==1
            contourf(det,N,spectrum_plot_transposed);
            shading flat;
            xlabel('time')
            ylabel('mode number')
        end         
    case 'spectrum'
        figure
        for ii=1:1:length(reslist)
            x_bar(ii)=(10^(9)*c)/reslist(ii);
        end
        if strfind(filename, 'coupledeq')==1
            bar(x_bar,spectrum_plot(ind,:),'r');
            xlabel('$\lambda$, nm','interpreter','latex');
            ylabel('amplitude, a.u.');
            title ('Optical spectrum');            
        else
            bar(x_bar,abs(fftshift(fft(power))),'r');
            xlabel('$\lambda$, nm','interpreter','latex');
            ylabel('amplitude, a.u.');
            title ('Optical spectrum');
        end
    case 'waveform'
        figure
        if strfind(filename, 'coupledeq')==1
            plot(fi,abs(fftshift(pulse)));
            xlabel('$\phi$','interpreter','latex')
            ylabel('$\Psi$','interpreter','latex')
            title ('Waveform');
        else
            plot(fi,power);
            xlabel('$\phi$','interpreter','latex')
            ylabel('$\Psi$','interpreter','latex')
            title ('Waveform');
        end
end
end