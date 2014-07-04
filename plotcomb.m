function plotcomb(filename,Snapshot,flag)

global c;
file = load(filename);
reslist=file.a_mode_list;
det=file.a_detuning_profile;
pump_profile=file.a_pump_profile;
initial_conditions=file.a_initial_conditions;
if strfind(filename, 'coupledeq')==1
    % get parameters from coupled eq file
    modes_number = file.a_modes_number; % number of modes
    Y = file.Y; % field amplitudes for every mode
    T = file.T; % time
    
%     det=linspace(detune_start,detune_end,size(T,1));
    fi=linspace(-pi,pi,size(Y,2));    

    N = (1-round(modes_number/2):1:round(modes_number/2)-1); % mode numbers
%     w = file.omega; % frequency corresponding to each mode
    spectrum_plot=file.spectrum_plot;
    spectrum_plot_transposed=file.spectrum_plot_transposed;
    %det = get_detuning(T); % detuning for given time
    int_pow_a = sum(Y.*conj(Y),2); % intracavity power by summing squared absolute values of field ampl. over each raw
%     fi = pi*N/(round(nModes/2)-1); % angular coordinate.
    % number of points in the range of [-pi;pi] or [0;2pi]
    max_amp=max(int_pow_a);
    % TODO: redo using indexes
    ind = find(det>=Snapshot,1);
    size_fi = max(size(fi));
    pulse = zeros(1,size_fi);    
    for jj=1:size_fi
        pulse(1,jj) = sum(Y(ind,:).*double(exp(-1i*N*fi(jj))),2);
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
            
            figure()
            plot(linspace(1,length(pump_profile),length(pump_profile)),pump_profile)
            title ('Pump profile');
            
            figure()
            plot(linspace(1,length(initial_conditions),length(initial_conditions)),abs(initial_conditions))
            title ('Initial_conditions');
            
            figure()
            plot(linspace(1,length(det),length(det)),det)
            title ('Detuning');
            
            figure()
            plot(N,reslist)
            title ('Eigen modes');
            
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