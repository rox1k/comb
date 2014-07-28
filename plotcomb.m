function plotcomb(filename,Snapshot,flag)

global c;

file = load(filename);
reslist=file.a_mode_list;
detuning=file.a_detuning_profile;
pump_profile=file.a_pump_profile;
initial_conditions=file.a_initial_conditions;

% get parameters from coupled eq file
ind=round(Snapshot*length(detuning));
if ind == 0
    ind = 1;
end
modes_number = file.a_modes_number; % number of modes
Y = file.Y; % field amplitudes for every mode
Y_wo_pump = file.Y_wo_pump;
% T = file.T; % time

%     det=linspace(detune_start,detune_end,size(T,1));
fi=linspace(-pi,pi,size(Y,2));    

N = (1-round(modes_number/2):1:round(modes_number/2)-1); % mode numbers
%     w = file.omega; % frequency corresponding to each mode
spectrum_plot=file.spectrum_plot;
spectrum_plot_transposed=file.spectrum_plot_transposed;
%det = get_detuning(T); % detuning for given time
int_pow_a = sum(Y.*conj(Y),2); % intracavity power by summing squared absolute values of field ampl. over each raw
power_wo_pump = sum(abs(Y_wo_pump ).^2,2);
%     fi = pi*N/(round(nModes/2)-1); % angular coordinate.
% number of points in the range of [-pi;pi] or [0;2pi]
max_amp=max(int_pow_a);

%     ind = find(det>=Snapshot,1);
size_fi = max(size(fi));
pulse = zeros(1,size_fi);    
for jj=1:size_fi
    pulse(1,jj) = sum(Y(ind,:).*exp(-1i*double(N*fi(jj))),2);
end
%plotting
switch flag
    case 'all'
        % plot for coupled modes equation data
      
        subplot(2,2,1:2)
        zoom on
%             set(plot(detuning,int_pow_a,'m.'),'MarkerSize',4)
        set(plot(1:length(detuning),int_pow_a,'m.'),'MarkerSize',4)
        hold on
        set(plot(1:length(detuning),power_wo_pump,'c.'),'MarkerSize',4)
%             xlabel('\zeta_0')
        xlabel('time')
        ylabel('$\sum |a_i|^2$','interpreter','latex')
        grid on
        line([ind ind],[0 max_amp], 'Color','b','LineStyle','--')
        hold off;
        xlim([1,length(detuning)]);

        %subplot(3,2,3:4)
        %contourf(det,N,spectrum_plot_transposed);
        %shading flat;
        %xlabel('time')
        %ylabel('mode number')
        %colormap('hsv')

        subplot(2,2,3)
        zoom off
        for ii=1:1:length(reslist)
            x_bar(ii)=(10^(9)*c)/reslist(ii);
        end
        stem(x_bar,spectrum_plot(ind,:),'MarkerSize',2);
        xlim([min(x_bar) max(x_bar)]);
        ylim([0 110]);
%         bar(x_bar,spectrum_plot(ind,:),'r');
        % ylabel('$ |a_i|$','interpreter','latex')
        xlabel('$\lambda$, nm','interpreter','latex');
        ylabel('Intensity, dB');
        title ('Optical spectrum');
        % set(gca, 'FontSize', 12);

        subplot(2,2,4)
%         plot(fi,abs(fftshift(pulse)));
        plot(fi,abs(pulse));
        xlim([min(fi) max(fi)]);
        xlabel('$\phi$','interpreter','latex')
        ylabel('$\Psi$','interpreter','latex')
        title ('Waveform');              

    % plot in separate figures
    case 'total_field'        
        figure
        set(plot(detuning,int_pow_a,'m.'),'MarkerSize',4)
        xlabel('\zeta_0')
        ylabel('$\sum |a_i|^2$','interpreter','latex')
        grid on

    case 'amps'
        figure
        
%             contourf(detuning,N,spectrum_plot_transposed);
        pcolor(detuning,N,spectrum_plot_transposed);
        shading flat;
        xlabel('time')
        ylabel('mode number')

    case 'spectrum'
        figure
        for ii=1:1:length(reslist)
            x_bar(ii)=(10^(9)*c)/reslist(ii);
        end
        
        stem(x_bar,spectrum_plot(ind,:),'MarkerSize',2);
        xlim([min(x_bar) max(x_bar)]);
        ylim([0 110]);
%         bar(x_bar,spectrum_plot(ind,:),'r');
        xlabel('$\lambda$, nm','interpreter','latex');
        ylabel('Intensity, dB');
        title ('Optical spectrum');     
    case 'waveform'
        figure
%         plot(fi,abs(fftshift(pulse)));
        plot(fi,abs(pulse));
        xlabel('$\phi$','interpreter','latex')
        ylabel('$\Psi$','interpreter','latex')
        title ('Waveform');
        xlim([min(fi) max(fi)]);
end