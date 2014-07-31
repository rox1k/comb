function plotcomb(filename,Snapshot,flag)

global c;

file = load(filename);

reslist=file.reslist;
% detuning_profile=file.detuning_profile;
modes_number = file.modes_number;
Y = file.Y; % field amplitudes for every mode
Y_wo_pump = file.Y_wo_pump;
plotsteps=file.plotsteps;
spectrum_plot=file.spectrum_plot;
spectrum_plot_transposed=file.spectrum_plot_transposed;

ind=round(Snapshot*plotsteps);
if ind == 0
    ind = 1;
end

phi=linspace(-pi,pi,size(Y,2));    
N = (1-round(modes_number/2):1:round(modes_number/2)-1); % mode numbers
%det = get_detuning(T); % detuning for given time
int_pow_a = sum(Y.*conj(Y),2); % intracavity power by summing squared absolute values of phield ampl. over each raw
power_wo_pump = sum(abs(Y_wo_pump ).^2,2);
%     phi = pi*N/(round(nModes/2)-1); % angular coordinate.
max_amp=max(int_pow_a);

%     ind = find(det>=Snapshot,1);
size_phi = max(size(phi));
pulse = zeros(1,size_phi);    
for jj=1:size_phi
    pulse(1,jj) = sum(Y(ind,:).*exp(-1i*double(N*phi(jj))),2);
end
%plotting
switch flag
    case 'all'
        subplot(2,2,1:2)
        set(plot(1:plotsteps,int_pow_a,'m.'),'MarkerSize',4)
        hold on
        set(plot(1:plotsteps,power_wo_pump,'c.'),'MarkerSize',4)
%             xlabel('\zeta_0')
        xlim([1,plotsteps]);
        xlabel('time')
        ylabel('$\sum |a_i|^2$','interpreter','latex')
        grid on
        line([ind ind],[0 max_amp], 'Color','b','LineStyle','--')
        hold off;        
        zoom on;

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
        stem(x_bar,spectrum_plot(ind,:),'MarkerSize',2);
        xlim([min(x_bar) max(x_bar)]);
%         ylim([0 110]);
%         bar(x_bar,spectrum_plot(ind,:),'r');
        % ylabel('$ |a_i|$','interpreter','latex')
        xlabel('$\lambda$, nm','interpreter','latex');
        ylabel('Intensity, dB');
        title ('Optical spectrum');
        % set(gca, 'FontSize', 12);

        subplot(2,2,4)
%         plot(phi,abs(fftshift(pulse)));
        plot(phi,abs(pulse),'r');
        xlim([min(phi) max(phi)]);
        xlabel('$\phi$','interpreter','latex')
        ylabel('$\Psi$','interpreter','latex')
        title ('Waveform');              

    % plot in separate figures
    case 'total_field'        
        figure
        set(plot(1:plotsteps,int_pow_a,'m.'),'MarkerSize',4)
        hold on
        set(plot(1:plotsteps,power_wo_pump,'c.'),'MarkerSize',4)
%             xlabel('\zeta_0')
        xlim([1,plotsteps]);
        xlabel('time')
        ylabel('$\sum |a_i|^2$','interpreter','latex')
        grid on

    case 'amps'
        figure
        
%             contourf(detuning_profile,N,spectrum_plot_transposed);
        pcolor(1:plotsteps,N,spectrum_plot_transposed);
        colorbar
%         colormap('hsv');
        shading flat;
        xlabel('time')
        ylabel('mode number')

    case 'spectrum'
        figure
        for ii=1:1:length(reslist)
            x_bar(ii)=(10^(9)*c)/reslist(ii);
        end
        
        stem(x_bar,spectrum_plot(ind,:),'MarkerSize',2,'r');
        xlim([min(x_bar) max(x_bar)]);
%         ylim([0 110]);
%         bar(x_bar,spectrum_plot(ind,:),'r');
        xlabel('$\lambda$, nm','interpreter','latex');
        ylabel('Intensity, dB');
        title ('Optical spectrum');     
    case 'waveform'
        figure
%         plot(phi,abs(fftshift(pulse)));
        plot(phi,abs(pulse),'r');
        xlabel('$\phi$','interpreter','latex')
        ylabel('$\Psi$','interpreter','latex')
        title ('Waveform');
        xlim([min(phi) max(phi)]);
end