clear all

axial_res=[6.3 28.7];      %micron, FWHM of gaussian

center_wavelength=0.78; %micron

second_sidelope=center_wavelength/2;

SNR=50:2:100;

SNR=SNR';

ratio_SNR=1-10.^(-1*SNR/20);

x=-10:0.0002:10;

for j=1:length(axial_res)

        

    sig=axial_res(j)/2/(2*log(2))^0.5;
    
    if j==1
        
        ratio_SNR=1-10.^(-1*SNR/20);
    elseif j==2
        ratio_SNR=1-10.^(-1*SNR/20)/0.578;
    end

    G_profile=gaussmf(x,[sig 0]);

    for p=1:length(SNR)
        deviation_lateral_postion_noise_level(p)=abs(x(find(G_profile>ratio_SNR(p),1,'first')));     %SD, ~Rq
    end

    deviation_lateral_postion_noise_level_record(:,j)=deviation_lateral_postion_noise_level';

end
plot(SNR,deviation_lateral_postion_noise_level_record);

xlabel('Signal to Noise Ratio (dB)');
ylabel('Measurement Error Caused by SNR (micron)');


deviation_lateral_postion_noise_level_sum=(deviation_lateral_postion_noise_level_record(:,1).^2+deviation_lateral_postion_noise_level_record(:,2).^2).^0.5;

plot(SNR,deviation_lateral_postion_noise_level_sum);

xlabel('Signal to Noise Ratio (dB)');
ylabel('Measurement Error Caused by SNR (micron)');

%plot(x,G_profile);