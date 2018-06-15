clear all

cd('D:\110516\');

axial_res=importdata('FWHM_PSF.txt');      %micron, FWHM of gaussian

Corresponding_spe_BW=importdata('FWHM.txt');      %micron, FWHM of gaussian

center_wavelength=0.78; %micron

second_sidelope=center_wavelength/2;

SNR=20:2:100;

SNR=SNR';

ratio_SNR=1-10.^(-1*SNR/20);

x=-10:0.0002:10;

for j=1:length(axial_res)

        

    sig=axial_res(j)/2/(2*log(2))^0.5;
        
        ratio_SNR=1-10.^(-1*SNR/20)/0.578;

    G_profile=gaussmf(x,[sig 0]);

    for p=1:length(SNR)
        deviation_lateral_postion_noise_level(p)=abs(x(find(G_profile>ratio_SNR(p),1,'first')));     %SD, ~Rq
    end

    deviation_lateral_postion_noise_level_record(:,j)=deviation_lateral_postion_noise_level';

end
imagesc(deviation_lateral_postion_noise_level_record,'xdata',Corresponding_spe_BW,'ydata',SNR);

xlabel('Bandwidth of Light Source (nm)');
ylabel('Signal to Noise Ratio (dB)');
%legend('5','15','25','35',length(Corresponding_spe_res));
plot(SNR,deviation_lateral_postion_noise_level_record);


%plot(x,G_profile);