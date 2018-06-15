clear all

axial_res=6;      %micron, FWHM of gaussian

center_wavelength=0.56; %micron

second_sidelope=center_wavelength/2;

SNR=10:5:100;

SNR=SNR';

ratio_SNR=1-10.^(-1*SNR/20);

x=-10:0.001:10;

sig=axial_res/2/(2*log(2))^0.5;

G_profile=gaussmf(x,[sig 0]);

for j=1:length(SNR)
    deviation_lateral_postion_noise_level(j)=abs(x(find(G_profile>ratio_SNR(j),1,'first')));     %SD, ~Rq
end

deviation_lateral_postion_noise_level=deviation_lateral_postion_noise_level';

plot(SNR,deviation_lateral_postion_noise_level);

%plot(x,G_profile);