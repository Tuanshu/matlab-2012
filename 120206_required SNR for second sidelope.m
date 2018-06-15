clear all

axial_res=1.4;      %micron, FWHM of gaussian

center_wavelength=0.56; %micron

second_sidelope=center_wavelength/2;

SNR=10:5:100;

SNR=SNR';

ratio_SNR=1-10.^(-1*SNR/20);

x=-10:0.001:10;

sig=axial_res/2/(2*log(2))^0.5;

G_profile=gaussmf(x,[sig 0]);

required_ratio=G_profile(find(x>second_sidelope,1,'first'));

required_SNR=20*log10(1/(1-required_ratio));


%plot(x,G_profile);