clear all

cd('D:\110516\');

axial_res=importdata('FWHM_PSF.txt');      %micron, FWHM of gaussian

Corresponding_spe_BW=importdata('FWHM.txt');      %micron, FWHM of gaussian

center_wavelength=0.976; %micron

second_sidelope=center_wavelength/2;

SNR=10:5:100;

SNR=SNR';

ratio_SNR=1-10.^(-1*SNR/20);

x=-10:0.001:10;

for p=1:length(axial_res)
sig=axial_res(p)/2/(2*log(2))^0.5;

G_profile=gaussmf(x,[sig 0]);

required_ratio=G_profile(find(x>second_sidelope,1,'first'));

required_SNR(p)=20*log10(1/(1-required_ratio));
end
plot(Corresponding_spe_BW,required_SNR);


xlabel('Bandwidth of Light Source (nm)');
ylabel('SNR_c_r_i_t_i_c_a_l (dB)');

%plot(x,G_profile);