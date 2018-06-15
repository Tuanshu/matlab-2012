clear all

%cd('D:\110516\');

axial_res_front=9.2;      %micron, FWHM of gaussian
axial_res_back=9.5;      %micron, FWHM of gaussian

SNR_front=74;
SNR_back=74-40;      %micron, FWHM of gaussian

%Corresponding_spe_BW=importdata('FWHM.txt');      %micron, FWHM of gaussian

center_wavelength=0.84; %micron

second_sidelope=center_wavelength/2;

SNR=(30:80)';

ratio_SNR=1-10.^(-1*SNR/20);

x=-20:0.0002:20;


sig_front=axial_res_front/2/(2*log(2))^0.5;

sig_back=axial_res_front/2/(2*log(2))^0.5;
        
ratio_SNR=1-10.^(-1*SNR/20);

G_profile_front=gaussmf(x,[sig_front 0]);
G_profile_back=gaussmf(x,[sig_back 0]);

for p=1:length(SNR)
    deviation_lateral_postion_noise_level_front(p)=abs(x(find(G_profile_front>ratio_SNR(p),1,'first')));     %SD, ~Rq
    deviation_lateral_postion_noise_level_back(p)=abs(x(find(G_profile_back>ratio_SNR(p),1,'first')));     %SD, ~Rq
end
for p=1:length(SNR)
    deviation_lateral_postion_noise_level_front_2D(:,p)=deviation_lateral_postion_noise_level_front;
    deviation_lateral_postion_noise_level_back_2D(p,:)=deviation_lateral_postion_noise_level_back;
end
deviation_lateral_postion_noise_level_total=deviation_lateral_postion_noise_level_front_2D+deviation_lateral_postion_noise_level_back_2D;

imagesc(deviation_lateral_postion_noise_level_total,'xdata',SNR,'ydata',SNR);
xlabel('SNR (Front)');
ylabel('SNR (Back)');

Error=deviation_lateral_postion_noise_level_total(SNR==SNR_back,SNR==SNR_front);
%legend('5','15','25','35',length(Corresponding_spe_res));


%plot(x,G_profile);