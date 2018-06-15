clear all

%% Setting


thickness_temp=2.68E-6;
n_should=1.5;

Starting_pixel_f_considered=250;
Ending_pixel_f_considered=350;

range_specified=1;

Max_Wavelength=1200;             %nm
Min_Wavelength=800;             %nm
N_f=8192*4;
N_t=N_f*8;
ROI_ratio=1/4;                  %only consider the first ROI_ratio data in TD

DC_cutoff=0;                 %in TD


%% Data Loading



cd('D:\110516\');
Data=importdata('120503_SLD Spectrum.txt');
Wavelength=Data(:,1);           %nm

C=3E8;

Frequency=C./(Wavelength*1E-9);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_New=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency_New=Frequency_New';

Position_Center_of_PSF=100E-6;  %mm

%Spectrum=Data(:,2)-Data(1,2);%-Data_R(:,2)-Data_S(:,2);
%Spectrum_Frequency=(Spectrum.*((Wavelength*1E-9).^2)/C)/max(Spectrum.*((Wavelength*1E-9).^2)/C);
%Spectrum_New=interp1(Frequency,Spectrum_Frequency,Frequency_New).*exp(i*4*pi.*Frequency_New/C.*(Position_Center_of_PSF));
%Spectrum_New(isnan(Spectrum_New))=0;
%Spectrum_New(Frequency_New<Min_Frequency)=0;
%%

FWHM=30;    %nm

Wavelength_Center=976;  %nm

Frequency_Center=C/(Wavelength_Center*1E-9);             %Hz
FWHM_Frequency=(FWHM*1E-9)*C/((Wavelength_Center*1E-9)^2);

Spectrum_New=gaussmf(Frequency_New,[FWHM_Frequency Frequency_Center]);

%plot(Frequency_New,Spectrum_New);

Wavelength_micron=C./Frequency_New*1E6;

%% To time domain

Spectrum_New((N_f+1):N_t)=0;

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;


Signal=fft(Spectrum_New);
Signal=Signal;
Spectrum_New=Spectrum_New(1:N_f);

%% After Water


Thickness_Water_Layer=(50:50:550).*1E-6;  %mm
Thickness_Water_Layer_micron=Thickness_Water_Layer.*1E6;
Data_Water=importdata('120503_Water Dispersion.txt');

Wavelength_Water=Data_Water(:,1);           %nm
%n_Water_Wavelength=Data_Water(:,2)+i*Data_Water(:,3); 

Frequency_Water=C./(Wavelength_Water*1E-6);


n_Water=interp1(Frequency_Water,Data_Water(:,2)+i*Data_Water(:,3),Frequency_New);%+i*interp1(Frequency_Water,Data_Water(:,3),Frequency_New);
n_Water(isnan(n_Water))=n_Water(find(n_Water>1,1,'first'));
dn=diff(n_Water);
dn(length(dn)+1)=dn(length(dn));
n_Water_Group=n_Water./(1+dn./n_Water);

n_Water_Group_Center=n_Water_Group(find(Frequency_Water<Frequency_Center,1,'first'));

[value_1 index_1]=max(Signal);
Position_Upper_peak=Position_micron(index_1);
OPD_Measured_micron(1:length(Thickness_Water_Layer_micron))=0;
for j=1:length(Thickness_Water_Layer)

Spectrum_After_Water=Spectrum_New.*exp(i*4*pi.*Frequency_New/C.*n_Water.*(Thickness_Water_Layer(j)));

Spectrum_After_Water(Frequency_New<Min_Frequency)=0;
Spectrum_After_Water(isnan(Spectrum_After_Water))=0;

%plot(Wavelength_micron,real(n_Water),Wavelength_micron,real(Spectrum_After_Water));

Spectrum_After_Water((N_f+1):N_t)=0;
Signal_After_Water=fft(Spectrum_After_Water);
Spectrum_After_Water=Spectrum_After_Water(1:N_f);

[value_2 index_2]=max(abs(Signal_After_Water));

OPD_Measured_micron(j)=Position_micron(index_2)-Position_Upper_peak;
Signal_After_Water_Record(:,j)=Signal_After_Water./abs(max(Signal));

end

plot(Wavelength_micron,Spectrum_After_Water);
plot(Thickness_Water_Layer_micron,OPD_Measured_micron);


%fitting=polyfit(Thickness_Water_Layer_micron,OPD_Measured_micron,1);

%OPD_Fitted=fitting(1).*Thickness_Water_Layer_micron+fitting(2);

%plot(Thickness_Water_Layer_micron,OPD_Measured_micron,Thickness_Water_Layer_micron,OPD_Fitted);


plot((Position_micron-Position_Upper_peak)./fitting(1),Signal_After_Water_Record);

xlabel('Cornea Thickness (micron)');
ylabel('PSF (norm.)');


plot(Thickness_Water_Layer_micron,(OPD_Measured_micron-OPD_Fitted)/fitting(1));

xlabel('Cornea Thickness (micron)');
ylabel('Measurement Error Caused by Dispersion (micron)');

%plot(Position_micron,real(Signal),Position_micron,real(Signal_After_Water));

%% FWHM calculation

for p=1:length(Thickness_Water_Layer)
    FWHM(p)=Position_micron(find(Signal_After_Water_Record(:,p)/max(Signal_After_Water_Record(:,p))>0.5,1,'last'))-Position_micron(find(Signal_After_Water_Record(:,p)/max(Signal_After_Water_Record(:,p))>0.5,1,'first'));
end
%% ROI
%Time=Time(1:round(length(Time)*ROI_ratio));
%Position=Position(1:round(length(Position)*ROI_ratio));
%Position_micron=Position_micron(1:round(length(Position_micron)*ROI_ratio));

%% Shift

%[maxvalue maxindex]=max(Signal_Envelope,[],1);
%needshift=-maxindex;

%Signal=circshift(Signal,needshift);
%Signal_Carrier=circshift(Signal_Carrier,needshift);
%Signal_Envelope=circshift(Signal_Envelope,needshift);


DC_cutoff=20;   %micron

%Signal_Bscan_Envelope_Cut=Signal_Bscan_Envelope(Position_micron>DC_cutoff)/max(Signal_Bscan_Envelope(Position_micron>DC_cutoff));
%Position_micron_cut=Position_micron(Position_micron>DC_cutoff);
%FWHM=Position_micron_cut(find(Signal_Bscan_Envelope_Cut>0.5,1,'last'))-Position_micron_cut(find(Signal_Bscan_Envelope_Cut>0.5,1,'first'));

%plot(Wavelength_micron,Spectrum_New);
%plot(Position_micron,real(Signal),Position_micron,abs(Signal));

%dlmwrite('Profile.txt',Profile,'delimiter','\t','newline','pc');

%dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');
%plot(Position_micron,Signal_Bscan_Envelope);