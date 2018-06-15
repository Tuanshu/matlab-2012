clear all

%% Setting
%%OO: 50ms ave 20 -> 1000ms -> 1/sec
%%APT: 0.001mm/sec

%%Position error: >0 mean the reference is away to the DC term

thickness_temp=2.68E-6;
n_should=1.5;

Position_Data=50;
Position_Ref=950;
Position_Ref_0=800;  %used for tilting

Starting_pixel_f_considered=250;
Ending_pixel_f_considered=350;

range_specified=1;

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=8192*16;
ROI_ratio=1/16;                  %only consider the first ROI_ratio data in TD

DC_cutoff=600;                 %in TD

array=1:1000;

Pixel_spacing=1;  %micron

Lateral_Position=Pixel_spacing*(array-array(1));

%% Data Loading

Signal_Bscan_Carrier(1:(N_t*ROI_ratio),1:length(array))=0;

Signal_Bscan_Envelope(1:(N_t*ROI_ratio),1:length(array))=0;

for jj=1:length(array)

cd('D:\120110_1mm green color filter_again\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata(sprintf('D%i.txt',array(jj)));
%Is_o=importdata(sprintf('%i',array(jj)));

cd('D:\120110_1mm green color filter_ref\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Ref=importdata(sprintf('D%i.txt',array(jj)));
%Is_o=importdata(sprintf('%i',array(jj)));


if jj==1
    
Wavelength=Data(:,1);           %nm

C=3E8;

Frequency=C./(Wavelength*1E-9);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_New=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency_New=Frequency_New';

end

Spectrum=Data(:,2);
Spectrum_Frequency=(Spectrum.*((Wavelength*1E-9).^2)/C)/max(Spectrum.*((Wavelength*1E-9).^2)/C);
Spectrum_New=interp1(Frequency,Spectrum_Frequency,Frequency_New);

Spectrum_New(isnan(Spectrum_New))=0;
Spectrum_New(Frequency_New<Min_Frequency)=0;
%plot(Frequency_New,Spectrum_New);


Spectrum_Ref=Ref(:,2);
Spectrum_Ref_Frequency=(Spectrum_Ref.*((Wavelength*1E-9).^2)/C)/max(Spectrum_Ref.*((Wavelength*1E-9).^2)/C);
Spectrum_Ref_New=interp1(Frequency,Spectrum_Ref_Frequency,Frequency_New);

Spectrum_Ref_New(isnan(Spectrum_Ref_New))=0;
Spectrum_Ref_New(Frequency_New<Min_Frequency)=0;

%% To time domain

Spectrum_New((N_f+1):N_t)=0;

Spectrum_Ref_New((N_f+1):N_t)=0;

if jj==1

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;

end

Signal=fft(Spectrum_New);
Spectrum_New=Spectrum_New(1:N_f);
Signal=Signal(1:round(length(Signal)*ROI_ratio));
Signal(1:DC_cutoff)=0;
Signal_Carrier=real(Signal);
Signal_Envelope=abs(Signal);


Signal_Ref=fft(Spectrum_Ref_New);
Spectrum_Ref_New=Spectrum_Ref_New(1:N_f);
Signal_Ref=Signal_Ref(1:round(length(Signal_Ref)*ROI_ratio));
Signal_Ref(1:DC_cutoff)=0;
Signal_Ref_Carrier=real(Signal_Ref);
Signal_Ref_Envelope=abs(Signal_Ref);

%% ROI

if jj==1

Time=Time(1:round(length(Time)*ROI_ratio));
Position=Position(1:round(length(Position)*ROI_ratio));
Position_micron=Position_micron(1:round(length(Position_micron)*ROI_ratio));

end

%% Profiling


[value_max index_max]=max(Signal_Envelope);
profile(jj)=Position_micron(index_max);

[value_Ref_max index_Ref_max]=max(Signal_Ref_Envelope);
profile_Ref(jj)=Position_micron(index_Ref_max);
%[maxvalue maxindex]=max(Signal_Envelope,[],1);
%needshift=-maxindex;

%Signal=circshift(Signal,needshift);
%Signal_Carrier=circshift(Signal_Carrier,needshift);
%Signal_Envelope=circshift(Signal_Envelope,needshift);

Signal_Bscan_Carrier(:,jj)=Signal_Carrier;
Signal_Bscan_Envelope(:,jj)=Signal_Envelope;
end
imagesc(10*log10(Signal_Bscan_Envelope),'xdata',Lateral_Position,'ydata',Position_micron);
profile_Sub=profile-profile_Ref;

Position_Error_Stage=profile_Ref(Position_Ref)-profile_Ref(Position_Data);

Height_Data=profile_Sub(Position_Data);
Height_Ref=profile_Sub(Position_Ref);
Height_Ref_0=profile_Sub(Position_Ref_0); %Height_Ref and Height_Ref2 should be the same

slope=(Height_Ref-Height_Ref_0)/(Position_Ref-Position_Ref_0);
Position_Error_Tilt=(Position_Ref-Position_Data)*slope;

for j=1:length(profile_Sub)
    profile_tilt(j)=profile_Sub(j)-slope*(j-1);
end

Real_Difference_1=profile_tilt(Position_Data)-profile_tilt(Position_Ref);
Real_Difference_2=profile(Position_Data)-(profile(Position_Ref)-Position_Error_Stage-Position_Error_Tilt);

plot(Lateral_Position,profile,Lateral_Position,profile_Ref);

plot(Lateral_Position,profile_Sub);

Position_Error=Position_Error_Stage+Position_Error_Tilt;

plot(Lateral_Position,profile_tilt);


Thickness=abs(mean(profile_Sub(1:20))-mean(profile_Sub((length(profile_Sub)-19):length(profile_Sub))));

dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');