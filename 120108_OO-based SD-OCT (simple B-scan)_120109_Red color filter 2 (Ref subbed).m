clear all

%% Setting
%%OO: 10ms ave 25 -> 250ms -> 4/sec
%%APT: 0.002mm/sec

thickness_temp=2.68E-6;
n_should=1.5;

Starting_pixel_f_considered=250;
Ending_pixel_f_considered=350;

range_specified=1;

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=8192*16;
ROI_ratio=1/16;                  %only consider the first ROI_ratio data in TD

DC_cutoff=600;                 %in TD

array=7:406;

Pixel_spacing=0.5;  %micron

Lateral_Position=Pixel_spacing*(array-array(1));


%% Data Loading

Signal_Bscan_Carrier(1:(N_t*ROI_ratio),1:length(array))=0;

Signal_Bscan_Envelope(1:(N_t*ROI_ratio),1:length(array))=0;

for jj=1:length(array)

cd('D:\120109_Red Color Filter 2\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata(sprintf('D%i.txt',array(jj)));
%Is_o=importdata(sprintf('%i',array(jj)));

cd('D:\120109_Red Color Filter 2_Ref\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Ref=importdata(sprintf('D%i.txt',array(jj)));

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

%% Shift

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
imagesc(10*log10(Signal_Bscan_Envelope));

profile_Sub=profile-profile_Ref;
plot(Lateral_Position,profile,Lateral_Position,profile_Ref);

plot(Lateral_Position,profile_Sub);

dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');