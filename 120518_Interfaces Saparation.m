clear all

%% Setting

Max_Wavelength=700;             %nm
Min_Wavelength=450;             %nm
N_f=8192;
N_t=8192*16;
ROI_ratio=1/4;                  %only consider the first ROI_ratio data in TD

DC_cutoff=1000;                 %in TD

array=14:15;

%% Data Loading

Signal_Bscan_Carrier(1:(N_t*ROI_ratio),1:length(array))=0;

Signal_Bscan_Envelope(1:(N_t*ROI_ratio),1:length(array))=0;

cd('D:\Ocean Optics\');
Data_1=importdata('Data_1.txt');
Data_2=importdata('Data_2.txt');

Wavelength=Data_1(:,1);           %nm

C=3E8;

Frequency=C./(Wavelength*1E-9);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_New=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency_New=Frequency_New';

Spectrum_1=Data_1(:,2)-Data_1(1,2);%-Data_R(:,2)-Data_S(:,2);
Spectrum_Frequency_1=(Spectrum_1.*((Wavelength*1E-9).^2)/C)/max(Spectrum_1.*((Wavelength*1E-9).^2)/C);
Spectrum_New_1=interp1(Frequency,Spectrum_Frequency_1,Frequency_New);

Spectrum_New_1(isnan(Spectrum_New_1))=0;
Spectrum_New_1(Frequency_New<Min_Frequency)=0;

Spectrum_2=Data_2(:,2)-Data_2(1,2);%-Data_R(:,2)-Data_S(:,2);
Spectrum_Frequency_2=(Spectrum_2.*((Wavelength*1E-9).^2)/C)/max(Spectrum_2.*((Wavelength*1E-9).^2)/C);
Spectrum_New_2=interp1(Frequency,Spectrum_Frequency_2,Frequency_New);

Spectrum_New_2(isnan(Spectrum_New_2))=0;
Spectrum_New_2(Frequency_New<Min_Frequency)=0;
%plot(Frequency_New,Spectrum_New);
Spectrum_New=[Spectrum_New_1 Spectrum_New_2];
%% To time domain

Spectrum_New((N_f+1):N_t,:)=0;



Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;


Signal=fft(Spectrum_New,[],1);
Spectrum_New=Spectrum_New(1:N_f,:);
Signal((size(Signal,1)/2+1):end,:)=0;


plot(Position_micron,Signal);
%% Shift
shift=-20:20;
amp=0.9:0.01:1.1;
MIN=99900000;
Range_Lower=25;     %micron
Range_Upper=35;     %micron
Range_Lower_index=find(Position_micron>Range_Lower,1,'first');
Range_Upper_index=find(Position_micron>Range_Upper,1,'first');
Weight_array=(1:size(Signal,1))';
Weight_array(Weight_array<Range_Lower_index)=0;
Weight_array(Weight_array>Range_Upper_index)=0;
Weight_array=(Weight_array-0.5*Range_Lower_index)/Range_Upper_index;

for p=1:length(shift)
    for q=1:length(amp)
        Signal_Diff=(Signal(:,1)-amp(q)*circshift(Signal(:,2),shift(p))).*Weight_array;
        if mean(abs(Signal_Diff(Range_Lower_index:Range_Upper_index)).^2)<MIN
            MIN=mean(abs(Signal_Diff(Range_Lower_index:Range_Upper_index)).^2);
            shift_best=shift(p);
            amp_best=amp(q);
        end
    end
end

plot(Position_micron,Signal(:,1),Position_micron,amp_best*circshift(Signal(:,2),shift_best));

plot(Position_micron,Signal(:,1)-amp_best*circshift(Signal(:,2),shift_best));
%% Amp Varying