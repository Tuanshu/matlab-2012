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

Wavelength_micron=(C./Frequency_New)*1E6;

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



Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;


Spectrum_for_FFT=Spectrum_New;
Spectrum_for_FFT((N_f+1):N_t,:)=0;
Signal_Old=fft(Spectrum_for_FFT);
Signal_Old((size(Signal_Old,1)/2+1):end,:)=0;
Window=(gaussmf(Position_micron,[2 20]));
Window(Position_micron>20)=1;
Window=[Window Window];
Signal_Windowed=Signal_Old.*Window;
plot(Position_micron,real(Signal_Windowed));
Spectrum_for_FFT=ifft(Signal_Windowed);
Wavelength_micron_for_FFT=Wavelength_micron;
Wavelength_micron_for_FFT((N_f+1):N_t,:)=0;
Wavelength_micron_for_FFT(Wavelength_micron_for_FFT==0)=1;


%% Shift
shift=-6:0.1:-4;       %micron
amp=0.3:0.1:0.9;
MIN=9999999999999999999999;
Range_Lower=17;     %micron
Range_Upper=20;     %micron
Range_Lower_index=find(Position_micron>Range_Lower,1,'first');
Range_Upper_index=find(Position_micron>Range_Upper,1,'first');
Weight_array=(1:size(Signal_Old,1))';
Weight_array(Weight_array<Range_Lower_index)=0;
Weight_array(Weight_array>Range_Upper_index)=0;
Weight_array=(Weight_array-0.5*Range_Lower_index)/Range_Upper_index;

for p=1:length(shift)
    for q=1:length(amp)
        Signal_Diff=fft(Spectrum_for_FFT(:,1))-fft(amp(q)*Spectrum_for_FFT(:,2).*exp(i*4*pi./Wavelength_micron_for_FFT.*shift(p)));
        if sum(real(Signal_Diff(Range_Lower_index:Range_Upper_index)).^2)<MIN
            MIN=sum(real(Signal_Diff(Range_Lower_index:Range_Upper_index)).^2);
            shift_best=shift(p);
            amp_best=amp(q);
        end
    end
end

plot(Position_micron,fft(Spectrum_for_FFT(:,1))-fft(amp_best*Spectrum_for_FFT(:,2).*exp(i*4*pi./Wavelength_micron_for_FFT.*shift_best)));

plot(Position_micron,fft(Spectrum_for_FFT(:,1)),Position_micron,fft(amp_best*Spectrum_for_FFT(:,2).*exp(i*4*pi./Wavelength_micron_for_FFT.*shift_best)));