clear all

%% Setting


thickness_temp=2.68E-6;
n_should=1.5;

Starting_pixel_f_considered=250;
Ending_pixel_f_considered=350;

range_specified=1;

Max_Wavelength=700;             %nm
Min_Wavelength=400;             %nm
N_f=8192;
N_t=8192*8;
ROI_ratio=1/4;                  %only consider the first ROI_ratio data in TD

DC_cutoff=20;                 %in TD

array=0:181;

current_sep=1900;

%% Data Loading

Signal_Bscan_Carrier(1:(N_t*ROI_ratio),1:length(array))=0;

Signal_Bscan_Envelope(1:(N_t*ROI_ratio),1:length(array))=0;
for jj=1:length(array)

cd('C:\TuanShu\120905_YPR Bscan_2\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata(sprintf('D%i.txt',array(jj)));

Wavelength=Data(:,1);           %nm

%Data=importdata('inter.txt')-importdata('sam.txt')-importdata('ref.txt')+importdata('bs.txt');
%cd('D:\120222\');
%Data_R=importdata('R1 Ref.txt');
%cd('D:\120222\R1 Sam\');
%Data_S=importdata(sprintf('D%i.txt',array(jj)));

if jj==1

C=3E8;

Frequency=C./(Wavelength*1E-9);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_New=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency_New=Frequency_New';

end

Spectrum=Data(:,2)-Data(1,2);%-Data_R(:,2)-Data_S(:,2);
Spectrum_Frequency=(Spectrum.*((Wavelength*1E-9).^2)/C)/max(Spectrum.*((Wavelength*1E-9).^2)/C);
Spectrum_New=interp1(Frequency,Spectrum_Frequency,Frequency_New);

Spectrum_New(isnan(Spectrum_New))=0;
Spectrum_New(Frequency_New<Min_Frequency)=0;
%plot(Frequency_New,Spectrum_New);

%% To time domain

Spectrum_New((N_f+1):N_t)=0;


if jj==1

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;

Time=Time(1:round(length(Time)*ROI_ratio));
Position=Position(1:round(length(Position)*ROI_ratio));
Position_micron=Position_micron(1:round(length(Position_micron)*ROI_ratio));
end

Signal=fft(Spectrum_New);

Spectrum_New=Spectrum_New(1:N_f);
Signal=Signal(1:round(length(Signal)*ROI_ratio));

Window=(gaussmf(Position_micron,[2 8]));
Window(Position_micron>8)=1;
%Signal=Signal.*Window;
%Signal(1:DC_cutoff)=0;
Signal_Carrier=real(Signal);
Signal_Envelope=abs(Signal);

%% ROI


%% Shift



Signal_Bscan_Carrier(:,jj)=Signal_Carrier;
Signal_Bscan_Envelope(:,jj)=Signal_Envelope;

%DC_cutoff=20;   %micron

Signal_Bscan_Envelope_Cut(:,jj)=Signal_Bscan_Envelope(Position_micron>DC_cutoff,jj)/max(Signal_Bscan_Envelope(Position_micron>DC_cutoff,jj));
Signal_Bscan_Carrier_Cut(:,jj)=Signal_Bscan_Carrier(Position_micron>DC_cutoff,jj)/max(Signal_Bscan_Envelope(Position_micron>DC_cutoff,jj));

[maxvalue maxindex]=max(Signal_Bscan_Envelope_Cut(1:current_sep,jj));
needshift=-maxindex+1000;
if abs(maxindex+200 -current_sep)<100
current_sep=maxindex+200;
end
profile(jj)=maxindex;
%Signal=circshift(Signal,needshift);
Signal_Bscan_Carrier_Cut(:,jj)=circshift(Signal_Bscan_Carrier_Cut(:,jj),needshift);
Signal_Bscan_Envelope_Cut(:,jj)=circshift(Signal_Bscan_Envelope_Cut(:,jj),needshift);


Position_micron_cut=Position_micron(Position_micron>DC_cutoff);
%FWHM=Position_micron_cut(find(Signal_Bscan_Envelope_Cut>0.5,1,'last'))-Position_micron_cut(find(Signal_Bscan_Envelope_Cut>0.5,1,'first'));
%FWHM_should=2*Position_micron(find((Signal_Bscan_Envelope./max(Signal_Bscan_Envelope))<0.5,1,'first'));
end

plot(Position_micron_cut,Signal_Bscan_Carrier_Cut);
xlabel('OPD (micron)');
ylabel('Interference Signal');

plot(Position_micron_cut,Signal_Bscan_Envelope_Cut);



plot(Wavelength,Spectrum);
xlabel('Wavelength (nm)');
ylabel('Spectrum');

imagesc(10*log10(Signal_Bscan_Envelope_Cut(800:1500,:)/max(max(Signal_Bscan_Envelope_Cut))),'xdata',(array*10)/1000,'ydata',Position_micron(800:1500)-Position_micron(1000));
colormap(gray);
xlabel('Lateral Position (mm)');
ylabel('Optical Path Difference (micron)');
dlmwrite('xdata.txt',(array*10)/1000,'delimiter','\t','newline','pc');

dlmwrite('ydata.txt',Position_micron(800:1500)-Position_micron(1000),'delimiter','\t','newline','pc');

dlmwrite('data.txt',10*log10(Signal_Bscan_Envelope_Cut(800:1500,:)/max(max(Signal_Bscan_Envelope_Cut))),'delimiter','\t','newline','pc');
%dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');
%plot(Position_micron,Signal_Bscan_Envelope);