clear all

%% Setting
%%OO: 20ms ave 5 -> 100ms -> 10/sec
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

array=0:999;

Pixel_spacing=1;  %micron

Lateral_Position=Pixel_spacing*(array-array(1));

%% Data Loading

Signal_Bscan_Carrier(1:(N_t*ROI_ratio),1:length(array))=0;

Signal_Bscan_Envelope(1:(N_t*ROI_ratio),1:length(array))=0;

for jj=1:length(array)

cd('D:\120213\5-2\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata(sprintf('D%i.txt',array(jj)));
%Is_o=importdata(sprintf('%i',array(jj)));

cd('D:\120213\5-R2\');
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
index(jj)=index_max;

[value_Ref_max index_Ref_max]=max(Signal_Ref_Envelope);
profile_Ref(jj)=Position_micron(index_Ref_max);

index_Ref(jj)=index_Ref_max;
%[maxvalue maxindex]=max(Signal_Envelope,[],1);
%needshift=-maxindex;

%Signal=circshift(Signal,needshift);
%Signal_Carrier=circshift(Signal_Carrier,needshift);
%Signal_Envelope=circshift(Signal_Envelope,needshift);

Signal_Bscan_Carrier(:,jj)=Signal_Carrier;
Signal_Bscan_Envelope(:,jj)=Signal_Envelope;
end

[max_value max_index]=max(Signal_Bscan_Envelope(:,2000));
XX=Position_micron(max_index);

imagesc(10*log10(Signal_Bscan_Envelope/max(max(Signal_Bscan_Envelope))),'xdata',Lateral_Position,'ydata',Position_micron-XX);

clear Signal_Bscan_Carrier
Signal_Bscan_Envelope_Shift=Signal_Bscan_Envelope;

for jj=1:length(array)
    Signal_Bscan_Envelope_Shift(:,jj)=circshift(Signal_Bscan_Envelope(:,jj),-index_Ref(jj)+3000);
end
    
[max_value max_index]=max(Signal_Bscan_Envelope_Shift(:,2000));
XX=Position_micron(max_index);

imagesc(10*log10(Signal_Bscan_Envelope_Shift/max(max(Signal_Bscan_Envelope_Shift))),'xdata',Lateral_Position,'ydata',Position_micron-XX);

profile_Sub=profile-profile_Ref;
index_Sub=index-index_Ref;


fitting=polyfit([1400:1800],index_Sub(1400:1800),1);


fitting_position=polyfit([1400:1800],profile_Sub(1400:1800),1);

fitted_curve=fitting(1).*[1:length(index_Sub)]+fitting(2);

fitted_curve_position=fitting_position(1).*[1:length(index_Sub)]+fitting_position(2);

profile_Tilt=profile_Sub-fitted_curve_position;

plot([1:length(index_Sub)],index_Sub,[1:length(index_Sub)],fitted_curve);

plot([1:length(index_Sub)],profile_Sub,[1:length(index_Sub)],fitted_curve_position);

Signal_Bscan_Envelope_Tilt=Signal_Bscan_Envelope;

for jj=1:length(array)
    Signal_Bscan_Envelope_Tilt(:,jj)=circshift(Signal_Bscan_Envelope_Shift(:,jj),round(-fitted_curve(jj)+3000));
end


[max_value max_index]=max(Signal_Bscan_Envelope_Tilt(:,2000));
XX=Position_micron(max_index);

imagesc(10*log10(Signal_Bscan_Envelope_Tilt/max(max(Signal_Bscan_Envelope_Tilt))),'xdata',Lateral_Position,'ydata',Position_micron-XX);

plot(Lateral_Position,profile,Lateral_Position,profile_Ref);

plot(Lateral_Position,profile_Sub);

Position_Error=Position_Error_Stage+Position_Error_Tilt;

plot(Lateral_Position,profile_tilt);


Thickness=abs(mean(profile_Sub(1:20))-mean(profile_Sub((length(profile_Sub)-19):length(profile_Sub))));

dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');