clear all

%% Setting
SR=40000;
Stage_Speed=2;  %mm/s

Stage_Speed_MTS=(Stage_Speed*1E-3);

C=3E8;

TH_first_peak=0.02;
TH_second_peak=0.005;
Delay_first_peak=80000; %the smaller one, pixel
Min_thickness=12000; %pixel

Max_thickness=16000; %pixel

cd('C:\');

qqq=1;
xxx=1;
yyy=1;

%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata('Data5.txt');

%%
N_t=length(Data);
N_f=N_t;

Time_Stage=1/SR:1/SR:(1/SR)*N_t;

Time=Time_Stage*Stage_Speed_MTS/C;
Position_micron=Time*C*1E6;

dTime=(Time(2)-Time(1))*2;      %*2 becuase of round trip

Frequency_Max=1/dTime;

Frequency=Frequency_Max/N_f:Frequency_Max/N_f:Frequency_Max;

Wavelength_micron=(C./Frequency)*1E6;

Spectrum=fft(Data,N_f);


Window1=(gaussmf(Frequency,[0.2E14 1.5E14]));
Window1(Frequency>1.5E14)=1;
Window2=(gaussmf(Frequency,[0.2E14 5.5E14]));
Window2(Frequency<5.5E14)=1;
Window=(Window1.*Window2)';
Spectrum=Window.*Spectrum;
Spectrum((round(length(Spectrum)/2)+1):end)=0;

plot(Frequency,Spectrum);

Data_New=ifft(Spectrum);
Data_New=Data_New(1:N_t);

Max_Wavelength_micron=1;

Min_Wavelength_micron=0.6;

Max_Wavelength_micron_index=find(Wavelength_micron<Max_Wavelength_micron,1,'first');

Min_Wavelength_micron_index=find(Wavelength_micron<Min_Wavelength_micron,1,'first');

Data_New(1:40)=Data_New(40);
Data_New((length(Data_New)-40):end)=Data_New((length(Data_New)-40));


plot(Position_micron,Data_New,Position_micron,abs(Data_New));
xlabel('OPD (micron)');
SD=std(real(Data_New(1E4:3E4)));

SNR=10*log10((max(abs(Data_New))^2)/(SD^2));
%SNR_even=10*log10((max(abs(max_array_even))^2)/(SD^2));
%SNR_odd=10*log10((max(abs(max_array_odd))^2)/(SD^2));
%plot(Wavelength_micron,Spectrum);
%xlabel('Wavelength (micron)');
%label('Interference Spectrum');
%cd('D:\120222\');

%dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');
%plot(Position_micron,Signal_Bscan_Envelope);