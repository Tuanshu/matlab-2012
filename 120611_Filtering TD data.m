clear all

%% Setting
SR=40000;
Stage_Speed=2;  %mm/s

Stage_Speed_MTS=(Stage_Speed*1E-3);

C=3E8;

cd('C:\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata('Data1.txt');

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


Window=(gaussmf(Frequency,[0.4E14 2.5E14]));
Window(Frequency>2.5E14)=1;
Window=Window';
Spectrum=Window.*Spectrum;
Spectrum((round(length(Spectrum)/2)+1):end)=0;

plot(Frequency,Spectrum);

Data_New=ifft(Spectrum);
Data_New=Data_New(1:N_t);

Max_Wavelength_micron=1;

Min_Wavelength_micron=0.6;

Max_Wavelength_micron_index=find(Wavelength_micron<Max_Wavelength_micron,1,'first');

Min_Wavelength_micron_index=find(Wavelength_micron<Min_Wavelength_micron,1,'first');

plot(Wavelength_micron(Max_Wavelength_micron_index:Min_Wavelength_micron_index),Spectrum(Max_Wavelength_micron_index:Min_Wavelength_micron_index));
plot(Position_micron,Data_New);



%cd('D:\120222\');

%dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');
%plot(Position_micron,Signal_Bscan_Envelope);