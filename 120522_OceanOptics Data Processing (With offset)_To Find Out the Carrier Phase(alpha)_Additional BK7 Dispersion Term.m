clear all

%% Setting

Number_of_Iteration=200;
thickness_temp=2.68E-6;
n_should=1.5;

Starting_pixel_f_considered=250;
Ending_pixel_f_considered=350;

range_specified=1;

Max_Wavelength=700;             %nm
Min_Wavelength=450;             %nm
N_f=8192;
N_t=8192*64;
ROI_ratio=1/32;                  %only consider the first ROI_ratio data in TD

DC_cutoff=0;                 %in TD

array=90;

%% Data Loading

Signal_Bscan_Carrier(1:(N_t*ROI_ratio),1:length(array))=0;

Signal_Bscan_Envelope(1:(N_t*ROI_ratio),1:length(array))=0;

for jj=1:length(array)

%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
if array(jj) <= 99
    cd('D:\120522\1S\');
    Data=importdata(sprintf('D%i.txt',array(jj)));
elseif array(jj) > 99
    cd('D:\120522\2S\');
    Data=importdata(sprintf('D%i.txt',array(jj)-100));
end
%cd('D:\120222\');
%Data_R=importdata('R1 Ref.txt');
%cd('D:\120222\R1 Sam\');
%Data_S=importdata(sprintf('D%i.txt',array(jj)));

if jj==1
    
Wavelength=Data(:,1);           %nm

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

%% Temporal Filtering

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;

Spectrum_New((N_f+1):N_t)=0;
Signal=fft(Spectrum_New);
Window=(gaussmf(Position_micron,[2 10]));
Window(Position_micron>10)=1;
Signal=Signal.*Window;
Signal((round(length(Signal)/2)+1):end)=0;
Spectrum_New=ifft(Signal);
Spectrum_New=Spectrum_New(1:N_f);

%% BK7


Wavelength_micron=(C./Frequency_New)*1E6;

% Assumed n1 = BK7
C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;


n_BK7=(C1*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C2)+C3*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C4)+C5*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C6)+1).^0.5;

%%
%Carrier_Function=Amp*cos(4*pi*Frequency*Thickness/C+Phi);

A= @(Amp,Thickness,Phi,Thickness_BK7) abs(Spectrum_New).*(Amp*exp(i*4*pi*Frequency_New*Thickness/C+Phi).*exp(i*4*pi*Frequency_New.*n_BK7*Thickness_BK7/C));
dAmp=0.0001;
dThickness=1E-12;
dPhi=0.0001;
dThickness_BK7=1E-12;
dA_dAmp= @(Amp,Thickness,Phi,Thickness_BK7) (A(Amp+dAmp,Thickness,Phi,Thickness_BK7)-A(Amp,Thickness,Phi,Thickness_BK7))/dAmp;
dA_dThickness= @(Amp,Thickness,Phi,Thickness_BK7) (A(Amp,Thickness+dThickness,Phi,Thickness_BK7)-A(Amp,Thickness,Phi,Thickness_BK7))/dThickness;
dA_dPhi= @(Amp,Thickness,Phi,Thickness_BK7) (A(Amp,Thickness,Phi+dPhi,Thickness_BK7)-A(Amp,Thickness,Phi,Thickness_BK7))/dPhi;
dA_dThickness_BK7= @(Amp,Thickness,Phi,Thickness_BK7) (A(Amp,Thickness,Phi,Thickness_BK7+dThickness_BK7)-A(Amp,Thickness,Phi,Thickness_BK7))/dThickness;

%Use the envelope of spectrum as the weight function of the averaging of
%delta_X
Current_Iteration=0;
ita_Amp=0.01;
ita_Thickness=0.01;
ita_Phi=0.01;
ita_Thickness_BK7=0.01;

Amp=1;
Thickness=18E-6;
Phi=0;
Thickness_BK7=0;
while (Current_Iteration<Number_of_Iteration)
    Current_Iteration=Current_Iteration+1;
    Spectral_Array_delta_Amp=((Spectrum_New)-A(Amp,Thickness,Phi,Thickness_BK7))./dA_dAmp(Amp,Thickness,Phi,Thickness_BK7);
    Spectral_Array_delta_Thickness=((Spectrum_New)-A(Amp,Thickness,Phi,Thickness_BK7))./dA_dThickness(Amp,Thickness,Phi,Thickness_BK7);
    Spectral_Array_delta_Phi=((Spectrum_New)-A(Amp,Thickness,Phi,Thickness_BK7))./dA_dPhi(Amp,Thickness,Phi,Thickness_BK7);
    Spectral_Array_delta_Thickness_BK7=((Spectrum_New)-A(Amp,Thickness,Phi,Thickness_BK7))./dA_dThickness_BK7(Amp,Thickness,Phi,Thickness_BK7);
    Spectral_Array_delta_Amp(isnan(Spectral_Array_delta_Amp))=0;
    Spectral_Array_delta_Thickness(isnan(Spectral_Array_delta_Thickness))=0;
    Spectral_Array_delta_Phi(isnan(Spectral_Array_delta_Phi))=0;
    Spectral_Array_delta_Thickness_BK7(isnan(Spectral_Array_delta_Thickness_BK7))=0;
    Spectral_Array_delta_Amp(isinf(Spectral_Array_delta_Amp))=0;
    Spectral_Array_delta_Thickness(isinf(Spectral_Array_delta_Thickness))=0;
    Spectral_Array_delta_Phi(isinf(Spectral_Array_delta_Phi))=0;
    Spectral_Array_delta_Thickness_BK7(isinf(Spectral_Array_delta_Thickness_BK7))=0;
    Amp=real(ita_Amp*sum(Spectral_Array_delta_Amp.*abs(Spectrum_New))./sum(abs(Spectrum_New)))+Amp;
    Thickness=real(ita_Thickness*sum(Spectral_Array_delta_Thickness.*abs(Spectrum_New))./sum(abs(Spectrum_New)))+Thickness;
    Phi=real(ita_Phi*sum(Spectral_Array_delta_Phi.*abs(Spectrum_New))./sum(abs(Spectrum_New)))+Phi;
    Thickness_BK7=real(ita_Thickness_BK7*sum(Spectral_Array_delta_Thickness_BK7.*abs(Spectrum_New))./sum(abs(Spectrum_New)))+Thickness_BK7;
end

plot(Frequency_New,A(Amp,Thickness,Phi,Thickness_BK7),Frequency_New,Spectrum_New);

end

%% To time domain