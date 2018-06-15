clear all;

P=0.001;            %(log(Watt), incident power)
%Sample_Reflectance=((1.366-1)/(1.366+1))^2;
%Sample_Reflectance=((1.376-1.366)/(1.376+1.366))^2;
Sample_Reflectance=1;
PsdBm=[-30:1:0]'; %incident to sample
wavelength=1.03;    %(micron, 暫時用center wavelength, 實際上應該是要用spectrum)
%Rslog=[-10:0.1:0];    %(log10(reflectance )of sample arm)
R=((1.5-1)/(1.5+1))^2;
%Rs=10.^Rslog;
Ps=(10.^(PsdBm/10))*1E-3;
Ps=0.3E-3;
%Ps=0.5*Rs*P;
%Pr=0.5*P;           %(power back from ref arm, assumed rr=1)
BW=37.5E3;              %Hz
ti=1/BW/2;
Ps=Ps/2;             %Beamsplitter;

Pr=Ps;

Ps=Sample_Reflectance*Ps;

qe=0.5;                             %(CCD quantum efficiency (electron to photon))
N=4096;                             %number of pixel


h=6.626E-34;        %plank constant
c=3E8;              %light speed
d_freq=(c/(1.20E-6))-(c/(1.40E-6));


SigTD=(2*qe*((Pr.*Ps).^0.5)/h/(c/(wavelength*10^-6)))*R;        %for one pixel after FFT, N is timed back asuumed FFT
shut_TD2=(qe*(Pr+Ps)/(h*(c/(wavelength*10^-6))))*BW;                   %noise is not timed back after FFT
dark_TD2=((qe*(1.7E-12)/h/(c/(wavelength*10^-6))*ti)^2)*BW;        %for one pixel after FFT, N is timed back asuumed FFT
excess_TD2=((qe/(h*(c/(wavelength*10^-6))))^2)*((Pr+Ps).^2)/d_freq*BW;

%SNRshut=10*log10((Sig.^2)./shut_2);
%SNRexcess=10*log10((Sig.^2)./excess_2);
%SNRdark=10*log10((Sig.^2)./dark_2);
SNRshut_TD=10*log10((SigTD.^2)./(shut_TD2))';
SNRdark_TD=10*log10((SigTD.^2)./(dark_TD2))';
SNRexcess_TD=10*log10((SigTD.^2)./(excess_TD2))';
SNRtotal_TD=10*log10((SigTD.^2)./(shut_TD2+excess_TD2+dark_TD2))';

plot(PsdBm,SNRtotal_TD,PsdBm,SNRdark_TD,PsdBm,SNRshut_TD,PsdBm,SNRexcess_TD);

xlabel('Incident Power on Sample (dBm)');
ylabel('SNR (dB)');
