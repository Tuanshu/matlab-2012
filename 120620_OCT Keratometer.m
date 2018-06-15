clear all
%%

cd('D:\Users\TuanShu');
Data=importdata('121115_Yb Fiber Output (300mA pump).txt');      % Data_2: the glass data 1

C=3E8;

Wavelength_Max=1.1; %micron
Wavelength_Min=0.9;

Lateral_Resolution=5;  %micron, convolution after calc

N_f=8192*2;

Max_Frequency=C/(Wavelength_Min*1E-6);

Wavelength=Data(:,1);           %nm
Frequency_Old=C./(Wavelength*1E-9);
Frequency=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency=Frequency';
Wavelength_Micron=(C./Frequency)*1E6;

Spectrum_Wavelength=Data(:,2)-Data(1,2);

Spectrum_Frequency=(Spectrum_Wavelength.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectrum=interp1(Frequency_Old,Spectrum_Frequency,Frequency); 
Spectrum(isnan(Spectrum))=0;

Spectrum_Normalized=Spectrum/max(Spectrum);

%%

Common_OPD=0;  %micron
M=0.2;

Wavelength_Max_Index=find(Wavelength_Micron<=Wavelength_Max,1,'first');

Wavelength_Min_Index=find(Wavelength_Micron<=Wavelength_Min,1,'first');

Spectrum_Normalized=Spectrum_Normalized(Wavelength_Max_Index:Wavelength_Min_Index);
Wavelength_Micron=Wavelength_Micron(Wavelength_Max_Index:Wavelength_Min_Index);

Pixel_size=5;          %micron
CCD_Half_Width=1024/2;      %pixel
CCD_Half_Height=768/2;    %pixel

plot(Wavelength_Micron,Spectrum_Normalized);

R_Reference=7.8;          %mm

R_Sample=7.8;           %mm

X=(1:CCD_Half_Width)*Pixel_size/M;
Y=(1:CCD_Half_Height)*Pixel_size/M;
One_X(1:length(X))=1;
One_Y(1:length(Y))=1;

X_Grid=(One_Y')*(X);
Y_Grid=(Y')*(One_X);

Distance_Reference=((R_Reference*1000)-(((R_Reference*1000)^2)-(X_Grid.^2)-((Y_Grid+5).^2)).^0.5)+Common_OPD;

%Distance_Reference=X_Grid/200;
Distance_Sample=(R_Sample*1000)-(((R_Sample*1000)^2)-(X_Grid.^2)-(Y_Grid.^2)).^0.5;
%Distance_Sample=(min(min(X_Grid))+max(max(X_Grid)))/400;

Inter=0;
for j=1:length(Wavelength_Micron)
    %Inter=Inter+(Spectrum_Normalized(j))+Spectrum_Normalized(j)*real(exp(i.*4.*pi./Wavelength_Micron(j).*Distance_Sample)./exp(i.*4.*pi./Wavelength_Micron(j).*Distance_Reference));
    Inter=Inter+(Spectrum_Normalized(j))+Spectrum_Normalized(j)*real(exp(i.*4.*pi./Wavelength_Micron(j).*(Distance_Sample-Distance_Reference)));
end

Inter_2=Inter(:,end:-1:1);
Inter_3=Inter(end:-1:1,end:-1:1);
Inter_4=Inter(end:-1:1,:);

Inter_full=[Inter_3 Inter_4;Inter_2 Inter];
filter_mactrx=ones(round(Lateral_Resolution/Pixel_size),round(Lateral_Resolution/Pixel_size));
Inter_full_Filtered=filter2(filter_mactrx,Inter_full,'same');
%FF=[1 0 -1; 1 0 -1; 1 0 -1];
%Inter_full=filter2(FF,Inter_full,'same');
imagesc(Inter_full_Filtered/max(max(Inter_full_Filtered)),'xdata',X,'ydata',Y);
xlabel('(Micron)');
ylabel('(Micron)');
colormap(gray);
axis equal