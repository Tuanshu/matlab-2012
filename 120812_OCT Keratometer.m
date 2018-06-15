clear all
%%

cd('D:\Users\TuanShu');
Data=importdata('130319_YbFiber(Clad)_Spectrum_dc.txt');      % Data_2: the glass data 1

C=3E8;

Wavelength_Max=1.1; %micron
Wavelength_Min=0.98;

Lateral_Resolution=5;  %micron, convolution after calc

N_f=8192*4;

Sise_Diameter=11000;     %micron

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

Common_OPD=100;  %micron
M=0.3;

Wavelength_Max_Index=find(Wavelength_Micron<=Wavelength_Max,1,'first');

Wavelength_Min_Index=find(Wavelength_Micron<=Wavelength_Min,1,'first');

Spectrum_Normalized=Spectrum_Normalized(Wavelength_Max_Index:Wavelength_Min_Index);
Wavelength_Micron=Wavelength_Micron(Wavelength_Max_Index:Wavelength_Min_Index);

Pixel_size=5;          %micron
CCD_Half_Width=1024/2;      %pixel
CCD_Half_Height=768/2;    %pixel

plot(Wavelength_Micron,Spectrum_Normalized);

R_Reference=999999;%7.79;          %mm

R_Sample=7.8;%7.8;           %mm

X=(1:CCD_Half_Width)*Pixel_size/M;
Y=(1:CCD_Half_Height)*Pixel_size/M;
One_X(1:length(X))=1;
One_Y(1:length(Y))=1;

X_Grid=(One_Y')*(X);
Y_Grid=(Y')*(One_X);

Distance_Reference=((R_Reference*1000)-(((R_Reference*1000)^2)-(X_Grid.^2)-((Y_Grid).^2)).^0.5)+Common_OPD;

%Distance_Reference=X_Grid/200;
Distance_Sample=(R_Sample*1000)-(((R_Sample*1000)^2)-(X_Grid.^2)-(Y_Grid.^2)).^0.5;

%Distance_Sample=(min(min(X_Grid))+max(max(X_Grid)))/400;
Mask(1:size(Distance_Sample,1),1:size(Distance_Sample,2))=1;
Mask(((X_Grid.^2)+(Y_Grid.^2))>(Sise_Diameter/2).^2)=0;
Distance_Sample=real(Distance_Sample);


Inter=0;
for j=1:length(Wavelength_Micron)
    %Inter=Inter+(Spectrum_Normalized(j))+Spectrum_Normalized(j)*real(exp(i.*4.*pi./Wavelength_Micron(j).*Distance_Sample)./exp(i.*4.*pi./Wavelength_Micron(j).*Distance_Reference));
    Inter=Inter+(Spectrum_Normalized(j))+Spectrum_Normalized(j)*(exp(i.*4.*pi./Wavelength_Micron(j).*(Distance_Sample-Distance_Reference)));
end
Inter=real(Inter);

K=6;
XXX=Distance_Sample(K,:)-Distance_Reference(K,:);
YYY=Inter(K,:);

%plot(YYY);

Inter=Inter.*Mask;
Inter(isnan(Inter))=0;
Inter_2=Inter(:,end:-1:1);
Inter_3=Inter(end:-1:1,end:-1:1);
Inter_4=Inter(end:-1:1,:);

Inter_full=[Inter_3 Inter_4;Inter_2 Inter];
Inter_full_Normalized=Inter_full/max(max(Inter_full));
%filter_mactrx=ones(round(Lateral_Resolution/Pixel_size),round(Lateral_Resolution/Pixel_size));
%Inter_full_Filtered=filter2(filter_mactrx,Inter_full,'same');
%FF=[1 0 -1; 1 0 -1; 1 0 -1];
%Inter_full=filter2(FF,Inter_full,'same');
imagesc(Inter_full_Normalized,'xdata',X,'ydata',Y);
[max_value max_index]=max(Inter);
xlabel('(Micron)');
ylabel('(Micron)');
colormap(gray);
caxis([0 1]);
axis equal

        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        


%%

        Lowpass=20;
        Longpass=250;
        Image_Temp=Inter;        
       
        [r c] = size(Image_Temp);
        [X Y] = meshgrid(1:c,1:r);
        X=X-1;
        Y=Y-1;
        NNN=360;
        R_max=150;
        r_wish=repmat([R_max/NNN:R_max/NNN:R_max],NNN,1);
        theta_wish=repmat([pi/2/NNN:pi/2/NNN:pi/2]',1,NNN);
        
        X_wish=r_wish.*cos(theta_wish);
        Y_wish=r_wish.*sin(theta_wish);

        image_Rtheta=interp2(X,Y,(Image_Temp),X_wish,Y_wish,'linear');
        
        imagesc(image_Rtheta);
        image_Rtheta=image_Rtheta-mean(mean(image_Rtheta(:,240:250)));

        FFT_R=fft(image_Rtheta,[],2);
        FFT_R(:,size(FFT_R,2):end)=0;    
        image_Rtheta_New=abs(ifft(FFT_R,[],2));
                imagesc(image_Rtheta_New);

xlabel('Radial Position (micron)');
ylabel('Angle (degree)');

        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        
        image_Rtheta_New(:,320:end)=0;

       
        X_grid=repmat(1:768,768,1);
        Y_grid=repmat([1:768]',1,768);
        
        [r_rt c_rt] = size(Image_Temp);
        [X_rt Y_rt] = meshgrid(1:c_rt,1:r_rt);
        
        r_XY=((X).^2+(Y).^2).^0.5;
        theta_XY=atan2((Y),(X))+pi;
        
        image_XY_New=interp2(r_wish,theta_wish,image_Rtheta_New,r_XY,theta_XY);
        
        
        imagesc(image_XY_New);
        

        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
