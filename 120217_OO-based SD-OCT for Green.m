%% Options

clear all


array=1:999;
lateral_index_reference_start=800;
lateral_index_reference_end=950;

lateral_index_sample_start=699;
lateral_index_sample_end=799;


range_specified=1;

Reference_Align=0;

Ratio=1;

Number_of_Loop=5;
Thickness=2.71.*1E-6;
Load_Previous_Result=0;
Position_Error=0.272*1E-6;%[0.019:0.001:0.024]*1E-6;%[0.02:0.0003:0.023]*1E-6;
n_should=1.7;
delta_n=0.1;
Wavelength_Center=540;

pixel_1=600;               % the end of 1
pixel_2=1800;%2000               % the end of 2
pixel_3=1120;               %for saperation

Wavelength_Considered_Min=480;          %nm
Wavelength_Considered_Max=600;

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=N_f*8;

Number_of_variable=3;           %Wavelength indep. variables


%% Global arrays generation

C=3E8;

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_Center=C/(Wavelength_Center*1E-9);
Frequency_Considered_Min=C/(Wavelength_Considered_Max*1E-9);             %Hz
Frequency_Considered_Max=C/(Wavelength_Considered_Min*1E-9);             %Hz


cd('D:\120213\5-1\');
Data=importdata(sprintf('D1.txt',1));      % Data_2: the glass data 1

Wavelength=Data(:,1);           %nm
Frequency_Old=C./(Wavelength*1E-9);
Frequency=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency=Frequency';
Wavelength_micron=(C./Frequency)*1E6;

Frequency_Center_Index=find(Frequency>Frequency_Center,1,'first');
Frequency_Considered_Min_Index=find(Frequency>Frequency_Considered_Min,1,'first');
Frequency_Considered_Max_Index=find(Frequency>Frequency_Considered_Max,1,'first');

% Time-domain

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;

%% Reference plane calibration and tilting angle generation
Profile(array)=0;
Profile_Calibration(array)=0;
Profile_New(array)=0;
for j=1:(lateral_index_reference_end-lateral_index_reference_start+1)

    cd('D:\120213\5-1\');
    Data_Reference=importdata(sprintf('D%i.txt',j+lateral_index_reference_start-1));      % Data_2: the glass data 1
    Spectrum_Reference_Old=Data_Reference(:,2);
    Spectrum_Reference_Frequency=(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
    Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency); 
    Spectrum_Reference(isnan(Spectrum_Reference))=0;
    Spectrum_Reference(Frequency<Min_Frequency)=0;
    Spectrum_Reference((N_f+1):N_t)=0;

    Signal_Reference=fft(Spectrum_Reference);
    Signal_Reference(1:pixel_1)=0;
    Signal_Reference((pixel_2+1):end)=0;
    [maxvalue_Reference maxindex_Reference]=max(abs(Signal_Reference));
    Profile(j+lateral_index_reference_start-1)=Position(maxindex_Reference);
    
    Spectrum_Reference=(ifft(Signal_Reference));
    Spectrum_Reference=2*(Spectrum_Reference(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform
    
    
    cd('D:\120213\5-R1\');
    Data_Reference_Calibration=importdata(sprintf('D%i.txt',j+lateral_index_reference_start-1));      % Data_2: the glass data 1
    Spectrum_Reference_Old_Calibration=Data_Reference_Calibration(:,2);
    Spectrum_Reference_Frequency_Calibration=(Spectrum_Reference_Old_Calibration.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
    Spectrum_Reference_Calibration=interp1(Frequency_Old,Spectrum_Reference_Frequency_Calibration,Frequency); 
    Spectrum_Reference_Calibration(isnan(Spectrum_Reference_Calibration))=0;
    Spectrum_Reference_Calibration(Frequency<Min_Frequency)=0;
    Spectrum_Reference_Calibration((N_f+1):N_t)=0;

    Signal_Reference_Calibration=fft(Spectrum_Reference_Calibration);
    Signal_Reference_Calibration(1:pixel_1)=0;
    Signal_Reference_Calibration((pixel_2+1):end)=0;
    [maxvalue_Reference_Calibration maxindex_Reference_Calibration]=max(abs(Signal_Reference_Calibration));
    Profile_Calibration(j+lateral_index_reference_start-1)=Position(maxindex_Reference_Calibration);
    
    Spectrum_Reference_Calibration=(ifft(Signal_Reference_Calibration));
    Spectrum_Reference_Calibration=2*(Spectrum_Reference_Calibration(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform
    
    
    Spectrum_Additional=exp(i*4*pi.*Frequency/C.*(Profile_Calibration(lateral_index_reference_start)-Profile_Calibration(j+lateral_index_reference_start-1)));        %總之都移到Profile_Reference(1)
    
    Spectrum_Reference_New=Spectrum_Reference.*Spectrum_Additional;
    Spectrum_Reference_New((N_f+1):N_t)=0;
    
    Signal_Reference_New=fft(Spectrum_Reference_New);
    Signal_Reference_New(1:pixel_1)=0;
    Signal_Reference_New((pixel_2+1):end)=0;
    [maxvalue_Reference_New maxindex_Reference_New]=max(abs(Signal_Reference_New));
    Profile_New(j+lateral_index_reference_start-1)=Position(maxindex_Reference_New);
    
end

% For tilting

A=polyfit(array(lateral_index_reference_start:lateral_index_reference_end),Profile_New(lateral_index_reference_start:lateral_index_reference_end),1);
Basedline=A(1)*array+A(2);
Profile_New(lateral_index_reference_start:lateral_index_reference_end)=Profile_New(lateral_index_reference_start:lateral_index_reference_end)-Basedline(lateral_index_reference_start:lateral_index_reference_end)+Basedline(lateral_index_reference_start);


%% Sample plane related calculation (at each position)
for j=1:(lateral_index_sample_end-lateral_index_sample_start+1)
    cd('D:\120213\5-1\');
    Data=importdata(sprintf('D%i.txt',j+lateral_index_sample_start-1));      % Data_2: the glass data 1
    Spectrum_Old=Data(:,2);
    Spectrum_Frequency=(Spectrum_Old.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
    Spectrum=interp1(Frequency_Old,Spectrum_Frequency,Frequency); 
    Spectrum(isnan(Spectrum))=0;
    Spectrum(Frequency<Min_Frequency)=0;
    Spectrum((N_f+1):N_t)=0;

    Signal=fft(Spectrum);
    Signal(1:pixel_1)=0;
    Signal((pixel_2+1):end)=0;
    [maxvalue maxindex]=max(abs(Signal));
    Profile(j+lateral_index_sample_start-1)=Position(maxindex);
    
    Spectrum=(ifft(Signal));
    Spectrum=2*(Spectrum(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform
    
    
    cd('D:\120213\5-R1\');
    Data_Calibration=importdata(sprintf('D%i.txt',j+lateral_index_sample_start-1));      % Data_2: the glass data 1
    Spectrum_Old_Calibration=Data_Calibration(:,2);
    Spectrum_Frequency_Calibration=(Spectrum_Old_Calibration.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
    Spectrum_Calibration=interp1(Frequency_Old,Spectrum_Frequency_Calibration,Frequency); 
    Spectrum_Calibration(isnan(Spectrum_Calibration))=0;
    Spectrum_Calibration(Frequency<Min_Frequency)=0;
    Spectrum_Calibration((N_f+1):N_t)=0;

    Signal_Calibration=fft(Spectrum_Calibration);
    Signal_Calibration(1:pixel_1)=0;
    Signal_Calibration((pixel_2+1):end)=0;
    [maxvalue_Calibration maxindex_Calibration]=max(abs(Signal_Calibration));
    Profile_Calibration(j+lateral_index_sample_start-1)=Position(maxindex_Calibration);   
    
    Spectrum_Calibration=(ifft(Signal_Calibration));
    Spectrum_Calibration=2*(Spectrum_Calibration(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform
    
    Spectrum_Additional=exp(i*4*pi.*Frequency/C.*(Profile_Calibration(lateral_index_reference_start)-Profile_Calibration(j+lateral_index_sample_start-1)));    %總之都移到Profile_Reference(1)
    Spectrum_Additional_Based=exp(i*4*pi.*Frequency/C.*(Basedline(lateral_index_reference_start)-Basedline(j+lateral_index_sample_start-1)));    
    Spectrum_New=Spectrum.*Spectrum_Additional.*Spectrum_Additional_Based;
    Spectrum_New((N_f+1):N_t)=0;
    
    Signal_New=fft(Spectrum_New);
    Signal_New(1:pixel_1)=0;
    Signal_New((pixel_2+1):end)=0;
    [maxvalue_New maxindex_New]=max(abs(Signal_New));
    Profile_New(j+lateral_index_sample_start-1)=Position(maxindex_New);

end



%Data_Spectroscopy=importdata('111118_Green 5-1.jws.txt');

Spectrum_Old=Data(:,2);
Spectrum_Reference_Old=Data_Reference(:,2);

Spectroscopy_Old=Data_Spectroscopy(end:-1:1,2)/100;
Wavelength_Spectroscopy=Data_Spectroscopy(end:-1:1,1);           %nm

Frequency_Spectroscopy_Old=C./(Wavelength_Spectroscopy*1E-9);
Spectrum_Frequency=(Spectrum_Old.*((Wavelength*1E-9).^2)/C);        %/max(Spectrum_Old.*((Wavelength*1E-9).^2)/C);  N more normalization
Spectrum_Reference_Frequency=(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectroscopy_Frequency=Spectroscopy_Old;                        %NOTE!

Frequency_Spectroscopy_Considered_Min_Index=find(Frequency_Spectroscopy_Old>Frequency_Considered_Min,1,'first');
Frequency_Spectroscopy_Considered_Max_Index=find(Frequency_Spectroscopy_Old>Frequency_Considered_Max,1,'first');


Spectrum=interp1(Frequency_Old,Spectrum_Frequency,Frequency);
Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency);
Spectroscopy=interp1(Frequency_Spectroscopy_Old,Spectroscopy_Frequency,Frequency);

Spectrum(isnan(Spectrum))=0;
Spectrum(Frequency<Min_Frequency)=0;

Spectrum_Reference(isnan(Spectrum))=0;
Spectrum_Reference(Frequency<Min_Frequency)=0;


Spectroscopy(isnan(Spectroscopy))=0;
Spectroscopy(Frequency<Min_Frequency)=0;


%% Spatial Filtering 

Spectrum((N_f+1):N_t)=0;
Spectrum_Reference((N_f+1):N_t)=0;

Signal=fft(Spectrum);
Signal_Reference=fft(Spectrum_Reference);

%% Data separation (Spatial Filtering)

Signal_Upper=Signal;
Signal_Lower=Signal;

%Signal_1=Signal;        %DC+self interference
%Signal_2=Signal;        %Upper
%Signal_3=Signal;        %Lower

%pixel_3=1600;               % the end of 3

%Signal_1((pixel_1+1):(length(Signal_1)-pixel_1))=0;

%Signal((pixel_2+1):(length(Signal)-pixel_2))=0;
Signal(1:pixel_1)=0;
Signal((pixel_2+1):(length(Signal)-pixel_2))=0;
Signal((length(Signal)-pixel_1+1):end)=0;

Signal((pixel_2+1):end)=0;


Signal_Upper(1:pixel_1)=0;
Signal_Upper((pixel_3+1):(length(Signal)-pixel_3))=0;
Signal_Upper((length(Signal)-pixel_1+1):end)=0;

Signal_Upper((pixel_3+1):end)=0;

Signal_Lower(1:pixel_3)=0;
Signal_Lower((pixel_2+1):(length(Signal)-pixel_2))=0;
Signal_Lower((length(Signal)-pixel_3+1):end)=0;

Signal_Lower((pixel_2+1):end)=0;

Signal_Reference(1:pixel_1)=0;
Signal_Reference((pixel_2+1):end)=0;
%Signal(((length(Signal)-pixel_1)+1):end)=0;

%Signal_3((pixel_3+1):(length(Signal_3)-pixel_3))=0;
%Signal_3(1:pixel_2)=0;
%Signal_3(((length(Signal_3)-pixel_2)+1):end)=0;

%plot(Position,Signal_1,Position,Signal_2,Position,Signal_3);
plot(Position,abs(Signal),Position,abs(Signal_Reference));



%% Again FD

%Spectrum_1=real(ifft(Signal_1));       %can take Real, since the imaginary part should be relatively small, and it came from error
%Spectrum_2=real(ifft(Signal_2));
%Spectrum_3=real(ifft(Signal_3));

Spectrum=(ifft(Signal));
Spectrum=2*(Spectrum(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform


Spectrum_Upper=(ifft(Signal_Upper));
Spectrum_Upper=2*(Spectrum_Upper(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform

Spectrum_Lower=(ifft(Signal_Lower));
Spectrum_Lower=2*(Spectrum_Lower(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform

Spectrum_Reference=(ifft(Signal_Reference));
Spectrum_Reference=2*(Spectrum_Reference(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform

%Spectrum_1=Spectrum_1(1:N_f);
%Spectrum_2=Spectrum_2(1:N_f);
%Spectrum_3=Spectrum_3(1:N_f);

%% Divided by Spectrum_Reference

[maxvalue maxindex]=max(abs(Signal(1:pixel_3)));
Position_1=Position(maxindex);
    
    %Thickness:the peak position of the second interface
[maxvalue maxindex]=max(abs(Signal_Reference));
Position_0=Position(maxindex);

Spectrum_Reference_New=Spectrum_Reference;
Spectrum_Reference_New=Spectrum_Reference_New.*exp(i*4*pi.*Frequency/C.*(Position_1-Position_0)); 
Spectrum_Reference_New((N_f+1):N_t)=0;
Signal_Reference_New=fft(Spectrum_Reference_New/2);

if Reference_Align == 1
    
    Spectrum_Reference=Spectrum_Reference_New(1:N_f);
end

plot(Position,abs(Signal),Position,abs(Signal_Reference),Position,abs(Signal_Reference_New));

dlmwrite('Position.txt',Position,'delimiter','\t','newline','pc','precision','%.12f');

dlmwrite('Signal_env.txt',abs(Signal)/max(abs(Signal_Reference_New)),'delimiter','\t','newline','pc','precision','%.12f');

dlmwrite('Signal_Reference_env.txt',abs(Signal_Reference_New)/max(abs(Signal_Reference_New)),'delimiter','\t','newline','pc','precision','%.12f');

dlmwrite('Signal_car.txt',real(Signal)/max(abs(Signal_Reference_New)),'delimiter','\t','newline','pc','precision','%.12f');

dlmwrite('Signal_Reference_car.txt',real(Signal_Reference_New)/max(abs(Signal_Reference_New)),'delimiter','\t','newline','pc','precision','%.12f');
%Thickness=2.468E-6;

% Additional Shift
Spectrum=Spectrum./Spectrum_Reference;%.*Spectrum_Additional;


%Spectrum((N_f+1):N_t)=0;
%Signal_Again=fft(Spectrum);
%Signal_Again(1:pixel_1)=0;
%Signal_Again((pixel_2+1):(length(Signal_Again)-pixel_2))=0;
%Signal_Again((length(Signal_Again)-pixel_1+1):end)=0;
%Spectrum=(ifft(Signal));
%Spectrum=Spectrum(1:N_f);
%Spectrum_Ratio=Spectrum_Lower./Spectrum_Upper;
%Spectrum=Spectrum-Spectrum_Reference;
%Spectrum=Spectrum-Spectrum_Sample;
plot(Wavelength_micron,Spectrum);
%% Theory - Sample Model (n1 - n - n2)

n2=1;

% Assumed n2 = BK7
C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;


n_bk7=(C1*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C2)+C3*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C4)+C5*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C6)+1).^0.5;

n_bk7=abs(n_bk7);

n_bk7(isnan(n_bk7))=0;
n1=n_bk7;

%% To generate the original spectrum, assuming reference is BK7, too

r_BK7=((n_bk7-1)./(n_bk7+1));
t_AN100=(2*(n_bk7)./(n_bk7+1));
%Spectrum_Original=Spectrum_Reference./(r_BK7).^2;

%% Fitting

n_final(1:length(Frequency),1:length(n_should))=1.5;
k_final(1:length(Frequency),1:length(n_should))=0;
Spectroscopy_Amplitude_Check(1:length(Spectrum),1:length(n_should))=0;
T_check(1:length(Spectrum),1:length(n_should))=0;

Spectrum_Check(1:length(Spectrum),1:length(n_should))=0;
Spectrum_Upper_Check(1:length(Spectrum),1:length(n_should))=0;
Spectrum_Lower_Check(1:length(Spectrum),1:length(n_should))=0;
Spectrum_Lower_Check_n(1:length(Spectrum),1:length(n_should))=0;
Spectrum_Lower_Check_k(1:length(Spectrum),1:length(n_should))=0;


for q=1:length(n_should)
Spectrum_Additional=exp(i*4*pi.*Frequency/C.*(Position_Error));
value_temp=10000000000;
n(1:length(Frequency),1)=1.5;
k(1:length(Frequency),1)=0;

%T_max_check(1:length(Spectrum),1:length(n_should))=0;

value_total_1(1:Number_of_Loop)=0;
value_total_2(1:Number_of_Loop)=0;

    n_o=(n_should-delta_n):0.0001:(n_should(q)+delta_n);
    k_o=0:0.0003:0.3;
    %k_o=k_o';
    delta_n_number=round(delta_n/(n_o(2)-n_o(1)))/2;
    %n_empty=n_o;
    %n_empty(:)=1;
    %k_empty=k_o;
    %k_empty(:)=1;
    %n_temp=k_empty*n_o;
    %k_temp=k_o*n_empty;
    n_temp=n_o;
    k_temp=k_o;
    
    n_should_index=find(n_o-n_should(q)>0,1,'first');
    
for p=1:Number_of_Loop
    
    % n Template Generation

    % For OCT    
% For Spectroscopy    
    %Spectrum_Amplitude_Check(1:length(Spectrum))=0;
    index_old=0;
    d0_Check(1:length(Spectrum))=0;
    dref_Check(1:length(Spectrum))=0;
    
    for j=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
            
        if j <= ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)
            index_now=Frequency_Center_Index+j-1;
        elseif j > ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)           %direction reverse!!
            index_now=-j+((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)+Frequency_Center_Index;
        end
            
        r1=(n1(index_now)-(n_temp+i*k(index_now)))./(n_temp+n1(index_now)+i*k(index_now));
        r1_r=((n_temp+i*k(index_now))-n1(index_now))./(n_temp+n1(index_now)+i*k(index_now));    
        t1=2*(n1(index_now))./(n_temp+n1(index_now)+i*k(index_now));
        t1_r=2*(n_temp+i*k(index_now))./(n_temp+n1(index_now)+i*k(index_now));
        t2=2*(n_temp+i*k(index_now))./(n_temp+n2+i*k(index_now));
        r2=((n_temp+i*k(index_now))-n2)./((n_temp+i*k(index_now))+n2);   
        d=exp(i*2*pi.*Frequency(index_now)/C.*(n_temp+i*k(index_now)).*Thickness);   %注意! -1*n!
        %d0=exp(i*2*pi.*Frequency(index_now)/C.*Thickness_0);   %注意! -1*n!
        %dref=exp(i*2*pi.*Frequency(index_now)/C.*(Thickness_0+Thickness));   %注意! -1*n!
        
        
        Spectrum_Upper_Temp=((r1)./r_BK7(index_now));                       %神說是cosine
        
        Spectrum_Lower_Temp=((t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now))*Ratio;                       %神說是cosine
        
        %Spectrum_Temp=((r1+t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now));                       %神說是cosine
        Spectrum_Temp=(Spectrum_Upper_Temp+Spectrum_Lower_Temp)*Spectrum_Additional(index_now);
        
        Merit=((Spectrum_Temp-(Spectrum(index_now))).^2); %(abs(abs(Spectrum_Temp)-abs(Spectrum(index_now))).^2)+
          
        if range_specified == 1
            if (index_now == Frequency_Center_Index) || (index_now == (Frequency_Center_Index -1))
                range_upper=min(n_should_index+delta_n_number,length(Merit));
                range_lower=max(n_should_index-delta_n_number,1);
            elseif (index_now > Frequency_Center_Index) || (index_now < (Frequency_Center_Index -1))
                range_upper=min(index_old+delta_n_number,length(Merit));
                range_lower=max(index_old-delta_n_number,1);
            end
            [value index]=min(Merit(range_lower:range_upper));
            index=range_lower+index-1;
            index_old=index;
        else
            [value index]=min(Merit);
        end
        
        value_total_1(p)=value_total_1(p)+abs(value);
        n(index_now)=n_temp(index);    %index1,index2=k,n  in situ saving all the solutions
            
    end
    
    for j=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
            
        if j <= ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)
            index_now=Frequency_Center_Index+j-1;
        elseif j > ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)           %direction reverse!!
            index_now=-j+((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)+Frequency_Center_Index;
        end
            
            % About the conditions
        
        r1=(n1(index_now)-(n(index_now)+i*k_temp))./(n(index_now)+n1(index_now)+i*k_temp);
        r1_r=((n(index_now)+i*k_temp)-n1(index_now))./(n(index_now)+n1(index_now)+i*k_temp);    
        t1=2*(n1(index_now))./(n(index_now)+n1(index_now)+i*k_temp);
        t1_r=2*(n(index_now)+i*k_temp)./(n(index_now)+n1(index_now)+i*k_temp);
        t2=2*(n(index_now)+i*k_temp)./(n(index_now)+n2+i*k_temp);
        r2=((n(index_now)+i*k_temp)-n2)./((n(index_now)+i*k_temp)+n2);   
        d=exp(i*2*pi.*Frequency(index_now)/C.*(n(index_now)+i*k_temp).*Thickness);   %注意! -1*n!
        
        d_n=exp(i*2*pi.*Frequency(index_now)/C.*(n(index_now)).*Thickness);   %注意! -1*n!
        d_k=exp(i*2*pi.*Frequency(index_now)/C.*(i*k_temp).*Thickness);   %注意! -1*n!
        %d0=exp(i*2*pi.*Frequency(index_now)/C.*Thickness_0);   %注意! -1*n!
        %dref=exp(i*2*pi.*Frequency(index_now)/C.*(Thickness_0+Thickness));   %注意! -1*n!
        Spectrum_Upper_Temp=((r1)./r_BK7(index_now));                       %神說是cosine
        
        Spectrum_Lower_Temp=((t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now))*Ratio;                       %神說是cosine
        Spectrum_Lower_Temp_n=((t1.*t1_r.*r2.*(d_n.^2))./r_BK7(index_now))*Ratio;                       %神說是cosine
        Spectrum_Lower_Temp_k=((t1.*t1_r.*r2.*(d_k.^2))./r_BK7(index_now))*Ratio;                       %神說是cosine        
        %Spectrum_Temp=((r1+t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now));                       %神說是cosine
        Spectrum_Temp=(Spectrum_Upper_Temp+Spectrum_Lower_Temp)*Spectrum_Additional(index_now);
        %Spectroscopy_Temp=abs((t_AN100(index_now)).*(t1.*t2.*d)./(1+r1_r.*r2.*(d.^2))).^2;
        %T_Temp=abs((t_AN100(index_now).^2).*(t1.*t2.*d)./(1+r1_r.*r2.*(d.^2))).^2;
        Merit=((Spectrum_Temp-(Spectrum(index_now))).^2);%+((Spectroscopy_Temp-(Spectroscopy(index_now))).^2); %(abs(abs(Spectrum_Temp)-abs(Spectrum(index_now))).^2)+
          

        [value index]=min(Merit);
        value_total_2(p)=value_total_2(p)+abs(value);
        Spectrum_Check(index_now,q)=Spectrum_Temp(index);
        Spectrum_Upper_Check(index_now,q)=Spectrum_Upper_Temp(index);
        Spectrum_Lower_Check(index_now,q)=Spectrum_Lower_Temp(index);
        Spectrum_Lower_Check_n(index_now,q)=Spectrum_Lower_Temp_n(index);
        Spectrum_Lower_Check_k(index_now,q)=Spectrum_Lower_Temp_k(index);
        k(index_now)=k_temp(index);
            
    end
   
end

    n_final(:,q)=n;
    k_final(:,q)=k;

end

%% Checking

%plot(Wavelength_micron,real(Spectrum),Wavelength_micron,real(Spectrum_Check),Wavelength_micron,real(Spectrum_Upper_Check),Wavelength_micron,real(Spectrum_Lower_Check));

%plot(Wavelength_micron,real(Spectrum),Wavelength_micron,real(Spectrum_Check),Wavelength_micron,real(Spectrum_Upper_Check),Wavelength_micron,real(Spectrum_Lower_Check),Wavelength_micron,real(Spectrum_Lower_Check_n),Wavelength_micron,real(Spectrum_Lower_Check_k));
%plot(Wavelength_micron,abs(Spectrum_Check(:,5)),Wavelength_micron,abs(Spectrum_Upper_Check(:,5)+Spectrum_Lower_Check_n(:,5)),Wavelength_micron,abs(Spectrum_Upper_Check(:,5)+Spectrum_Lower_Check_k(:,5)));
%plot(Wavelength_micron,abs(Spectrum_Check(:,5)),Wavelength_micron,real(Spectrum_Check(:,5)),Wavelength_micron,real(Spectrum_Upper_Check+Spectrum_Lower_Check_n(:,5)),Wavelength_micron,real(Spectrum_Upper_Check+Spectrum_Lower_Check_k(:,5)),Wavelength_micron,abs(Spectrum_Upper_Check+Spectrum_Lower_Check_n(:,5)),Wavelength_micron,abs(Spectrum_Upper_Check+Spectrum_Lower_Check_k(:,5)));

%plot(Wavelength_micron,angle(Spectrum),Wavelength_micron,angle(Spectrum_Check));
%plot(Wavelength_micron,Spectrum_Check,Wavelength_micron,Spectrum);

%plot(Wavelength_micron,Spectroscopy_Check,Wavelength_micron,Spectroscopy);
plot(Wavelength_micron,n_final);
plot(Wavelength_micron(Wavelength_micron<0.6),k_final(Wavelength_micron<0.6,:));
%plot(Wavelength_micron(Wavelength_micron<0.6),k_final(Wavelength_micron<0.6,5));
%plot(1:Number_of_Loop,value_total_1,1:Number_of_Loop,value_total_2);
dlmwrite('n_should.txt',n_should,'delimiter','\t','newline','pc','precision','%.12f');

dlmwrite('n_final.txt',n_final,'delimiter','\t','newline','pc','precision','%.12f');

dlmwrite('k_final.txt',k_final,'delimiter','\t','newline','pc','precision','%.12f');

dlmwrite('Spectrum_abs.txt',abs(Spectrum),'delimiter','\t','newline','pc','precision','%.12f');

dlmwrite('Spectrum_real.txt',real(Spectrum),'delimiter','\t','newline','pc','precision','%.12f');

dlmwrite('Frequency.txt',Frequency,'delimiter','\t','newline','pc','precision','%.12f');

dlmwrite('Wavelength_micron.txt',Wavelength_micron,'delimiter','\t','newline','pc','precision','%.12f');
