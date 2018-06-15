%% Options
%% to try with modeling
clear all

center_wavelength_spectral_shapping=625;    %nm
bandwidth_spectral_shapping=20;

index_calculation=0;
index_initial_boundary_reference=2150;
index_initial_boundary_sample=2150;

array=1:999;

index_used_reference_plane=800;

index_edge_range_left=760;

index_edge_range_right=780;


lateral_index_calibration_start=550;
lateral_index_calibration_end=950;

lateral_index_reference_start=800;
lateral_index_reference_end=950;

lateral_index_sample_start=550;
lateral_index_sample_end=550;


range_specified=1;

Reference_Align=0;

Ratio=1;
Ratio2=1;
Number_of_Loop=5;
Thickness=2.4.*1E-6;
Load_Previous_Result=0;
Position_Error=0E-6;%[0.019:0.001:0.024]*1E-6;%[0.02:0.0003:0.023]*1E-6;

A0=3;     %roughly n0^2


%%
n_should=1.7;
delta_n=0.1;
Wavelength_Center=600;

pixel_1=1400;               % DC cutoff
%pixel_2=1800;%2000               % the end of 2
%pixel_3=1120;               %for saperation

Wavelength_Considered_Min=560;          %nm
Wavelength_Considered_Max=640;

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=N_f*4;

Number_of_variable=3;           %Wavelength indep. variables


%% Global arrays generation

C=3E8;

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_Center=C/(Wavelength_Center*1E-9);
Frequency_Considered_Min=C/(Wavelength_Considered_Max*1E-9);             %Hz
Frequency_Considered_Max=C/(Wavelength_Considered_Min*1E-9);             %Hz


cd('D:\120222\R2-3\');
Data=importdata('D1.txt');      % Data_2: the glass data 1

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

%% Spectral shapping function genreation

center_frequency_spectral_shapping=C/(600*1E-9);    %nm
bandwidth_frequency_spectral_shapping=(bandwidth_spectral_shapping*1E-9)*C/(center_wavelength_spectral_shapping*1E-9)^2;

spectral_shapping_function=gaussmf(Frequency,[bandwidth_frequency_spectral_shapping center_frequency_spectral_shapping]);

%% Theory - Sample Model (n1 - n - n2)

n2=1;

% Assumed n1 = BK7
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
%% Calibration plane - only use its maximum position

Profile_Sample_Reference(array)=0;
Profile_Calibration_Shifted(array)=0;
Profile_Calibration(array)=0;
index_Profile_Calibration(array)=0;
Distance(array)=0;
index_curernt_boundary=index_initial_boundary_reference;
DD=0;
for j=1:(lateral_index_calibration_end-lateral_index_calibration_start+1)
    Data=importdata(sprintf('D%i.txt',j+lateral_index_calibration_start-1));      % Data_2: the glass data 1
    Spectrum_Old=Data(:,2)-mean(Data(1:1000,2));
    Spectrum_Frequency=(Spectrum_Old.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
    Spectrum=interp1(Frequency_Old,Spectrum_Frequency,Frequency); 
    Spectrum(isnan(Spectrum))=0;
    Spectrum=Spectrum.*spectral_shapping_function;
    Spectrum(Frequency<Min_Frequency)=0;
    Spectrum((N_f+1):N_t)=0;
    Signal=fft(Spectrum);
    Signal(1:pixel_1)=0;
    Signal(round(length(Signal)/2):end)=0;
    
    Signal_Sample_Reference=Signal;
    Signal_Sample_Reference((index_curernt_boundary+1):end)=0;
    
    %Signal_Sample_Reference(1:2700)=0;
    
    %Signal_Sample_Reference(3200:end)=0;
    
    [maxvalue maxindex_sample_reference]=max(abs(Signal_Sample_Reference));
    Profile_Sample_Reference(j+lateral_index_calibration_start-1)=Position(maxindex_sample_reference);

    Signal_Calibration=Signal;
    Signal_Calibration(1:index_curernt_boundary)=0;
    
    [maxvalue maxindex_calibration]=max(abs(Signal_Calibration));
    
    if abs((maxindex_sample_reference+maxindex_calibration)/2-index_curernt_boundary) < 100
        index_curernt_boundary=round((maxindex_sample_reference+maxindex_calibration)/2);
    end 

    if j == 1
        maxindex_calibration_previous=maxindex_calibration;
    end
    
    if abs(maxindex_calibration_previous-(maxindex_calibration))>1000
        maxindex_calibration=maxindex_calibration_previous;
    elseif abs(maxindex_calibration_previous+DD-(maxindex_calibration))>25 && ((j+lateral_index_calibration_start-1) >index_edge_range_left) && ((j+lateral_index_calibration_start-1) <index_edge_range_right)
        DD=(maxindex_calibration_previous+DD)-(maxindex_calibration);
        edge_index=lateral_index_calibration_start+j-1;
    end

    Profile_Calibration_Shifted(j+lateral_index_calibration_start-1)=Position(maxindex_calibration+DD);
    Profile_Calibration(j+lateral_index_calibration_start-1)=Position(maxindex_calibration);
    index_Profile_Calibration(j+lateral_index_calibration_start-1)=maxindex_calibration;
    
    maxindex_calibration_previous=maxindex_calibration;
        %Spectrum_Additional=exp(i*4*pi.*Frequency/C.*(Profile_Calibration(
        %lateral_index_reference_start)-Profile_Calibration(j+lateral_index_reference_start-1)));        %總之都移到Profile_Reference(1)
        
    if (j+lateral_index_calibration_start-1) ==index_used_reference_plane
        index_used_reference_boundary=index_curernt_boundary;
    end
        
end
    
%% Reference plane calibration and tilting angle generation

A=polyfit(array(lateral_index_reference_start:lateral_index_reference_end),Profile_Calibration_Shifted(lateral_index_reference_start:lateral_index_reference_end)-Profile_Sample_Reference(lateral_index_reference_start:lateral_index_reference_end),1);
Distance_fit=A(1)*array+A(2);   %the distance between upper interface and additional interface!!
%Profile_New(lateral_index_reference_start:lateral_index_reference_end)=Profile_New(lateral_index_reference_start:lateral_index_reference_end)-Basedline(lateral_index_reference_start:lateral_index_reference_end)+Basedline(lateral_index_reference_start);
Profile_Sample_Reference=Profile_Sample_Reference-Profile_Calibration_Shifted+Profile_Calibration_Shifted(lateral_index_reference_start)+Distance_fit;

%% Generate the Reference Spectrum
Profile_New(array)=0;

Data=importdata(sprintf('D%i.txt',index_used_reference_plane));      % Data_2: the glass data 1
Spectrum_Old=Data(:,2)-mean(Data(1:1000,2));
Spectrum_Frequency=(Spectrum_Old.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectrum=interp1(Frequency_Old,Spectrum_Frequency,Frequency); 
Spectrum(isnan(Spectrum))=0;
Spectrum=Spectrum.*spectral_shapping_function;
Spectrum(Frequency<Min_Frequency)=0;
Spectrum((N_f+1):N_t)=0;

Signal=fft(Spectrum);
Signal(1:pixel_1)=0;
Signal(round(length(Signal)/2):end)=0;
    
Signal_Sample_Reference=Signal;
Signal_Sample_Reference((index_used_reference_boundary+1):end)=0;
    
Spectrum_Used_Reference=ifft(Signal_Sample_Reference);
Spectrum_Used_Reference=2*Spectrum_Used_Reference(1:N_f);
        
%Spectrum_Additional=exp(i*4*pi.*Frequency/C.*(Profile_Calibration_Shifted(lateral_index_reference_start)-Profile_Calibration_Shifted(index_used_reference_plane)));    %總之都移到Profile_Reference(1)
%Spectrum_Additional_Based=exp(i*4*pi.*Frequency/C.*(Distance_fit(index_used_reference_plane)));      

Spectrum_Additional_Error=exp(i*4*pi.*Frequency/C.*(Position_Error));    

Spectrum_Used_Reference=Spectrum_Used_Reference;%.*Spectrum_Additional.*Spectrum_Additional_Based.*Spectrum_Additional_Error;  

%Spectrum_Used_Reference(N_f+1:N_t)=0;

%Signal_Used_Reference=fft(Spectrum_Used_Reference);
%[maxvalue_Used_Reference maxindex_Used_Reference]=max(abs(Signal_Used_Reference));    

%Profile_New(index_used_reference_plane)=Position(maxindex_Used_Reference);

%% Sample plane related calculation (at each position)
index_curernt_boundary=index_initial_boundary_sample;
Profile_Old(array)=0;
for j=1:(lateral_index_sample_end-lateral_index_sample_start+1)
    Data=importdata(sprintf('D%i.txt',j+lateral_index_sample_start-1));      % Data_2: the glass data 1
    Spectrum_Old=Data(:,2)-mean(Data(1:1000,2));
    Spectrum_Frequency=(Spectrum_Old.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
    Spectrum=interp1(Frequency_Old,Spectrum_Frequency,Frequency); 
    Spectrum(isnan(Spectrum))=0;
    Spectrum=Spectrum.*spectral_shapping_function;
    Spectrum(Frequency<Min_Frequency)=0;
    Spectrum_Ori=Spectrum;
    Spectrum((N_f+1):N_t)=0;

    Signal=fft(Spectrum);
    Signal(1:pixel_1)=0;
    Signal(round(length(Signal)/2):end)=0;
    [maxvalue maxindex]=max(abs(Signal(1:index_curernt_boundary)));
    Profile(j+lateral_index_sample_start-1)=Position(maxindex);
    
    Spectrum=(ifft(Signal));
    Spectrum=2*(Spectrum(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform
    
    Signal_Calibration=Signal;
    Signal_Calibration(1:index_curernt_boundary)=0;
    Spectrum_Calibration=(ifft(Signal_Calibration));
    
    Spectrum_Calibration=2*(Spectrum_Calibration(1:N_f));  
    
    [maxvalue_Calibration maxindex_Calibration]=max(abs(Signal_Calibration));  
    
    if abs((maxindex+maxindex_Calibration)/2-index_curernt_boundary) < 100
        index_curernt_boundary=round((maxindex+maxindex_Calibration)/2);
    end 
    
    Signal_Sample=Signal;
    Signal_Sample((index_curernt_boundary+1):end)=0;
    Spectrum_Sample=(ifft(Signal_Sample));
    Spectrum_Sample=Spectrum_Sample(1:N_f);
    
    %% This part is only for drawing, for the real calculation, only
    %% Spectrum_Additional_based should be used, since the
    %% Spectrum_Additional term will be null in ther afterward calculation 
    %Spectrum_Additional=exp(i*4*pi.*Frequency/C.*(Profile_Calibration_Shifted(lateral_index_reference_start)-Profile_Calibration_Shifted(j+lateral_index_sample_start-1)));    %總之都移到Profile_Reference(1)
    %Spectrum_Additional_Based=exp(i*4*pi.*Frequency/C.*(Distance_fit(j+lateral_index_sample_start-1)));    
    %Spectrum_New=Spectrum.*Spectrum_Additional.*Spectrum_Additional_Based;
    
    %Spectrum_Additional=exp(i*4*pi.*Frequency/C.*(Profile_Calibration(lateral_index_reference_start)-Profile_Calibration(j+lateral_index_sample_start-1)));    %總之都移到Profile_Reference(1)
    %Spectrum_Additional_Based=exp(i*4*pi.*Frequency/C.*(-Basedline(j+lateral_index_sample_start-1)));    
    %Spectrum_Sample_New=Spectrum.*Spectrum_Additional_Based;
    %Spectrum_Calibration_New=Spectrum_Calibration.*Spectrum_Additional_Based;
    %Spectrum_New=Spectrum_Sample;%.*Spectrum_Additional.*Spectrum_Additional_Based;
    %Spectrum_New((N_f+1):N_t)=0;
    %Signal_New=fft(Spectrum_New);
    %Signal_New(1:pixel_1)=0;
    %Signal_New(round(length(Signal_New)/2):end)=0;
    %B_scan(:,j)=Signal_New(1:round(length(Signal_New)/4));
    %[maxvalue_New maxindex_New]=max(abs(Signal_New));    
%    [maxvalue_Old maxindex_Old]=max(abs(Signal_Sample));
    %Profile_New(j+lateral_index_sample_start-1)=Position(maxindex_New);
    %Spectrum_New=Spectrum_New(1:N_f);
%    Profile_Old(j+lateral_index_sample_start-1)=Position(maxindex_Old);
    
    
%% Filtering the Merit Function
%Spectrum_Used_Reference_Shifted=Spectrum_Used_Reference.*exp(i*4*pi.*Frequency/C.*(-30E-6));
%Spectrum_Devided=Spectrum_New./Spectrum_Used_Reference_Shifted;
%Spectrum_Devided(N_f+1:N_t)=0;
%Signal_Devided=fft(Spectrum_Devided);
%Signal_Devided(1:1460)=0;
%Signal_Devided(round(length(Signal_Devided)/2:end))=0;
%Spectrum_Devided_New=ifft(Signal_Devided);
%Spectrum_Devided_New=Spectrum_Devided_New(1:N_f);
%Spectrum_Devided_New=Spectrum_Devided_New.*exp(i*4*pi.*Frequency/C.*(30E-6));

%Spectrum_Devided=Spectrum_New./Spectrum_Used_Reference;
    
%% Start the n k fitting    

if index_calculation==1   
    
    
    e=e1+i*e2;
    e1=A0-A1*Frequency./((Frequency.^2)+A2^2);
    e2=A3./((Frequency.^2)+A2^2);
    %n=e1.^0.5;      %approxiamation
    %k=e2./2./n;
        
    n=((e1+(e1.^2+e2.^2)^0.5)^0.5)/(2^0.5);    
    k=((e1-(e1.^2+e2.^2)^0.5)^0.5)/(2^0.5);
    
    while (q<Number_of_Loop)
        e=e1+i*e2;
        e1=A0-A1*Frequency./((Frequency.^2)+A2^2);
        e2=A3./((Frequency.^2)+A2^2);
        %n=e1.^0.5;      %approxiamation
        %k=e2./2./n;  
        n=((e1+(e1.^2+e2.^2)^0.5)^0.5)/(2^0.5);    
        k=((e1-(e1.^2+e2.^2)^0.5)^0.5)/(2^0.5);
        
        r1=(n1-(n+i*k))./(n+n1+i*k);
        r1_r=((n+i*k)-n1(index_now))./(n+n1(index_now)+i*k);    
        t1=2*(n1(index_now))./(n+n1(index_now)+i*k);
        t1_r=2*(n+i*k)./(n+n1(index_now)+i*k);
        t2=2*(n+i*k)./(n+n2+i*k);
        t2_r=2*(n2)./(n+n2+i*k);
        r2=((n+i*k)-n2)./((n+i*k)+n2);   
        d=exp(i*2*pi.*Frequency/C.*(n+i*k).*Thickness);   %注意! -1*n!
        d0=exp(i*2*pi.*Frequency/C.*(Distance_fit(j)-Thickness));   %注意! -1*n!

        Spectrum_Upper_Temp=((r1)./r_BK7(index_now));                       %神說是cosine
        Spectrum_Lower_Temp=((t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now));                         %神說是cosine
        Spectrum_Temp=Spectrum_Upper_Temp+Spectrum_Lower_Temp;
        Error=(((Spectrum_Temp)-(Spectrum_Devided)).^2); %(abs(abs(Spectrum_Temp)-abs(Spectrum(index_now))).^2)+
        
        Error_Total=sum(abs(((Spectrum_Temp)-(Spectrum_Devided)).^2));
        
        Error_Difference=Error_Total-Error_Total_Previous;

        Error_Total_Previous=Error_Total;
    
        if Error_Difference > 0                    % Start to diverge
            
            if Current_Scanning_Parameter == 1              % Reset privious change
                A0=A0-Direction*delta_A0;
            elseif Current_Scanning_Parameter == 2
                A1=A1-Direction*delta_A1;
            elseif Current_Scanning_Parameter == 3
                A2=A2-Direction*delta_A2;
            elseif Current_Scanning_Parameter == 3
                A3=A3-Direction*delta_A3;
            elseif Current_Scanning_Parameter == 4
                A4=A4-Direction*delta_A4;
            end
                
            if Current_Scanning_Parameter < 4       %5: for last parameter Efficiency_3
                Current_Scanning_Parameter=Current_Scanning_Parameter+1;
            else
                Direction=Direction*(-1);
                Current_Scanning_Parameter=1;
            end
        end
        
        if Current_Scanning_Parameter == 1              % Reset privious change
           A0=A0-Direction*delta_A0;
        elseif Current_Scanning_Parameter == 2
           A1=A1-Direction*delta_A1;
        elseif Current_Scanning_Parameter == 3
           A2=A2-Direction*delta_A2;
        elseif Current_Scanning_Parameter == 4
           A3=A3-Direction*delta_A3;
        end  
    end
end    
end


%% Checking


%dlmwrite('Wavelength_micron.txt',Wavelength_micron,'delimiter','\t','newline','pc','precision','%.12f');
if index_calculation == 1
plot(Wavelength_micron(Wavelength_micron<0.8),n_final(Wavelength_micron<0.8,:));
plot(Wavelength_micron(Wavelength_micron<0.8),k_final(Wavelength_micron<0.8,:));

plot(Wavelength_micron(Wavelength_micron<0.8),Spectrum_Devided(Wavelength_micron<0.8,:),Wavelength_micron(Wavelength_micron<0.8),d_n_check(Wavelength_micron<0.8,:).^2);


dlmwrite('Wavelength_micron.txt',Wavelength_micron,'delimiter','\t','newline','pc');

dlmwrite('n_final.txt',n_final,'delimiter','\t','newline','pc');

dlmwrite('k_final.txt',k_final,'delimiter','\t','newline','pc');

dlmwrite('T_abs_check.txt',T_abs_check,'delimiter','\t','newline','pc');

end