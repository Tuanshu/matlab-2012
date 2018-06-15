%% Options

clear all

index_initial_boundary_reference=2750;
index_initial_boundary_sample=2950;

array=1:999;
lateral_index_reference_start=900;
lateral_index_reference_end=900;

lateral_index_sample_start=600;
lateral_index_sample_end=600;

Distance_Cor=0;

range_specified=1;

Reference_Align=0;

Ratio=1;
Ratio2=0.25;
Number_of_Loop=5;
Thickness=2.7.*1E-6;
Load_Previous_Result=0;
Position_Error=0.272*1E-6;%[0.019:0.001:0.024]*1E-6;%[0.02:0.0003:0.023]*1E-6;
n_should=1.7;
delta_n=0.1;
Wavelength_Center=540;

pixel_1=1250;               % DC cutoff
%pixel_2=1800;%2000               % the end of 2
%pixel_3=1120;               %for saperation

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


cd('D:\120220\G5\');
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

r_BK7=((1-n_bk7)./(n_bk7+1));
t_AN100=(2*(n_bk7)./(n_bk7+1));
%Spectrum_Original=Spectrum_Reference./(r_BK7).^2;

%% Reference plane calibration and tilting angle generation
Profile(array)=0;
Profile_Calibration(array)=0;
Profile_New(array)=0;
Distance(array)=0;
index_curernt_boundary=index_initial_boundary_reference;
for j=1:(lateral_index_reference_end-lateral_index_reference_start+1)
    cd('D:\120220\G5\');
    Data_Reference=importdata(sprintf('D%i.txt',j+lateral_index_reference_start-1));      % Data_2: the glass data 1
    Spectrum_Reference_Old=Data_Reference(:,2)-mean(Data_Reference(1:1000,2));
    Spectrum_Reference_Frequency=(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
    Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency); 
    Spectrum_Reference(isnan(Spectrum_Reference))=0;
    Spectrum_Reference(Frequency<Min_Frequency)=0;
    Spectrum_Reference((N_f+1):N_t)=0;

    Signal_Reference=fft(Spectrum_Reference);
    Signal_Reference(1:pixel_1)=0;
    Signal_Reference(round(length(Signal_Reference)/2):end)=0;
    
    [maxvalue_Reference maxindex_Reference]=max(abs(Signal_Reference(1:index_curernt_boundary)));
    Profile(j+lateral_index_reference_start-1)=Position(maxindex_Reference);

    Signal_Reference_Calibration=Signal_Reference;
    Signal_Reference_Calibration(1:index_curernt_boundary)=0;
    
    [maxvalue_Reference_Calibration maxindex_Reference_Calibration]=max(abs(Signal_Reference_Calibration));
    Profile_Calibration(j+lateral_index_reference_start-1)=Position(maxindex_Reference_Calibration);
    
    if j == round((lateral_index_reference_end-lateral_index_reference_start)/2)
        Signal_Reference_Calibration_Center=Signal_Reference_Calibration;
        Spectrum_Reference_Calibration_Center=(ifft(Signal_Reference_Calibration_Center));
        Spectrum_Reference_Calibration_Center=2*(Spectrum_Reference_Calibration_Center(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform
    end
    
    if abs((maxindex_Reference+maxindex_Reference_Calibration+index_curernt_boundary)/2-index_curernt_boundary) < 100
        index_curernt_boundary=round((maxindex_Reference+maxindex_Reference_Calibration+index_curernt_boundary)/2);
    end 

    Distance(j+lateral_index_reference_start-1)=Profile_Calibration(j+lateral_index_reference_start-1)-Profile(j+lateral_index_reference_start-1);
    Distance(j+lateral_index_reference_start-1)=Distance(j+lateral_index_reference_start-1)+Distance_Cor;
    %Spectrum_Additional=exp(i*4*pi.*Frequency/C.*(Profile_Calibration(lateral_index_reference_start)-Profile_Calibration(j+lateral_index_reference_start-1)));        %總之都移到Profile_Reference(1)
    
    %Spectrum_Reference_New=Spectrum_Reference.*Spectrum_Additional;
    %Spectrum_Reference_New((N_f+1):N_t)=0;
    
    %Signal_Reference_New=fft(Spectrum_Reference_New);
    %Signal_Reference_New(1:pixel_1)=0;
    %Signal_Reference_New(round(length(Signal_Reference_New)/2):end)=0;
    %Signal_Reference_New((pixel_2+1):end)=0;
    %[maxvalue_Reference_New maxindex_Reference_New]=max(abs(Signal_Reference_New));
    %Profile_New(j+lateral_index_reference_start-1)=Position(maxindex_Reference_New);
    
end

% For tilting

A=polyfit(array(lateral_index_reference_start:lateral_index_reference_end),Distance(lateral_index_reference_start:lateral_index_reference_end),1);
Distance_fit=A(1)*array+A(2);   %the distance between upper interface and additional interface!!
%Profile_New(lateral_index_reference_start:lateral_index_reference_end)=Profile_New(lateral_index_reference_start:lateral_index_reference_end)-Basedline(lateral_index_reference_start:lateral_index_reference_end)+Basedline(lateral_index_reference_start);

%% Sample plane related calculation (at each position)
index_curernt_boundary=index_initial_boundary_sample;
n_sum=0;
k_sum=0;
for j=1:(lateral_index_sample_end-lateral_index_sample_start+1)
    Data=importdata(sprintf('D%i.txt',j+lateral_index_sample_start-1));      % Data_2: the glass data 1
    Spectrum_Old=Data(:,2)-mean(Data(1:1000,2));
    Spectrum_Frequency=(Spectrum_Old.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
    Spectrum=interp1(Frequency_Old,Spectrum_Frequency,Frequency); 
    Spectrum(isnan(Spectrum))=0;
    Spectrum(Frequency<Min_Frequency)=0;
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
    
    [maxvalue_Calibration maxindex_Calibration]=max(abs(Signal_Calibration));
    Profile_Calibration(j+lateral_index_sample_start-1)=Position(maxindex_Calibration);   
    
    if abs((maxindex+maxindex_Calibration+index_curernt_boundary)/2-index_curernt_boundary) < 100
        index_curernt_boundary=round((maxindex+maxindex_Calibration+index_curernt_boundary)/2);
    end 
    
    Signal_Sample=Signal;
    Signal_Calibration=Signal;
    
    Signal_Sample((index_curernt_boundary+1):end)=0;
    Signal_Calibration(1:index_curernt_boundary)=0;
    
    
    Spectrum_Sample=(ifft(Signal_Sample));
    Spectrum_Sample=2*(Spectrum_Sample(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform
    
    
    Spectrum_Calibration=(ifft(Signal_Calibration));
    Spectrum_Calibration=2*(Spectrum_Calibration(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform
    
    %% This part is only for drawing, for the real calculation, only
    %% Spectrum_Additional_based should be used, since the
    %% Spectrum_Additional term will be null in ther afterward calculation 
    %Spectrum_Additional=exp(i*4*pi.*Frequency/C.*(Profile_Calibration(lateral_index_reference_start)-Profile_Calibration(j+lateral_index_sample_start-1)));    %總之都移到Profile_Reference(1)
    %Spectrum_Additional_Based=exp(i*4*pi.*Frequency/C.*(Basedline(lateral_index_reference_start)-Basedline(j+lateral_index_sample_start-1)));    
    %Spectrum_New=Spectrum.*Spectrum_Additional.*Spectrum_Additional_Based;
    
    %Spectrum_Additional=exp(i*4*pi.*Frequency/C.*(Profile_Calibration(lateral_index_reference_start)-Profile_Calibration(j+lateral_index_sample_start-1)));    %總之都移到Profile_Reference(1)
    %Spectrum_Additional_Based=exp(i*4*pi.*Frequency/C.*(-Basedline(j+lateral_index_sample_start-1)));    
    %Spectrum_Sample_New=Spectrum.*Spectrum_Additional_Based;
    %Spectrum_Calibration_New=Spectrum_Calibration.*Spectrum_Additional_Based;
    Spectrum_New=Spectrum_Sample./Spectrum_Calibration;
    
    %Signal_New=fft(Spectrum_New);
    %Signal_New(1:pixel_1)=0;
    %Signal_New(round(length(Signal_New)/2):end)=0;
    %[maxvalue_New maxindex_New]=max(abs(Signal_New));
    %Profile_New(j+lateral_index_sample_start-1)=Position(maxindex_New);
    
%% Start the n k fitting    
    
    n_final(1:length(Frequency),1:length(n_should))=1.5;
    k_final(1:length(Frequency),1:length(n_should))=0;
    T_check(1:length(Frequency),1:length(n_should))=0;  
    Spectrum_Check(1:length(Frequency),1:length(n_should))=0;
    Spectrum_Temp_check(1:length(Frequency))=0;
    
    for q=1:length(n_should)
    
        value_temp=10000000000;
        n(1:length(Frequency),1)=1.5;
        k(1:length(Frequency),1)=0;

        %T_max_check(1:length(Spectrum),1:length(n_should))=0;

        value_total_1(1:Number_of_Loop)=0;
        value_total_2(1:Number_of_Loop)=0;

        n_o=(n_should-delta_n):0.0001:(n_should(q)+delta_n);
        k_o=-0.3:0.0006:0.3;
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
    
            index_old=0;
            d0_Check(1:length(Spectrum))=0;
            dref_Check(1:length(Spectrum))=0;
    
            for w=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)

                if w <= ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)
                    index_now=Frequency_Center_Index+w-1;
                elseif w > ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)           %direction reverse!!
                    index_now=-w+((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)+Frequency_Center_Index;
                end

                % About the conditions

                r1=(n1(index_now)-(n(index_now)+i*k_temp))./(n(index_now)+n1(index_now)+i*k_temp);
                r1_r=((n(index_now)+i*k_temp)-n1(index_now))./(n(index_now)+n1(index_now)+i*k_temp);    
                t1=2*(n1(index_now))./(n(index_now)+n1(index_now)+i*k_temp);
                t1_r=2*(n(index_now)+i*k_temp)./(n(index_now)+n1(index_now)+i*k_temp);
                t2=2*(n(index_now)+i*k_temp)./(n(index_now)+n2+i*k_temp);
                t2_r=2*(n2)./(n(index_now)+n2+i*k_temp);
                r2=((n(index_now)+i*k_temp)-n2)./((n(index_now)+i*k_temp)+n2);   
                d=exp(i*2*pi.*Frequency(index_now)/C.*(n(index_now)+i*k_temp).*Thickness);   %注意! -1*n!
                d0=exp(i*2*pi.*Frequency(index_now)/C.*(Distance_fit(j)-Thickness));   %注意! -1*n!

                d_to_ref=exp(i*2*pi.*Frequency(index_now)/C.*(Distance_fit(j)-Distance_fit(round((lateral_index_reference_start+lateral_index_reference_end)/2))));   %注意! -1*n!
                
                %Spectrum_Upper_Temp=((r1)./(r_BK7(index_now).*t1.*t1_r.*t2.*t2_r.*(d.^2).*(d0.^2)));                       %神說是cosine
                %Spectrum_Lower_Temp=((r2.*t1.*t1_r.*(d.^2))./(r_BK7(index_now).*t1.*t1_r.*t2.*t2_r.*(d.^2).*(d0.^2)))*Ratio;                       %神說是cosine
                %Spectrum_Temp=Spectrum_Upper_Temp+Spectrum_Lower_Temp;
                
                T_abs=abs(t1.*t1_r.*t2.*t2_r.*(d.^2));
                
                Merit=((T_abs-abs(Spectrum_Calibration(index_now)/Spectrum_Reference_Calibration_Center(index_now))*Ratio2).^2);%+((Spectroscopy_Temp-(Spectroscopy(index_now))).^2); %(abs(abs(Spectrum_Temp)-abs(Spectrum(index_now))).^2)+


                [value index]=min(Merit);
                value_total_2(p)=value_total_2(p)+abs(value);
                k(index_now)=k_temp(index);

            end

            
            for w=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
            
                if w <= ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)
                    index_now=Frequency_Center_Index+w-1;
                elseif w > ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)           %direction reverse!!
                    index_now=-w+((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)+Frequency_Center_Index;
                end
            
                r1=(n1(index_now)-(n_temp+i*k(index_now)))./(n_temp+n1(index_now)+i*k(index_now));
                r1_r=((n_temp+i*k(index_now))-n1(index_now))./(n_temp+n1(index_now)+i*k(index_now));    
                t1=2*(n1(index_now))./(n_temp+n1(index_now)+i*k(index_now));
                t1_r=2*(n_temp+i*k(index_now))./(n_temp+n1(index_now)+i*k(index_now));
                t2=2*(n_temp+i*k(index_now))./(n_temp+n2+i*k(index_now));
                t2_r=2*(n2)./(n_temp+n2+i*k(index_now));
                r2=((n_temp+i*k(index_now))-n2)./((n_temp+i*k(index_now))+n2);   
                d=exp(i*2*pi.*Frequency(index_now)/C.*(n_temp+i*k(index_now)).*Thickness);   %注意! -1*n!
                d0=exp(i*2*pi.*Frequency(index_now)/C.*(Distance_fit(j)-Thickness));   %注意! -1*n!
                Spectrum_Upper_Temp=((r1)./(r_BK7(index_now).*t1.*t1_r.*t2.*t2_r.*(d.^2).*(d0.^2)));                       %神說是cosine
                Spectrum_Lower_Temp=((r2)./(r_BK7(index_now).*t2.*t2_r.*(d0.^2)));                       %神說是cosine
                Spectrum_Temp=Spectrum_Upper_Temp+Spectrum_Lower_Temp;
                Merit=(angle(Spectrum_Temp)-angle(Spectrum_New(index_now)*Ratio)).^2; %(abs(abs(Spectrum_Temp)-abs(Spectrum(index_now))).^2)+
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
                Spectrum_Upper_Temp_check(index_now)=((r1(index))./(r_BK7(index_now).*t1(index).*t1_r(index).*t2(index).*t2_r(index).*(d(index).^2).*(d0.^2)));                       %神說是cosine
                Spectrum_Lower_Temp_check(index_now)=((r2(index).*t1(index).*t1_r(index).*(d(index).^2))./(r_BK7(index_now).*t1(index).*t1_r(index).*t2(index).*t2_r(index).*(d(index).^2).*(d0.^2)));     
                Spectrum_Temp_check(index_now)=Spectrum_Upper_Temp_check(index_now)+Spectrum_Lower_Temp_check(index_now);
            end

        end

        n_final(:,q)=n;
        k_final(:,q)=k;

    end
    n_sum=n_sum+n_final;
    k_sum=k_sum+k_final;
end

n_ave=n_sum/(lateral_index_sample_end-lateral_index_sample_start+1);
k_ave=k_sum/(lateral_index_sample_end-lateral_index_sample_start+1);
%% Fitting


%% Checking

r1_check=(n1-(n_final+i*k_final))./(n_final+n1+i*k_final);
r1_r_check=((n_final+i*k_final)-n1)./(n_final+n1+i*k_final);    
t1_check=2*(n1)./(n_final+n1+i*k_final);
t1_r_check=2*(n_final+i*k_final)./(n_final+n1(index_now)+i*k_final);
t2_check=2*(n_final+i*k_final)./(n_final+n2+i*k_final);
t2_r_check=2*(n2)./(n_final+n2+i*k_final);
r2_check=((n_final+i*k_final)-n2)./((n_final+i*k_final)+n2);   
d_check=exp(i*2*pi.*Frequency/C.*(n_final+i*k_final).*Thickness);   %注意! -1*n!
d0_check=exp(i*2*pi.*Frequency/C.*(Distance_fit(j)-Thickness));   %注意! -1*n!

Spectrum_Upper_check=((r1_check)./(r_BK7.*t1_check.*t1_r_check.*t2_check.*t2_r_check.*(d_check.^2).*(d0_check.^2)));                       %神說是cosine
Spectrum_Lower_check=((r2_check.*t1_check.*t1_r_check.*(d_check.^2))./(r_BK7.*t1_check.*t1_r_check.*t2_check.*t2_r_check.*(d_check.^2).*(d0_check.^2)));                       %神說是cosine
Spectrum_check=(Spectrum_Upper_check+Spectrum_Upper_check);

T_abs_check=abs(t1_check.*t1_r_check.*t2_check.*t2_r_check.*(d_check.^2));
                
%dlmwrite('Wavelength_micron.txt',Wavelength_micron,'delimiter','\t','newline','pc','precision','%.12f');
plot(Wavelength_micron(Wavelength_micron<0.6),n_ave(Wavelength_micron<0.6,:));
plot(Wavelength_micron(Wavelength_micron<0.6),k_ave(Wavelength_micron<0.6,:));