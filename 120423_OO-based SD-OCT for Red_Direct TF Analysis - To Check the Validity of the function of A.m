%% Options

clear all
index_calculation=1;
TF_Analysis=1;  %1: direct in FD, 2: transfer to TD first

array=1:999;
index_used_reference_plane=100;

lateral_index_calibration_start=50;
lateral_index_calibration_end=700;

lateral_index_reference_start=50;
lateral_index_reference_end=250;

lateral_index_sample_start=650;
lateral_index_sample_end=650;


range_specified=1;

Reference_Align=0;

Ratio=1;
Ratio2=1;
Number_of_Loop=1;
Thickness=2.5.*1E-6;
Load_Previous_Result=0;
Position_Upper_Interface=10.82*1E-6;%[0.019:0.001:0.024]*1E-6;%[0.02:0.0003:0.023]*1E-6;
%%
n_should=1.8;
delta_n=0.4;
Wavelength_Center=600;

Scanning='Thickness';

if strcmp(Scanning,'n_should')
    Parameter=n_should;
elseif strcmp(Scanning,'Thickness')
    Parameter=Thickness;
elseif strcmp(Scanning,'Position_Upper_Interface')
    Parameter=Position_Upper_Interface;
    
elseif strcmp(Scanning,'Ratio')
    Parameter=Ratio;
end

pixel_1=700;               % DC cutoff
%pixel_2=1800;%2000               % the end of 2
%pixel_3=1120;               %for saperation

Wavelength_Considered_Min=590;          %nm
Wavelength_Considered_Max=600;

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


cd('D:\120418\Filter_inter\');
Data=importdata('D0.txt');      % Data_2: the glass data 1

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

r_BK7=((n_bk7-1)./(n_bk7+1));
t_AN100=(2*(n_bk7)./(n_bk7+1));
%Spectrum_Original=Spectrum_Reference./(r_BK7).^2;

%% Data Loading

Filter_inter=0;
cd('D:\120418\filter_inter\');
for j=0:39
    Filter_inter=Filter_inter+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Filter_inter=Filter_inter/40;

Filter_sam=0;
cd('D:\120418\filter_sam\');
for j=0:39
    Filter_sam=Filter_sam+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Filter_sam=Filter_sam/40;


ref=0;
cd('D:\120418\ref\');
for j=0:39
    ref=ref+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
ref=ref/40;


Glass_inter=0;
cd('D:\120418\Glass_inter\');
for j=0:39
    Glass_inter=Glass_inter+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Glass_inter=Glass_inter/40;


Glass_sam=0;
cd('D:\120418\Glass_sam\');
for j=0:39
    Glass_sam=Glass_sam+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Glass_sam=Glass_sam/40;




bs=0;
cd('D:\120418\bs\');
for j=0:39
    bs=bs+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
bs=bs/40;


Spectrum_Sample_Wavelength=Filter_inter(:,2)-Filter_sam(:,2)-ref(:,2)+bs(:,2);
%Spectrum_Sample_Wavelength=ref(:,2);
Spectrum_Reference_Wavelength=Glass_inter(:,2)-Glass_sam(:,2)-ref(:,2)+bs(:,2);

Spectrum_Sample_Wavelength=Spectrum_Sample_Wavelength-Spectrum_Sample_Wavelength(1);

Spectrum_Reference_Wavelength=Spectrum_Reference_Wavelength-Spectrum_Reference_Wavelength(1);

plot(Wavelength,Spectrum_Sample_Wavelength);

xlabel('Wavelength (micron)');
ylabel('Spectral Power (a.u.)');

%% Spectrum generation

Gauss_Window=1;

Spectrum_Sample_Frequency=(Spectrum_Sample_Wavelength.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectrum_Sample=interp1(Frequency_Old,Spectrum_Sample_Frequency,Frequency,'spline'); 
Spectrum_Sample(isnan(Spectrum_Sample))=0;
Spectrum_Sample(Frequency<Min_Frequency)=0;
Spectrum_Sample(N_f+1:N_t)=0;
Signal_Sample=fft(Spectrum_Sample).*(1-gaussmf(Position_micron,[Gauss_Window 0]));  
Signal_Sample(round(length(Signal_Sample)/2)+1:end)=0;
Spectrum_Sample=(ifft(Signal_Sample));
Spectrum_Sample=2*Spectrum_Sample(1:N_f);

Spectrum_Reference_Frequency=(Spectrum_Reference_Wavelength.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency,'spline'); 
Spectrum_Reference(isnan(Spectrum_Reference))=0;
Spectrum_Reference(Frequency<Min_Frequency)=0;
Spectrum_Reference(N_f+1:N_t)=0;
Signal_Reference=fft(Spectrum_Reference).*(1-gaussmf(Position_micron,[Gauss_Window 0]));  

Signal_Reference(round(length(Signal_Reference)/2)+1:end)=0;
Spectrum_Reference=(ifft(Signal_Reference));
Spectrum_Reference=2*Spectrum_Reference(1:N_f);

Phase_Sample=angle(Spectrum_Sample);
Phase_Reference=angle(Spectrum_Reference);


for w=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
    index_now=w+Frequency_Considered_Min_Index;
    if (Phase_Sample(index_now)-Phase_Sample(index_now-1))>2*pi*0.9
        Phase_Sample(index_now:end)=Phase_Sample(index_now:end)-2*pi;
    elseif (Phase_Sample(index_now)-Phase_Sample(index_now-1))<-2*pi*0.9
        Phase_Sample(index_now:end)=Phase_Sample(index_now:end)+2*pi;
    end
    if (Phase_Reference(index_now)-Phase_Reference(index_now-1))>2*pi*0.9
        Phase_Reference(index_now:end)=Phase_Reference(index_now:end)-2*pi;
    elseif (Phase_Reference(index_now)-Phase_Reference(index_now-1))<-2*pi*0.9
        Phase_Reference(index_now:end)=Phase_Reference(index_now:end)+2*pi;
    end
end

                
Spectrum_Devided=Spectrum_Sample./Spectrum_Reference;

plot(Position_micron,Signal_Sample,Position_micron,Signal_Reference);

xlabel('Optical path Difference (micron)');
ylabel('Interference Signal');

plot(Wavelength_micron(Wavelength_micron<0.80),real(Spectrum_Reference(Wavelength_micron<0.80,:)),Wavelength_micron(Wavelength_micron<0.80),real(Spectrum_Sample(Wavelength_micron<0.80,:)));


plot(Wavelength_micron(Wavelength_micron<0.80),real(Spectrum_Devided(Wavelength_micron<0.80,:)),Wavelength_micron(Wavelength_micron<0.80),abs(Spectrum_Devided(Wavelength_micron<0.80,:)));

xlabel('Wavelength (micron)');
ylabel('A');

%% to generate the phase spectrum

%% Start the n k fitting    

if index_calculation==1   
    
    n_final(1:length(Frequency),1:length(Parameter))=1.5;
    k_final(1:length(Frequency),1:length(Parameter))=0;
    Spectrum_Upper_Record(1:length(Frequency),1:length(Parameter))=0;
    Spectrum_Lower_Record(1:length(Frequency),1:length(Parameter))=0;
    T_check(1:length(Frequency),1:length(Parameter))=0;  
    Spectrum_Check(1:length(Frequency),1:length(Parameter))=0;
    for q=1:length(Parameter)
        
    if strcmp(Scanning,'n_should')
        n_should_Now=n_should(q);
        Thickness_Now=Thickness;
        Ratio_Now=Ratio;
        Position_Error_Now=Position_Upper_Interface;
    elseif strcmp(Scanning,'Thickness')
        n_should_Now=n_should;
        Thickness_Now=Thickness(q);
        Ratio_Now=Ratio;
        Position_Error_Now=Position_Upper_Interface;
    elseif strcmp(Scanning,'Position_Upper_Interface')
        n_should_Now=n_should;
        Thickness_Now=Thickness;
        Ratio_Now=Ratio;
        Position_Error_Now=Position_Upper_Interface(q);
    elseif strcmp(Scanning,'Ratio')
        n_should_Now=n_should;
        Thickness_Now=Thickness;
        Ratio_Now=Ratio(q);
        Position_Error_Now=Position_Upper_Interface;
    end
        value_temp=10000000000;
        n(1:length(Frequency),1)=1.5;
        k(1:length(Frequency),1)=0;
        Spectrum_Upper(1:length(Frequency),1)=0;
        Spectrum_Lower(1:length(Frequency),1)=0;
        %T_max_check(1:length(Spectrum),1:length(n_should))=0;

        value_total_1(1:Number_of_Loop)=0;
        value_total_2(1:Number_of_Loop)=0;

        n_o=(n_should_Now-delta_n):0.001:(n_should_Now+delta_n);
        k_o=-0.1:0.001:0.3;
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
    
        n_should_index=find(n_o-n_should_Now>0,1,'first');
    
        
        A=(1./r_BK7) .* ((n1-(n+i*k))./(n1+(n+i*k))+(exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now) ./ (1-(((n+i*k)-n1)./((n+i*k)+n1)) .* (((n+i*k)-n2)./((n+i*k)+n2)) .* exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now))) .* ((2*n1)./(n1+(n+i*k))) .* ((2*(n+i*k))./(n1+(n+i*k))) .* (((n+i*k)-n2)./((n+i*k)+n2)));
        dA_dn=(1./r_BK7) .* (-2*n1./(n1+(n+i*k)).^2 + (exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now) ./ (1-(((n+i*k)-n1)./((n+i*k)+n1)) .* (((n+i*k)-n2)./((n+i*k)+n2)) .* exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now))) .* ((2*n1)./(n1+(n+i*k))) .* ((2*(n+i*k))./(n1+(n+i*k))) .* (((n+i*k)-n2)./((n+i*k)+n2)) .* (-1./(n1+(n+i*k)) + n1./ (n1+(n+i*k)) ./ (n+i*k) + 2*n2 ./ ((n+i*k)+n2) ./ ((n+i*k)-n2) +i*4*pi*Thickness_Now.*Frequency/C + (exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now) ./ (1-(((n+i*k)-n1)./((n+i*k)+n1)) .* (((n+i*k)-n2)./((n+i*k)+n2)) .* exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now))) .* (((n+i*k)-n1)./((n+i*k)+n1)) .* (((n+i*k)-n2)./((n+i*k)+n2)) .*(2*n1 ./ ((n+i*k)+n1) ./ ((n+i*k)-n1) + 2*n2 ./ ((n+i*k)+n2) ./ ((n+i*k)-n2) + i*4*pi*Thickness_Now.*Frequency/C)));
        dA_dn_Reduced=(1./r_BK7) .* (-2*n1./(n1+(n+i*k)).^2 + (exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now) .* ((2*n1)./(n1+(n+i*k))) .* ((2*(n+i*k))./(n1+(n+i*k))) .* (((n+i*k)-n2)./((n+i*k)+n2)) .* (-1./(n1+(n+i*k)) + n1./ (n1+(n+i*k)) ./ (n+i*k) + 2*n2 ./ ((n+i*k)+n2) ./ ((n+i*k)-n2) +i*4*pi*Thickness_Now.*Frequency/C)));
        dA_dn_Lower_Reduced=(1./r_BK7) .* (exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now) .* ((2*n1)./(n1+(n+i*k))) .* ((2*(n+i*k))./(n1+(n+i*k))) .* (((n+i*k)-n2)./((n+i*k)+n2)) .* (-1./(n1+(n+i*k)) + n1./ (n1+(n+i*k)) ./ (n+i*k) + 2*n2 ./ ((n+i*k)+n2) ./ ((n+i*k)-n2) +i*4*pi*Thickness_Now.*Frequency/C));
        %dA_dn_Lower_Patial=exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now).*((2*n1)./(n1+(n+i*k))) .* ((2*(n+i*k))./(n1+(n+i*k))) .* (((n+i*k)-n2)./((n+i*k)+n2)) .* (-1./(n1+(n+i*k)) + n1./ (n1+(n+i*k)) ./ (n+i*k) + 2*n2 ./ ((n+i*k)+n2) ./ ((n+i*k)-n2)+i*4*pi*Thickness_Now.*Frequency/C);
        dA_dn_Lower_Patial=exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now).*((2*n1)./(n1+(n+i*k))) .* ((2*(n+i*k))./(n1+(n+i*k))) .* (((n+i*k)-n2)./((n+i*k)+n2)) .* (-1./(n1+(n+i*k)) + n1./ (n1+(n+i*k)) ./ (n+i*k) + 2*n2 ./ ((n+i*k)+n2) ./ ((n+i*k)-n2)+i*4*pi*Thickness_Now.*Frequency/C);
        dA_dn_Lower_Patial_Exp=exp(i*4*pi.*Frequency/C.*(n+i*k)*Thickness_Now).*i*4*pi*Thickness_Now.*Frequency/C;
        dA_dn_Upper=(1./r_BK7) .* (-2*n1./(n1+(n+i*k)).^2);
        
        for p=1:Number_of_Loop
    
            index_old=0;
            d0_Check(1:length(Frequency))=0;
            dref_Check(1:length(Frequency))=0;
    
            for w=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
            
                if w <= ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)
                    index_now=Frequency_Center_Index+w-1;
                elseif w > ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)           %direction reverse!!
                    index_now=-w+((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)+Frequency_Center_Index;
                end
            
                r1=(n1(index_now)-(n(index_now)+i*k(index_now)))./(n(index_now)+n1(index_now)+i*k(index_now));
                r1_r=((n(index_now)+i*k(index_now))-n1(index_now))./(n(index_now)+n1(index_now)+i*k(index_now));    
                t1=2*(n1(index_now))./(n(index_now)+n1(index_now)+i*k(index_now));
                t1_r=2*(n(index_now)+i*k(index_now))./(n(index_now)+n1(index_now)+i*k(index_now));
                t2=2*(n(index_now)+i*k(index_now))./(n(index_now)+n2+i*k(index_now));
                t2_r=2*(n2)./(n(index_now)+n2+i*k(index_now));
                r2=((n(index_now)+i*k(index_now))-n2)./((n(index_now)+i*k(index_now))+n2);   
                d=exp(i*2*pi.*Frequency(index_now)/C.*(n(index_now)+i*k(index_now)).*Thickness_Now);   %注意! -1*n!
                d0=exp(i*2*pi.*Frequency(index_now)/C.*(Position_Error_Now));   %注意! -1*n!
                Spectrum_Upper_Temp=((r1)./r_BK7(index_now));         
                Spectrum_Lower_Temp=((t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now))./(1-r1_r.*(d.^2).*r2);                         %神說是cosine     
                Spectrum_Lower_Temp_Reduced_1(index_now)=(t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now);                         %神說是cosine     
                Spectrum_Lower_Temp_Patial_1(index_now)=(t1.*t1_r.*r2.*(d.^2));                   
                Spectrum_Lower_Temp_Exp_1(index_now)=((d.^2));   
                %Spectrum_Lower_Temp_Patial_1(index_now)=(t1_r);   
                Spectrum_Temp_1(index_now)=Spectrum_Upper_Temp+Spectrum_Lower_Temp;%+Spectrum_Multi_Temp;
                Spectrum_Temp_Reduced_1(index_now)=Spectrum_Upper_Temp+Spectrum_Lower_Temp_Reduced_1(index_now);%+Spectrum_Multi_Temp;
                Spectrum_Upper_Temp_1(index_now)=Spectrum_Upper_Temp;
                QQ(w)=A(index_now)-Spectrum_Temp_1(index_now);

            end
            
            delta_n=0.000000000001;
            delta_k=0;
            n=n+delta_n;
            k=k+delta_k;
            
            for w=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
            
                if w <= ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)
                    index_now=Frequency_Center_Index+w-1;
                elseif w > ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)           %direction reverse!!
                    index_now=-w+((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)+Frequency_Center_Index;
                end
            
                r1=(n1(index_now)-(n(index_now)+i*k(index_now)))./(n(index_now)+n1(index_now)+i*k(index_now));
                r1_r=((n(index_now)+i*k(index_now))-n1(index_now))./(n(index_now)+n1(index_now)+i*k(index_now));    
                t1=2*(n1(index_now))./(n(index_now)+n1(index_now)+i*k(index_now));
                t1_r=2*(n(index_now)+i*k(index_now))./(n(index_now)+n1(index_now)+i*k(index_now));
                t2=2*(n(index_now)+i*k(index_now))./(n(index_now)+n2+i*k(index_now));
                t2_r=2*(n2)./(n(index_now)+n2+i*k(index_now));
                r2=((n(index_now)+i*k(index_now))-n2)./((n(index_now)+i*k(index_now))+n2);   
                d=exp(i*2*pi.*Frequency(index_now)/C.*(n(index_now)+i*k(index_now)).*Thickness_Now);   %注意! -1*n!
                d0=exp(i*2*pi.*Frequency(index_now)/C.*(Position_Error_Now));   %注意! -1*n!
                Spectrum_Upper_Temp=((r1)./r_BK7(index_now));         
                Spectrum_Lower_Temp=((t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now))./(1-r1_r.*(d.^2).*r2);                         %神說是cosine        
                Spectrum_Lower_Temp_Reduced_2=(t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now);                         %神說是cosine     
                Spectrum_Lower_Temp_Patial_2=(t1.*t1_r.*r2.*(d.^2));                         %神說是cosine     
                Spectrum_Lower_Temp_Exp_2=(d.^2);        
                %Spectrum_Lower_Temp_Patial_2=(t1_r);                         %神說是cosine     
                Spectrum_Temp_2=Spectrum_Upper_Temp+Spectrum_Lower_Temp;%+Spectrum_Multi_Temp;
                Spectrum_Temp_Reduced_2=Spectrum_Upper_Temp+Spectrum_Lower_Temp_Reduced_2;%+Spectrum_Multi_Temp;
                
                QQ2(w)=dA_dn(index_now)-(Spectrum_Temp_2-Spectrum_Temp_1(index_now))/delta_n;
                QQ3(w)=dA_dn_Reduced(index_now)-(Spectrum_Temp_Reduced_2-Spectrum_Temp_Reduced_1(index_now))/delta_n;
                QQ4(w)=dA_dn_Lower_Reduced(index_now)-(Spectrum_Lower_Temp_Reduced_2-Spectrum_Lower_Temp_Reduced_1(index_now))/delta_n;
                QQ5(w)=dA_dn_Lower_Patial(index_now)-(Spectrum_Lower_Temp_Patial_2-Spectrum_Lower_Temp_Patial_1(index_now))/delta_n;
                
                QQ6(w)=dA_dn_Lower_Patial_Exp(index_now)-(Spectrum_Lower_Temp_Exp_2-Spectrum_Lower_Temp_Exp_1(index_now))/delta_n;
            end

        end
    end

    
end    

%plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered
%_Max_Index),QQ,Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),QQ2,Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),QQ3,Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),QQ4);
plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),QQ2);