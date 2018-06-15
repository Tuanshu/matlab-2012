%% Options

clear all
index_calculation=1;
TF_Analysis=1;  %1: direct in FD, 2: transfer to TD first

Averaging=1;    %Number of frame averaging
Sample_Path='D:\120820_YPR_Face Up_2\YPR\';
Reference_Path='D:\120820_YPR_Face Up_2\glass\';
Spectroscopy_Path='D:\120808_YPR\';

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

Ratio_Upper2Reference=1+[-0.2:0.05:0.2];
Position_Upper2Reference=(3+[-0.5:0.05:0.5]*1E-6);
Number_of_Loop=100;
Thickness=(4.4+[-0.2:0.04:0.2]).*1E-6;% ([-0.1:0.1:0.1]+2.55).*1E-6;
Load_Previous_Result=0;
Ratio_Lower2Upper=1+[-0.2:0.04:0.2];  %([-1:0.02:1]+10.82)*1E-6;%[0.019:0.001:0.024]*1E-6;%[0.02:0.0003:0.023]*1E-6;
%%
n_should=1.6;
delta_n=0.4;
Wavelength_Center=550;

pixel_1=700;               % DC cutoff
%pixel_2=1800;%2000               % the end of 2
%pixel_3=1120;               %for saperation

Center_Wavelength_micron=0.58;
Wavelength_Considered_Min=500;          %nm
Wavelength_Considered_Max=650;

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


cd(sprintf('%sinter\\',Sample_Path));
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

%% T 


cd(sprintf('%s\\',Spectroscopy_Path));
Data_Spectroscopy=importdata('1208030331_5micron 1.jws.txt');

Spectrum_Spectroscopy_Old=Data_Spectroscopy(end:-1:1,2)/100;
Wavelength_Spectroscopy=Data_Spectroscopy(end:-1:1,1);     
Frequency_Spectroscopy=C./(Wavelength_Spectroscopy*1E-9);

T=interp1(Frequency_Spectroscopy,Spectrum_Spectroscopy_Old,Frequency,'spline'); 

%% Theory - Sample Model (n1 - n - n2)

n1=1;

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
n2=n_bk7;

%% To generate the original spectrum, assuming reference is BK7, too

r_BK7=((1-n_bk7)./(n_bk7+1));
t_AN100=(2*(n_bk7)./(n_bk7+1));
%Spectrum_Original=Spectrum_Reference./(r_BK7).^2;

%% Data Loading

Filter_inter=0;
cd(sprintf('%sinter\\',Sample_Path));
for j=0:(Averaging-1)
    Filter_inter=Filter_inter+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Filter_inter=Filter_inter/Averaging;

Filter_sam=0;
cd(sprintf('%ssam\\',Sample_Path));
for j=0:(Averaging-1)
    Filter_sam=Filter_sam+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Filter_sam=Filter_sam/Averaging;


Filter_ref=0;
cd(sprintf('%sref\\',Sample_Path));
for j=0:(Averaging-1)
    Filter_ref=Filter_ref+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Filter_ref=Filter_ref/Averaging;


Filter_bs=0;
cd(sprintf('%sbs\\',Sample_Path));
for j=0:(Averaging-1)
    Filter_bs=Filter_bs+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Filter_bs=Filter_bs/Averaging;


Glass_inter=0;
cd(sprintf('%sinter\\',Reference_Path));
for j=0:(Averaging-1)
    Glass_inter=Glass_inter+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Glass_inter=Glass_inter/Averaging;

Glass_sam=0;
cd(sprintf('%ssam\\',Reference_Path));
for j=0:(Averaging-1)
    Glass_sam=Glass_sam+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Glass_sam=Glass_sam/Averaging;


Glass_ref=0;
cd(sprintf('%sref\\',Reference_Path));
for j=0:(Averaging-1)
    Glass_ref=Glass_ref+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Glass_ref=Glass_ref/Averaging;


Glass_bs=0;
cd(sprintf('%sbs\\',Reference_Path));
for j=0:(Averaging-1)
    Glass_bs=Glass_bs+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Glass_bs=Glass_bs/Averaging;

Spectrum_Sample_Wavelength=Filter_inter(:,2);%-Filter_sam(:,2)-Filter_ref(:,2)+Filter_bs(:,2);
Spectrum_Reference_Wavelength=Glass_ref(:,2);%Glass_inter(:,2);%-Glass_sam(:,2)-Glass_ref(:,2)+Glass_bs(:,2);
%Spectrum_Sample_Wavelength=ref(:,2);

Spectrum_Sample_Wavelength=Spectrum_Sample_Wavelength-Spectrum_Sample_Wavelength(1);
Spectrum_Reference_Wavelength=Spectrum_Reference_Wavelength-Spectrum_Reference_Wavelength(1);

%plot(Wavelength,Spectrum_Sample_Wavelength);

%xlabel('Wavelength (micron)');
%ylabel('Spectral Power (a.u.)');

%% Spectrum generation

Separation_Position=30.2;   %(micron)
Gauss_Window_sig=1;
Gauss_Window_C=2;

Gauss_Window_Profile=1-gaussmf(Position_micron,[Gauss_Window_sig Gauss_Window_C]);
Gauss_Window_Profile(Position_micron<Gauss_Window_C)=0;

Spectrum_Sample_Frequency=(Spectrum_Sample_Wavelength.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectrum_Sample=interp1(Frequency_Old,Spectrum_Sample_Frequency,Frequency,'spline'); 
Spectrum_Sample(isnan(Spectrum_Sample))=0;
Spectrum_Sample(Frequency<Min_Frequency)=0;
Spectrum_Sample(N_f+1:N_t)=0;
Signal_Sample=fft(Spectrum_Sample);%.*Gauss_Window_Profile;  
Signal_Sample(round(length(Signal_Sample)/2)+1:end)=0;


Spectrum_Reference_Frequency=(Spectrum_Reference_Wavelength.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency,'spline'); 
Spectrum_Reference(isnan(Spectrum_Reference))=0;
Spectrum_Reference(Frequency<Min_Frequency)=0;
Spectrum_Reference(N_f+1:N_t)=0;
Signal_Reference=fft(Spectrum_Reference);%.*Gauss_Window_Profile;  
Signal_Reference(round(length(Signal_Reference)/2)+1:end)=0;


Signal_Sample_1=Signal_Sample;
Signal_Sample_2=Signal_Sample;

Signal_Sample_1(Position_micron>Separation_Position)=0;
Signal_Sample_2(Position_micron<=Separation_Position)=0;

Spectrum_Sample=(ifft(Signal_Sample));
Spectrum_Sample=2*Spectrum_Sample(1:N_f);


Spectrum_Reference=(ifft(Signal_Reference));
Spectrum_Reference=2*Spectrum_Reference(1:N_f);


Spectrum_Sample_1=(ifft(Signal_Sample_1));
Spectrum_Sample_1=2*Spectrum_Sample_1(1:N_f);

Spectrum_Sample_2=(ifft(Signal_Sample_2));
Spectrum_Sample_2=2*Spectrum_Sample_2(1:N_f);

Spectrum_Devided=Spectrum_Sample./Spectrum_Reference;
Spectrum_Devided_2=Spectrum_Sample_2./Spectrum_Sample_1;

Phase_Sample=angle(Spectrum_Sample);

plot(Position_micron,Signal_Sample_1,Position_micron,Signal_Sample_2);
plot(Position_micron,Signal_Sample,Position_micron,Signal_Reference);


plot(Wavelength_micron(Wavelength_micron<1),Spectrum_Sample_1(Wavelength_micron<1),Wavelength_micron(Wavelength_micron<1),Spectrum_Sample_2(Wavelength_micron<1));

for w=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
    index_now=w+Frequency_Considered_Min_Index;
    if (Phase_Sample(index_now)-Phase_Sample(index_now-1))>2*pi*0.9
        Phase_Sample(index_now:end)=Phase_Sample(index_now:end)-2*pi;
    elseif (Phase_Sample(index_now)-Phase_Sample(index_now-1))<-2*pi*0.9
        Phase_Sample(index_now:end)=Phase_Sample(index_now:end)+2*pi;
    end
end

%% Start the n k fitting    

if index_calculation==1   
    Merit_Best=999999999999999999999;
    for p=1:length(Thickness)
        for q=1:length(Ratio_Lower2Upper)
            for w=1:length(Ratio_Upper2Reference)
                for m=1:length(Position_Upper2Reference)
                    

                    Thickness_Now=Thickness(p);
                    Ratio_Lower2Upper_Now=Ratio_Lower2Upper(q);
                    Ratio_Upper2Reference_Now=Ratio_Upper2Reference(w);
                    Position_Upper2Reference_Now=Position_Upper2Reference(m);
                            
                    Frequency_Considered=Frequency(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
                    Aexp=Spectrum_Devided(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);   %注意! -1*n!;
                    T_Considered=T(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
                    %Weight_Function=(abs(Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index)).*T_Considered)./max((abs(Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index).*T_Considered)));
                    Wavelength_micron_Considered=Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
                    n1_Considered=n1(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
                    n2_Considered=n2;
                    r_BK7_Considered=r_BK7(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
                    t_AN100_Considered=t_AN100(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
                    n(1:length(Frequency_Considered),1)=n_should;%1.4+(Wavelength_micron_Considered-0.6)*2.8-(Wavelength_micron_Considered-0.7)*1.4;
                    %A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now) (exp(i*4*pi.*Frequency_Considered/C.*(Ratio_Lower2Upper_Now))./r_BK7_Considered./Ratio_Upper2Reference_Now) .* (Position_Upper2Reference_Now.*(n1_Considered-n)./(n1_Considered+n)+(exp(i*4*pi.*Frequency_Considered/C.*n*Thickness_Now) ./ (1-((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/C.*n*Thickness_Now))) .* ((2*n1_Considered)./(n1_Considered+n)) .* ((2*n)./(n1_Considered+n)) .* ((n-n2_Considered)./(n+n2_Considered)));
                    t1= @ (n) (2*n1_Considered)./(n1_Considered+n); 
                    t1_r= @ (n) (2*n)./(n2_Considered+n);
                    t2= @ (n) (2*n)./(n+n2_Considered); 
                    r1= @ (n) (n1_Considered-n)./(n1_Considered+n);
                    r1_r= @ (n) (n-n1_Considered)./(n+n1_Considered);
                    r2= @ (n) (n-n2_Considered)./(n+n2_Considered);
                    A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now) Ratio_Lower2Upper_Now./r1(n).*(t1(n).*t1_r(n).*r2(n).*exp(i*4*pi.*Frequency_Considered.*n.*Thickness_Now/C)./(1-r1_r(n).*r2(n).*exp(i*4*pi*Frequency_Considered.*n.*Thickness_Now/C)));
                    B = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference) exp(i*4*pi.*Frequency_Considered.*Position_Upper2Reference/C).*Position_Upper2Reference./r_BK7_Considered.*(r_1+Ratio_Lower2Upper_Now.*(t1.*t1_r.*r2.*exp(i*2.*Frequency_Considered.*n.*Thickness_Now/C)./(1-r1_r.*r2.*exp(i*2*Frequency_Considered.*n.*Thickness_Now/C))));
                    
                    %T_abs= @ (n,Thickness_Now) abs((2.*n2_Considered./(n2_Considered+n1_Considered)).*(2*n1_Considered./(n+n1_Considered)).*(2*n./(n+n2_Considered)).*exp(i*2*pi.*Frequency_Considered/C.*n.*Thickness_Now)./(1-((n-n1_Considered)./(n+n1_Considered)).*((n-n2_Considered)./(n+n2_Considered)).*exp(i*4*pi.*Frequency_Considered/C.*n.*Thickness_Now))).^2;
                    T_abs = @ (n,Thickness_Now) abs(t_AN100_Considered.*t1(n).*t2(n).*exp(i*4*pi*Frequency_Considered.*n*Thickness_Now/C)./(1-r1_r(n).*r2(n).*exp(i*4*pi*Frequency_Considered.*n*Thickness_Now/C))).^2;


                    %dA_dn = @ (n) (exp(i*4*pi.*Frequency_Considered/C.*(Ratio_Lower2Upper))./r_BK7_Considered./Ratio_Upper2Reference_Now) .* (Position_Upper2Reference_Now.*(-2*n1_Considered./(n1_Considered+n).^2) + (exp(i*4*pi.*Frequency_Considered/C.*n*Thickness_Now) ./ (1-((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/C.*n*Thickness_Now))) .* ((2*n1_Considered)./(n1_Considered+n)) .* ((2*n)./(n1_Considered+n)) .* ((n-n2_Considered)./(n+n2_Considered)) .* (-1./(n1_Considered+n) + n1_Considered./ (n1_Considered+n) ./ n + 2*n2_Considered ./ (n+n2_Considered) ./ (n-n2_Considered) +i*4*pi*Thickness_Now.*Frequency_Considered/C + (exp(i*4*pi.*Frequency_Considered/C.*n*Thickness_Now) ./ (1-((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/C.*n*Thickness_Now))) .* ((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .*(2*n1_Considered ./ (n+n1_Considered) ./ (n-n1_Considered) + 2*n2_Considered ./ (n+n2_Considered) ./ (n-n2_Considered) + i*4*pi*Thickness_Now.*Frequency_Considered/C)));
                    Current_Loop=1;
                    Number_of_Loop_Checking=25;
                    Center_Wavelength_micron_Index=find(Wavelength_micron_Considered<Center_Wavelength_micron,1,'first');     

                    dn=1E-10;
                    qqqq=1;
                    while (Current_Loop<Number_of_Loop)
                        if mod(Current_Loop,Number_of_Loop_Checking)==0
                            CONT_left=1;
                            CONT_right=1;
                            while(CONT_left || CONT_right)
                                index_needshift_left=find(abs(diff(n(1:Center_Wavelength_micron_Index)))>0.008,1,'last');
                                index_needshift_right=find(abs(diff(n((Center_Wavelength_micron_Index+1):end)))>0.008,1,'first')+Center_Wavelength_micron_Index;
                                CONT_left=0;
                                CONT_right=0;
                                if (index_needshift_left>5)
                                    CONT_left=1;
                                    n(1:index_needshift_left)=n(1:index_needshift_left)-n(index_needshift_left)+n(index_needshift_left+1);
                                elseif (index_needshift_right<(length(n)-5)) 
                                    CONT_right=1;
                                    n((index_needshift_right+1):end)=n((index_needshift_right+1):end)-n(index_needshift_right+1)+n(index_needshift_right);
                                end

                            end

                        end
                        
                        
                        n=n+0.1*((Aexp-A(n,Thickness_Now,Ratio_Lower2Upper_Now))./((A(n+dn,Thickness_Now,Ratio_Lower2Upper_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now))/dn));
                        Current_Loop=Current_Loop+1;
                        n(isnan(n))=n_should;
                        
                        
                        if Current_Loop<58 && Current_Loop>50
                            n_Record(:,qqqq)=n;
                            qqqq=qqqq+1;
                        end
                        
                    end
                        Current_Progress=100*((p-1)+((q-1)+((w-1)+(m-1)/length(Position_Upper2Reference))/length(Ratio_Upper2Reference))/length(Ratio_Lower2Upper))/length(Thickness);
                        disp(sprintf('%f%%',Current_Progress));
                        %Merit=sum(Weight_Function.*(T_Considered-T_abs(n,Thickness_Now)).^2);
                        Merit=sum((T_Considered-T_abs(n,Thickness_Now)).^2);
                        if Merit < Merit_Best
                            Current_Loop_Record=Current_Loop;
                            Merit_Best=Merit;
                            Thickness_Best=Thickness_Now;
                            Ratio_Lower2Upper_Best=Ratio_Lower2Upper_Now;
                            Ratio_Upper2Reference_Best=Ratio_Upper2Reference_Now;
                            Position_Upper2Reference_Best=Position_Upper2Reference_Now;
                            T_Best=T_abs(n,Thickness_Now);
                            n_Best=n;
                        end

                end
            end
        end
    end
end
plot(Wavelength_micron_Considered,real(n_Best),Wavelength_micron_Considered,imag(n_Best));

plot(Wavelength_micron_Considered,T_Best,Wavelength_micron_Considered,T_Considered);
%% Time Frequency Analysis

t_sub=2.*n2./(n2+n1);
r1_Best=(n1_Considered-(n_Best))./(n_Best+n1_Considered);
r1_r_Best=((n_Best)-n1_Considered)./(n_Best+n1_Considered);    
t1_Best=2*(n1_Considered)./(n_Best+n1_Considered);
t1_r_Best=2*(n_Best)./(n_Best+n1_Considered);
t2_Best=2*(n_Best)./(n_Best+n2);
t2_r_Best=2*(n2)./(n_Best+n2);
r2_Best=((n_Best)-n2)./((n_Best)+n2);   
d_Best=exp(i*2*pi.*Frequency_Considered/C.*(n_Best).*Thickness_Best);   %注意! -1*n!
    %d0_Best=exp(i*2*pi.*Frequency/C.*(Distance_fit(j)-Thickness_Now));
    %%注意! -1*n!
%Spectrum_Reference_Considered=Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
%Spectrum_Sample_Considered=Spectrum_Sample(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
%Spectrum_Upper_Temp_Best=(exp(i*4*pi.*Frequency_Considered/C.*(Ratio_Lower2Upper_Best))./r_BK7_Considered./Ratio_Upper2Reference_Best) .* (Position_Upper2Reference_Best.*(n1_Considered-n_Best)./(n1_Considered+n_Best));
%Spectrum_Lower_Temp_Best=(exp(i*4*pi.*Frequency_Considered/C.*(Ratio_Lower2Upper_Best))./r_BK7_Considered./Ratio_Upper2Reference_Best) .* ((exp(i*4*pi.*Frequency_Considered/C.*n_Best*Thickness_Best) ./ (1-((n_Best-n1_Considered)./(n_Best+n1_Considered)) .* ((n_Best-n2_Considered)./(n_Best+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/C.*n_Best*Thickness_Best))) .* ((2*n1_Considered)./(n1_Considered+n_Best)) .* ((2*n_Best)./(n1_Considered+n_Best)) .* ((n_Best-n2_Considered)./(n_Best+n2_Considered)));
%Spectrum_Sample_Upper_Temp_Best=Spectrum_Upper_Temp_Best.*Spectrum_Reference_Considered;
%Spectrum_Sample_Lower_Temp_Best=Spectrum_Lower_Temp_Best.*Spectrum_Reference_Considered;
%plot(Wavelength_micron_Considered,Spectrum_Sample_Upper_Temp_Best,Wavelength_micron_Considered,Spectrum_Sample_Lower_Temp_Best,Wavelength_micron_Considered,Spectrum_Sample_Considered);
%plot(Wavelength_micron_Considered,Spectrum_Upper_Temp_Best,Wavelength_micron_Considered,Spectrum_Lower_Temp_Best,Wavelength_micron_Considered,Aexp);


%Spectrum_Sample_Upper_Temp_Best_fit=interp1(Frequency_Considered,Spectrum_Sample_Upper_Temp_Best,Frequency); 
%Spectrum_Sample_Lower_Temp_Best_fit=interp1(Frequency_Considered,Spectrum_Sample_Lower_Temp_Best,Frequency); 
%Spectrum_Sample_Upper_Temp_Best_fit(isnan(Spectrum_Sample_Upper_Temp_Best_fit))=0;
%Spectrum_Sample_Lower_Temp_Best_fit(isnan(Spectrum_Sample_Lower_Temp_Best_fit))=0;

%Spectrum_Sample_Upper_Temp_Best_fit(N_f+1:N_t)=0;
%Spectrum_Sample_Lower_Temp_Best_fit(N_f+1:N_t)=0;
%Signal_Sample_Upper_Temp_Best_fit=fft(Spectrum_Sample_Upper_Temp_Best_fit);
%Signal_Sample_Lower_Temp_Best_fit=fft(Spectrum_Sample_Lower_Temp_Best_fit);

%plot(Position_micron,Signal_Sample_Upper_Temp_Best_fit,Position_micron,Signal_Sample_Lower_Temp_Best_fit,Position_micron,abs(2*Signal_Sample));
%xlabel('Wavelength (micron)');

%cd('D:\120524\');

%dlmwrite('SAM7 n.txt',real(n_Best),'delimiter','\t','newline','pc');
%dlmwrite('SAM7 k.txt',imag(n_Best),'delimiter','\t','newline','pc');
%dlmwrite('SAM7 T.txt',T_Best,'delimiter','\t','newline','pc');
%dlmwrite('Considered T.txt',T_Considered,'delimiter','\t','newline','pc');
%dlmwrite('Considered Wavelength.txt',Wavelength_micron_Considered,'delimiter','\t','newline','pc');

plot(Wavelength_micron_Considered,real(n_Best),Wavelength_micron_Considered,imag(n_Best));
plot(Wavelength_micron_Considered,A(n_Best,Thickness_Best,Ratio_Lower2Upper_Best),Wavelength_micron_Considered,Aexp);
%plot(Wavelength_micron_Considered,Spectrum_Sample_Upper_Temp_Best,Wavelength_micron_Considered,Spectrum_Sample_Lower_Temp_Best,Wavelength_micron_Considered,Spectrum_Sample_Considered);
plot(Wavelength_micron_Considered,T_Best,Wavelength_micron_Considered,T_Considered);