%% Options

clear all
index_calculation=1;
TF_Analysis=1;  %1: direct in FD, 2: transfer to TD first

Averaging=1;    %Number of frame averaging
Sample_Path='D:\Users\TuanShu\120824_YPR_Down\YPR\';
Reference_Path='D:\Users\TuanShu\120824_YPR_Down\glass\';
Spectroscopy_Path='D:\Users\TuanShu\';

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

Ratio_Upper2Reference=1.1;
Position_Upper2Reference=(15.85)*1E-6;
Number_of_Loop=75;
Thickness=(4.4).*1E-6;% ([-0.1:0.1:0.1]+2.55).*1E-6; In this file, it is the initial guess of thickness
Load_Previous_Result=0;
Ratio_Lower2Upper=1.1;  %([-1:0.02:1]+10.82)*1E-6;%[0.019:0.001:0.024]*1E-6;%[0.02:0.0003:0.023]*1E-6;
%%
n_should=1.7;
delta_n=0.2;
Wavelength_Center=550;

pixel_1=700;               % DC cutoff
%pixel_2=1800;%2000               % the end of 2
%pixel_3=1120;               %for saperation

Center_Wavelength_micron=0.56;
Wavelength_Considered_Min=490;          %nm
Wavelength_Considered_Max=600;

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=N_f*8;

Number_of_variable=3;           %Wavelength indep. variables


%% Global arrays generation

c=3E8;

Max_Frequency=c/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=c/(Max_Wavelength*1E-9);             %Hz

Frequency_Center=c/(Wavelength_Center*1E-9);
Frequency_Considered_Min=c/(Wavelength_Considered_Max*1E-9);             %Hz
Frequency_Considered_Max=c/(Wavelength_Considered_Min*1E-9);             %Hz


cd(sprintf('%sinter\\',Sample_Path));
Data=importdata('D0.txt');      % Data_2: the glass data 1

Wavelength=Data(:,1);           %nm
Frequency_Old=c./(Wavelength*1E-9);
Frequency=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency=Frequency';
Wavelength_micron=(c./Frequency)*1E6;

Frequency_Center_Index=find(Frequency>Frequency_Center,1,'first');
Frequency_Considered_Min_Index=find(Frequency>Frequency_Considered_Min,1,'first');
Frequency_Considered_Max_Index=find(Frequency>Frequency_Considered_Max,1,'first');

% Time-domain

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=c*Time;
Position_micron=Position*1E6;

%% T 


cd(sprintf('%s\\',Spectroscopy_Path));
Data_Spectroscopy=importdata('1208030331_5micron 1.jws.txt');

Spectrum_Spectroscopy_Old=Data_Spectroscopy(end:-1:1,2)/100;
Wavelength_Spectroscopy=Data_Spectroscopy(end:-1:1,1);     
Frequency_Spectroscopy=c./(Wavelength_Spectroscopy*1E-9);

T=interp1(Frequency_Spectroscopy,Spectrum_Spectroscopy_Old,Frequency,'spline'); 

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
n1=n_bk7;%/1.516781257666726*1.520;

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

Spectrum_Sample_Wavelength=Filter_inter(:,2)-Filter_sam(:,2)-Filter_ref(:,2)+Filter_bs(:,2);
Spectrum_Reference_Wavelength=Glass_inter(:,2)-Glass_sam(:,2)-Glass_ref(:,2)+Glass_bs(:,2);
%Spectrum_Sample_Wavelength=ref(:,2);

Spectrum_Sample_Wavelength=Spectrum_Sample_Wavelength-Spectrum_Sample_Wavelength(1);
Spectrum_Reference_Wavelength=Spectrum_Reference_Wavelength-Spectrum_Reference_Wavelength(1);

%plot(Wavelength,Spectrum_Sample_Wavelength);

%xlabel('Wavelength (micron)');
%ylabel('Spectral Power (a.u.)');

%% Spectrum generation

Separation_Position=37.62;   %(micron)
Gauss_Window_sig=0.2;
Gauss_Window_C=29;

Gauss_Window_Profile=1-gaussmf(Position_micron,[Gauss_Window_sig Gauss_Window_C]);
Gauss_Window_Profile(Position_micron<Gauss_Window_C)=0;


Gauss_Window_C_2=48;

Gauss_Window_Profile_2=1-gaussmf(Position_micron,[Gauss_Window_sig Gauss_Window_C_2]);
Gauss_Window_Profile_2(Position_micron>Gauss_Window_C_2)=0;

Spectrum_Sample_Frequency=(Spectrum_Sample_Wavelength.*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);
Spectrum_Sample=interp1(Frequency_Old,Spectrum_Sample_Frequency,Frequency,'spline'); 
Spectrum_Sample(isnan(Spectrum_Sample))=0;
Spectrum_Sample(Frequency<Min_Frequency)=0;
Spectrum_Sample(N_f+1:N_t)=0;
Signal_Sample=fft(Spectrum_Sample).*Gauss_Window_Profile.*Gauss_Window_Profile_2;  
Signal_Sample(round(length(Signal_Sample)/2)+1:end)=0;

Gauss_Window_C_ref=13;

Gauss_Window_Profile_ref=1-gaussmf(Position_micron,[Gauss_Window_sig Gauss_Window_C_ref]);
Gauss_Window_Profile_ref(Position_micron<Gauss_Window_C_ref)=0;


Gauss_Window_C_2_ref=23;

Gauss_Window_Profile_2_ref=1-gaussmf(Position_micron,[Gauss_Window_sig Gauss_Window_C_2_ref]);
Gauss_Window_Profile_2_ref(Position_micron>Gauss_Window_C_2_ref)=0;

Spectrum_Reference_Frequency=(Spectrum_Reference_Wavelength.*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);
Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency,'spline'); 
Spectrum_Reference(isnan(Spectrum_Reference))=0;
Spectrum_Reference(Frequency<Min_Frequency)=0;
Spectrum_Reference(N_f+1:N_t)=0;
Signal_Reference=fft(Spectrum_Reference).*Gauss_Window_Profile_ref.*Gauss_Window_Profile_2_ref;  
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

Spectrum_Devided=Spectrum_Sample_1./Spectrum_Reference;
Spectrum_Devided_2=Spectrum_Sample_2./Spectrum_Reference;
Spectrum_Devided_3=Spectrum_Sample_2./Spectrum_Sample_1;

Phase_Sample=angle(Spectrum_Sample);

plot(Position_micron,Signal_Sample_1,Position_micron,Signal_Sample_2);
plot(Position_micron,Signal_Sample_1,Position_micron,Signal_Sample_2,Position_micron,Signal_Reference);
xlabel('Position (micron)');
ylabel('Interference Power (a.u.)');


plot(Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Devided(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Devided_2(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),(abs(Spectrum_Devided_2(Wavelength_micron<1)./abs(Spectrum_Devided(Wavelength_micron<1)))));
xlabel('Wavelength (micron)');
ylabel('Interference Power spectral density (a.u.)');


Frequency_Considered=Frequency(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
Aexp=abs(Spectrum_Devided(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index));   %注意! -1*n!;
Bexp=abs(Spectrum_Devided_2(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index));   %注意! -1*n!;
Cexp=angle(Spectrum_Devided_3(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index));   %注意! -1*n!;
T_Considered=T(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
%Weight_Function=(abs(Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index)).*T_Considered)./max((abs(Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index).*T_Considered)));
Wavelength_micron_Considered=Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
Center_Wavelength_micron_Index=find(Wavelength_micron_Considered<Center_Wavelength_micron,1,'first');  
n1_Considered=n1(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
n2_Considered=n2;
r_BK7_Considered=r_BK7(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
t_AN100_Considered=t_AN100(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);


Cexp=unwrap(Cexp);

plot(Wavelength_micron_Considered,Aexp,Wavelength_micron_Considered,Bexp);
xlabel('Wavelength (micron)');
ylabel('Normalized Inteference Spectrum');


%Cexp=Cexp-Cexp(Center_Wavelength_micron_Index);

    Phase_diff=Cexp(Center_Wavelength_micron_Index);
    NN=fix(Phase_diff/(2*pi));
    Cexp=Cexp-NN*2*pi;
    
    
plot(Wavelength_micron_Considered,Cexp);
xlabel('Wavelength (micron)');
ylabel('Continuous Phase Spectrum');


%%

t1= @ (n) (2*n1_Considered)./(n1_Considered+n); 
t1_r= @ (n) (2*n)./(n2_Considered+n);
t2= @ (n) (2*n)./(n+n2_Considered); 
r1= @ (n) (n1_Considered-n)./(n1_Considered+n);
r1_r= @ (n) (n-n1_Considered)./(n+n1_Considered);
r2= @ (n) (n-n2_Considered)./(n+n2_Considered);
                    %A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now) exp(i*4*pi.*Frequency_Considered.*Position_Upper2Reference_Now/c).*Ratio_Upper2Reference_Now./r_BK7_Considered.*(r1(n)+Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(i*2.*Frequency_Considered.*n.*Thickness_Now/c)./(1-r1_r(n).*r2(n).*exp(i*2*Frequency_Considered.*n.*Thickness_Now/c))));
                    %A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now) abs(exp(i*4*pi.*Frequency_Considered.*Position_Upper2Reference_Now/c).*Ratio_Upper2Reference_Now./r_BK7_Considered.*(r1(n)));
                    %B = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now) abs(exp(i*4*pi.*Frequency_Considered.*Position_Upper2Reference_Now/c).*Ratio_Upper2Reference_Now./r_BK7_Considered.*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(i*2.*Frequency_Considered.*n.*Thickness_Now/c))));
A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) abs(Ratio_Upper2Reference_Now./r_BK7_Considered.*(r1(n)));
B = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) abs(Ratio_Upper2Reference_Now./r_BK7_Considered.*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c))));
C = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) unwrap(angle(1./r1(n).*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c)))));
                   
det_array = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9) Q1.*Q5.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 -Q3.*Q5.*Q7 -Q2.*Q4.*Q9 -Q1.*Q6.*Q8;
                    %A = @ (n,Ratio_Upper2Reference_Now,Position_Upper2Reference) exp(i*4*pi.*Frequency_Considered.*Position_Upper2Reference_Now/c).*Ratio_Upper2Reference_Now./r_BK7_Considered.*r1(n);
                    
                    %A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now) Ratio_Lower2Upper_Now./r1(n).*(t1(n).*t1_r(n).*r2(n).*exp(i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c)./(1-r1_r(n).*r2(n).*exp(i*4*pi*Frequency_Considered.*n.*Thickness_Now/c)));                    
                    %T_abs= @ (n,Thickness_Now) abs((2.*n2_Considered./(n2_Considered+n1_Considered)).*(2*n1_Considered./(n+n1_Considered)).*(2*n./(n+n2_Considered)).*exp(i*2*pi.*Frequency_Considered/c.*n.*Thickness_Now)./(1-((n-n1_Considered)./(n+n1_Considered)).*((n-n2_Considered)./(n+n2_Considered)).*exp(i*4*pi.*Frequency_Considered/c.*n.*Thickness_Now))).^2;
T_abs = @ (n,Thickness_Now) abs(t_AN100_Considered.*t1(n).*t2(n).*exp(i*2*pi*Frequency_Considered.*n.*Thickness_Now/c)./(1-r1_r(n).*r2(n).*exp(i*4*pi*Frequency_Considered.*n.*Thickness_Now/c))).^2;          %Note, 120827 remove 2 from 分子4*pi > 2*pi
                    
dn=1E-10;
dt=1E-16;

n(1:length(Frequency_Considered),1)=n_should;

Thickness_Now(1:length(Frequency_Considered),1)=Thickness;
%plot(Wavelength_micron_Considered,Aexp,Wavelength_micron_Considered,A(n,Thickness,Ratio_Lower2Upper(1),Ratio_Upper2Reference(1)));
%plot(Wavelength_micron_Considered,Bexp,Wavelength_micron_Considered,B(n,Th
%ickness,Ratio_Lower2Upper(1),Ratio_Upper2Reference(1)));
%plot(Wavelength_micron_Considered,Cexp,Wavelength_micron_Considered,C(n,35E-6,Ratio_Lower2Upper(1),Ratio_Upper2Reference(1)));


%% Start the n k fitting    
if index_calculation==1   
    Merit_Best=999999999999999999999;
    for p=1:length(Thickness)
        for q=1:length(Ratio_Lower2Upper)
            for w=1:length(Ratio_Upper2Reference)
                for m=1:length(Position_Upper2Reference)
                    

                    %Thickness_Now=Thickness(p);
                    Ratio_Lower2Upper_Now=Ratio_Lower2Upper(q);
                    Ratio_Upper2Reference_Now=Ratio_Upper2Reference(w);
                    Position_Upper2Reference_Now=Position_Upper2Reference(m);
                            
                    n(1:length(Frequency_Considered),1)=n_should;%1.4+(Wavelength_micron_Considered-0.6)*2.8-(Wavelength_micron_Considered-0.7)*1.4;
                    %A = @
                    %(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now) (exp(i*4*pi.*Frequency_Considered/c.*(Ratio_Lower2Upper_Now))./r_BK7_Considered./Ratio_Upper2Reference_Now) .* (Position_Upper2Reference_Now.*(n1_Considered-n)./(n1_Considered+n)+(exp(i*4*pi.*Frequency_Considered/c.*n*Thickness_Now) ./ (1-((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/c.*n*Thickness_Now))) .* ((2*n1_Considered)./(n1_Considered+n)) .* ((2*n)./(n1_Considered+n)) .* ((n-n2_Considered)./(n+n2_Considered)));
                    %dA_dn = @ (n) (exp(i*4*pi.*Frequency_Considered/c.*(Ratio_Lower2Upper))./r_BK7_Considered./Ratio_Upper2Reference_Now) .* (Position_Upper2Reference_Now.*(-2*n1_Considered./(n1_Considered+n).^2) + (exp(i*4*pi.*Frequency_Considered/c.*n*Thickness_Now) ./ (1-((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/c.*n*Thickness_Now))) .* ((2*n1_Considered)./(n1_Considered+n)) .* ((2*n)./(n1_Considered+n)) .* ((n-n2_Considered)./(n+n2_Considered)) .* (-1./(n1_Considered+n) + n1_Considered./ (n1_Considered+n) ./ n + 2*n2_Considered ./ (n+n2_Considered) ./ (n-n2_Considered) +i*4*pi*Thickness_Now.*Frequency_Considered/c + (exp(i*4*pi.*Frequency_Considered/c.*n*Thickness_Now) ./ (1-((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/c.*n*Thickness_Now))) .* ((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .*(2*n1_Considered ./ (n+n1_Considered) ./ (n-n1_Considered) + 2*n2_Considered ./ (n+n2_Considered) ./ (n-n2_Considered) + i*4*pi*Thickness_Now.*Frequency_Considered/c)));
                    Current_Loop=1;
                    Number_of_Loop_Checking=25;   
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
                        dA_dn_real = (A(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dB_dn_real = (B(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dC_dn_real = (C(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                    
                        dA_dn_imag = (A(n+i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dB_dn_imag = (B(n+i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dC_dn_imag = (C(n+i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        
                        
                        dA_dthickness = (A(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        dB_dthickness = (B(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        dC_dthickness = (C(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        
                        C_Temp=C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        %C_Temp=C_Temp-C_Temp(Center_Wavelength_micron_Index);
                        
                        
                            Phase_diff=C_Temp(Center_Wavelength_micron_Index);
                            NN=fix(Phase_diff/(2*pi));
                            C_Temp=C_Temp-NN*2*pi;
                        
                        delta_A=Aexp-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        delta_B=Bexp-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        delta_C=Cexp-C_Temp;
                        
                        %dn_real=(delta_A.*dB_dn_imag-dA_dn_imag.*delta_B)./(dA_dn_real.*dB_dn_imag-dA_dn_imag.*dB_dn_real);
                        %dn_imag=(dA_dn_real.*delta_B-delta_A.*dB_dn_real)./(dA_dn_real.*dB_dn_imag-dA_dn_imag.*dB_dn_real);
                        
                        dn_real = det_array(delta_A,dA_dn_imag,dA_dthickness,delta_B,dB_dn_imag,dB_dthickness,delta_C,dC_dn_imag,dC_dthickness)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        dn_imag = det_array(dA_dn_real,delta_A,dA_dthickness,dB_dn_real,delta_B,dB_dthickness,dC_dn_real,delta_C,dC_dthickness)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        dthickness = det_array(dA_dn_real,dA_dn_imag,delta_A,dB_dn_real,dB_dn_imag,delta_B,dC_dn_real,dC_dn_imag,delta_C)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        
                        % SLOWER!!! why??? dn_real=(Aexp-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now).*(B(n+i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn-(A(n+i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn.*Bexp-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))./((A(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn.*(B(n+i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn-(A(n+i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn.*(B(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn);
                        %dn_imag=((A(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn.*Bexp-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-Aexp-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now).*(B(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn)./((A(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn.*(B(n+i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn-(A(n+i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn.*(B(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn);
                        
                        n=n+0.1*(dn_real+i*dn_imag);
                        Thickness_Now=Thickness_Now+0.1*((dthickness));
                        
                        %n=n+0.1*((Aexp-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now))./((A(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now))/dn))+          0.1*((Bexp-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now))./((B(n+  i*  dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now,Position_Upper2Reference_Now))/  (i*dn)   ));
                        Current_Loop=Current_Loop+1;
                        n(isnan(n))=n_should;
                        n(isinf(n))=n_should;
                        
                        if Current_Loop>70
                            n_Record(:,qqqq)=n;
                            Thickness_Record(:,qqqq)=Thickness_Now;
                            qqqq=qqqq+1;
                        end
                        
                    end
                        Current_Progress=100*((p-1)+((q-1)+((w-1)+(m-1)/length(Position_Upper2Reference))/length(Ratio_Upper2Reference))/length(Ratio_Lower2Upper))/length(Thickness);
                        disp(sprintf('%f%%',Current_Progress));
                        %Merit=sum(Weight_Function.*(T_Considered-T_abs(n,Thickness_Now)).^2);
                        Merit=sum((T_Considered-T_abs(n,Thickness_Now)).^2);
                        %Merit=std(dthickness);
                        if Merit < Merit_Best
                            Current_Loop_Record=Current_Loop;
                            Merit_Best=Merit;
                            Thickness_Best=mean(Thickness_Now);
                            Thickness_std=std(Thickness_Now);
                            Thickness_Record=Thickness_Now;
                            dthickness_Best=dthickness;
                            Ratio_Lower2Upper_Best=Ratio_Lower2Upper_Now;
                            Ratio_Upper2Reference_Best=Ratio_Upper2Reference_Now;
                            Position_Upper2Reference_Best=Position_Upper2Reference_Now;
                            T_Best=T_abs(n,Thickness_Now);
                            C_Best=C_Temp;
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
d_Best=exp(i*2*pi.*Frequency_Considered/c.*(n_Best).*Thickness_Best);   %注意! -1*n!
    %d0_Best=exp(i*2*pi.*Frequency/c.*(Distance_fit(j)-Thickness_Now));
    %%注意! -1*n!
%Spectrum_Reference_Considered=Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
%Spectrum_Sample_Considered=Spectrum_Sample(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
%Spectrum_Upper_Temp_Best=(exp(i*4*pi.*Frequency_Considered/c.*(Ratio_Lower2Upper_Best))./r_BK7_Considered./Ratio_Upper2Reference_Best) .* (Position_Upper2Reference_Best.*(n1_Considered-n_Best)./(n1_Considered+n_Best));
%Spectrum_Lower_Temp_Best=(exp(i*4*pi.*Frequency_Considered/c.*(Ratio_Lower2Upper_Best))./r_BK7_Considered./Ratio_Upper2Reference_Best) .* ((exp(i*4*pi.*Frequency_Considered/c.*n_Best*Thickness_Best) ./ (1-((n_Best-n1_Considered)./(n_Best+n1_Considered)) .* ((n_Best-n2_Considered)./(n_Best+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/c.*n_Best*Thickness_Best))) .* ((2*n1_Considered)./(n1_Considered+n_Best)) .* ((2*n_Best)./(n1_Considered+n_Best)) .* ((n_Best-n2_Considered)./(n_Best+n2_Considered)));
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
%plot(Wavelength_micron_Considered,Spectrum_Sample_Upper_Temp_Best,Wavelength_micron_Considered,Spectrum_Sample_Lower_Temp_Best,Wavelength_micron_Considered,Spectrum_Sample_Considered);
%plot(Wavelength_micron_Considered,T_Best,Wavelength_micron_Considered,T_Considered);
plot(Wavelength_micron_Considered,Cexp,Wavelength_micron_Considered,C_Best);
plot(Wavelength_micron_Considered,Aexp,Wavelength_micron_Considered,A(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best));
plot(Wavelength_micron_Considered,Bexp,Wavelength_micron_Considered,B(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best));
plot(Wavelength_micron_Considered,Thickness_Best+dthickness_Best-mean(dthickness_Best));
xlabel('Wavelength (micron)');
ylabel('Calculated Thickness (m)');

plot(Wavelength_micron_Considered*1000,real(n_Best),Wavelength_micron_Considered*1000,n1_Considered);
xlabel('Wavelength (nm)');
ylabel('n');
plot(Wavelength_micron_Considered*1000,imag(n_Best));
xlabel('Wavelength (nm)');
ylabel('k');

plot(Wavelength_micron_Considered,T_Best,Wavelength_micron_Considered,T_Considered);
xlabel('Wavelength (micron)');
ylabel('T');


Thickness_Best
Thickness_std
dlmwrite('n.txt',real(n_Best),'delimiter','\t','newline','pc');
dlmwrite('k.txt',imag(n_Best),'delimiter','\t','newline','pc');
dlmwrite('T_Best.txt',T_Best,'delimiter','\t','newline','pc');

dlmwrite('T_Considered.txt',T_Considered,'delimiter','\t','newline','pc');
dlmwrite('Wavelength_micron_Considered.txt',Wavelength_micron_Considered,'delimiter','\t','newline','pc');


Signal_Sample_1_Real=real(Signal_Sample_1);
Signal_Sample_2_Real=real(Signal_Sample_2);

dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

dlmwrite('Signal_Sample_1_Real.txt',Signal_Sample_1_Real,'delimiter','\t','newline','pc');
dlmwrite('Signal_Sample_2_Real.txt',Signal_Sample_2_Real,'delimiter','\t','newline','pc');

%plot(Wavelength_micron_Considered,A(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best),Wavelength_micron_Considered,Aexp);
%plot(Wavelength_micron_Considered,B(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best),Wavelength_micron_Considered,Bexp);
%plot(Wavelength_micron_Consi