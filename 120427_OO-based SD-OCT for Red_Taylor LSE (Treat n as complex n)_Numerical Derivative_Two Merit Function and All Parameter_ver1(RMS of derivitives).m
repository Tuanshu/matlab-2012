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

Ratio_Reference2Sample=1;
Ratio_Upper2Lower=1;
Number_of_Loop=100;
Thickness=2.45E-6;% ([-0.1:0.1:0.1]+2.55).*1E-6;
Load_Previous_Result=0;
Position_Upper_Interface=10.82E-6;  %([-1:0.02:1]+10.82)*1E-6;%[0.019:0.001:0.024]*1E-6;%[0.02:0.0003:0.023]*1E-6;
%%
n_should=1.8;
delta_n=0.4;
Wavelength_Center=600;

Scanning_1='Thickness';
Scanning_2='Position_Upper_Interface';


if strcmp(Scanning_1,'n_should')
    Parameter_1=n_should;
elseif strcmp(Scanning_1,'Thickness')
    Parameter_1=Thickness;
elseif strcmp(Scanning_1,'Position_Upper_Interface')
    Parameter_1=Position_Upper_Interface;
elseif strcmp(Scanning_1,'Ratio_Reference2Sample')
    Parameter_1=Ratio_Reference2Sample;
elseif strcmp(Scanning_1,'Ratio_Upper2Lower')
    Parameter_1=Ratio_Upper2Lower;
    
end

if strcmp(Scanning_2,'n_should')
    Parameter_2=n_should;
elseif strcmp(Scanning_2,'Thickness')
    Parameter_2=Thickness;
elseif strcmp(Scanning_2,'Position_Upper_Interface')
    Parameter_2=Position_Upper_Interface;
elseif strcmp(Scanning_2,'Ratio_Reference2Sample')
    Parameter_2=Ratio_Reference2Sample;
elseif strcmp(Scanning_2,'Ratio_Upper2Lower')
    Parameter_2=Ratio_Upper2Lower;
end

pixel_1=700;               % DC cutoff
%pixel_2=1800;%2000               % the end of 2
%pixel_3=1120;               %for saperation

Wavelength_Considered_Min=480;          %nm
Wavelength_Considered_Max=660;

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

%% T 


cd('D:\120418\');
Data_Spectroscopy=importdata('red.jws.txt');

Spectrum_Spectroscopy_Old=Data_Spectroscopy(end:-1:1,2)/100;
Wavelength_Spectroscopy=Data_Spectroscopy(end:-1:1,1);     
Frequency_Spectroscopy=C./(Wavelength_Spectroscopy*1E-9);

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

%plot(Wavelength,Spectrum_Sample_Wavelength);

%xlabel('Wavelength (micron)');
%ylabel('Spectral Power (a.u.)');

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

%plot(Position_micron,Signal_Sample,Position_micron,Signal_Reference);

%xlabel('Optical path Difference (micron)');
%ylabel('Interference Signal');

%plot(Wavelength_micron(Wavelength_micron<0.80),real(Spectrum_Reference(Wavelength_micron<0.80,:)),Wavelength_micron(Wavelength_micron<0.80),real(Spectrum_Sample(Wavelength_micron<0.80,:)));


%plot(Wavelength_micron(Wavelength_micron<0.80),real(Spectrum_Devided(Wavelength_micron<0.80,:)),Wavelength_micron(Wavelength_micron<0.80),abs(Spectrum_Devided(Wavelength_micron<0.80,:)));

%xlabel('Wavelength (micron)');
%ylabel('A');

%% to generate the phase spectrum


%% Start the n k fitting    

if index_calculation==1   
            
            Frequency_Considered=Frequency(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
            Aexp=Spectrum_Devided(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);   %注意! -1*n!;
            T_Considered=T(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
            Wavelength_micron_Considered=Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
            n1_Considered=n1(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
            n2_Considered=n2;
            r_BK7_Considered=r_BK7(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
            n(1:length(Frequency_Considered),1)=1.6;%1.4+(Wavelength_micron_Considered-0.6)*2.8-(Wavelength_micron_Considered-0.7)*1.4;
            A = @ (n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower) (exp(i*4*pi.*Frequency_Considered/C.*(Ratio_Reference2Sample))./r_BK7_Considered./Ratio_Reference2Sample) .* (Ratio_Upper2Lower.*(n1_Considered-n)./(n1_Considered+n)+(exp(i*4*pi.*Frequency_Considered/C.*n*Thickness) ./ (1-((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/C.*n*Thickness))) .* ((2*n1_Considered)./(n1_Considered+n)) .* ((2*n)./(n1_Considered+n)) .* ((n-n2_Considered)./(n+n2_Considered)));
            T_abs= @ (n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower) abs((2.*n2_Considered./(n2_Considered+n1_Considered)).*(2*n1_Considered./(n+n1_Considered)).*(2*n./(n+n2_Considered)).*exp(i*2*pi.*Frequency_Considered/C.*n.*Thickness)./(1-((n-n1_Considered)./(n+n1_Considered)).*((n-n2_Considered)./(n+n2_Considered)).*exp(i*4*pi.*Frequency_Considered/C.*n.*Thickness))).^2;
            
            
            %dA_dn = @ (n) (exp(i*4*pi.*Frequency_Considered/C.*(Ratio_Reference2Sample))./r_BK7_Considered./Ratio_Reference2Sample) .* (Ratio_Upper2Lower.*(-2*n1_Considered./(n1_Considered+n).^2) + (exp(i*4*pi.*Frequency_Considered/C.*n*Thickness) ./ (1-((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/C.*n*Thickness))) .* ((2*n1_Considered)./(n1_Considered+n)) .* ((2*n)./(n1_Considered+n)) .* ((n-n2_Considered)./(n+n2_Considered)) .* (-1./(n1_Considered+n) + n1_Considered./ (n1_Considered+n) ./ n + 2*n2_Considered ./ (n+n2_Considered) ./ (n-n2_Considered) +i*4*pi*Thickness.*Frequency_Considered/C + (exp(i*4*pi.*Frequency_Considered/C.*n*Thickness) ./ (1-((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .* exp(i*4*pi.*Frequency_Considered/C.*n*Thickness))) .* ((n-n1_Considered)./(n+n1_Considered)) .* ((n-n2_Considered)./(n+n2_Considered)) .*(2*n1_Considered ./ (n+n1_Considered) ./ (n-n1_Considered) + 2*n2_Considered ./ (n+n2_Considered) ./ (n-n2_Considered) + i*4*pi*Thickness.*Frequency_Considered/C)));
            Current_Loop=1;
            Number_of_Loop_Checking=25;
            Center_Wavelength_micron=0.56;
            Center_Wavelength_micron_Index=find(Wavelength_micron_Considered<Center_Wavelength_micron,1,'first');     
            
            Parameter_Array=[Thickness Position_Upper_Interface Ratio_Reference2Sample Ratio_Upper2Lower];
            
            dn=1E-10;
            dThickness=1E-16;
            dPosition_Upper_Interface=1E-16;
            dRatio_Reference2Sample=1E-10;
            dRatio_Upper2Lower=1E-10;            
            while (Current_Loop<Number_of_Loop)
                if mod(Current_Loop,Number_of_Loop_Checking)==0
                    CONT_left=1;
                    CONT_right=1;
                    while(CONT_left || CONT_right)
                        index_needshift_left=find(abs(diff(n(1:Center_Wavelength_micron_Index)))>0.01,1,'last');
                        index_needshift_right=find(abs(diff(n((Center_Wavelength_micron_Index+1):end)))>0.01,1,'first')+Center_Wavelength_micron_Index;
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
                n=n+0.1*((Aexp-A(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower))./((A(n+dn,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower)-A(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower))/dn));
                Current_Loop=Current_Loop+1;
                Thickness=Thickness+0.1*(((sum(T_Considered-T_abs(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower)).^2).^0.5)./((sum((T_abs(n,Thickness+dThickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower)-T_abs(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower))/dThickness).^2).^0.5));
                Position_Upper_Interface=Position_Upper_Interface+0.1*(((sum(T_Considered-T_abs(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower)).^2).^0.5)./((sum((T_abs(n,Thickness,Position_Upper_Interface+dPosition_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower)-T_abs(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower))/dPosition_Upper_Interface).^2).^0.5));
                Ratio_Reference2Sample=Ratio_Reference2Sample+0.1*(((sum(T_Considered-T_abs(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower)).^2).^0.5)./((sum((T_abs(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample+dRatio_Reference2Sample,Ratio_Upper2Lower)-T_abs(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower))/dRatio_Reference2Sample).^2).^0.5));
                Ratio_Upper2Lower=Ratio_Upper2Lower+0.1*(((sum(T_Considered-T_abs(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower)).^2).^0.5)./((sum((T_abs(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower+dRatio_Upper2Lower)-T_abs(n,Thickness,Position_Upper_Interface,Ratio_Reference2Sample,Ratio_Upper2Lower))/dRatio_Upper2Lower).^2).^0.5));
                n(isnan(n))=0;
                if isnan(Thickness)
                    Thickness=0;
                elseif isnan(Position_Upper_Interface)
                    Position_Upper_Interface=0;
                elseif isnan(Ratio_Reference2Sample)
                    Ratio_Reference2Sample=0;
                elseif isnan(Ratio_Upper2Lower)
                    Ratio_Upper2Lower=0;
                end
            
            end
               
end
plot(Wavelength_micron_Considered,real(n),Wavelength_micron_Considered,imag(n));

%% Time Frequency Analysis
backup=0;
if backup == 1

clear TFD
% Remark: for 500~600nm, it's 6E14 to 5E14 Hz, 1E14 Hz difference
Frequency_TF_Min=4E14;
Frequency_TF_Max=6.5E14;
Spectral_Resolution=10E13;  %Use in method 1

Position_TF_micron_min=5;
Position_TF_micron_max=25;
Spatial_Resolution=0.25;  %Use in method 2, micron

Frequency_TF=Frequency_TF_Min:Spectral_Resolution:Frequency_TF_Max;
Position_TF=Position_TF_micron_min:Spatial_Resolution:Position_TF_micron_max;

Frequency_TF_Min_index=find(Frequency>Frequency_TF_Min,1,'first');
Frequency_TF_Max_index=find(Frequency>Frequency_TF_Max,1,'first');
Position_TF_micron_min_index=find(Position_micron>Position_TF_micron_min,1,'first');
Position_TF_micron_max_index=find(Position_micron>Position_TF_micron_max,1,'first');

Wavelength_TF_micron=C./Frequency_TF*1E6;

Spectrum_Input=Spectrum_Devided;

if TF_Analysis == 1             % Direct calculation in FD

TFD(1:length(Position_micron),1:length(Frequency_TF))=0;

for j=1:length(Frequency_TF)
    Spectrum_Windowed=Spectrum_Input.*gaussmf(Frequency,[Spectral_Resolution Frequency_TF(j)]);
    Spectrum_Windowed(N_f+1:N_t)=0;
    TFD(:,j)=fft(Spectrum_Windowed);  
end

imagesc(abs(TFD(Position_TF_micron_min_index:Position_TF_micron_max_index,:)),'xdata',Wavelength_TF_micron,'ydata',Position_micron(Position_TF_micron_min_index:Position_TF_micron_max_index));

elseif TF_Analysis == 2         % Calculation in TD

TFD(1:length(Position_TF),1:length(Frequency))=0;    
    
Spectrum_Input_Temp=Spectrum_Input;
Spectrum_Input_Temp(N_f+1:N_t)=0;
Signal_Input=fft(Spectrum_Input_Temp);

for j=1:length(Position_TF)

    Spectrum_Windowed=ifft(Signal_Input.*gaussmf(Position_micron,[Spatial_Resolution Position_TF(j)]));
    TFD(j,:)=Spectrum_Windowed(1:N_f);
end    

imagesc(abs(TFD(:,Frequency_TF_Min_index:Frequency_TF_Max_index)),'xdata',Wavelength_micron(Frequency_TF_Min_index:Frequency_TF_Max_index),'ydata',Position_TF);

end
%% Fitting

if index_calculation == 1
    
for p=1:size(k_final,2)
    for q=1:size(k_final,1)
        if q>1
            if k_final(q,p)>0.08
                k_final(q,p)=k_final(q-1,p);
            end
        end
    end
end
    
%% Checking
T_abs_check(1:N_f,1:length(Parameter))=0;
for j=1:length(Parameter)
    if index_calculation == 1        
    t_sub=2.*n2./(n2+n1);
    r1_check=(n1(index_now)-(n_final(:,j)+i*k_final(:,j)))./(n_final(:,j)+n1(index_now)+i*k_final(:,j));
    r1_r_check=((n_final(:,j)+i*k_final(:,j))-n1(index_now))./(n_final(:,j)+n1(index_now)+i*k_final(:,j));    
    t1_check=2*(n1(index_now))./(n_final(:,j)+n1(index_now)+i*k_final(:,j));
    t1_r_check=2*(n_final(:,j)+i*k_final(:,j))./(n_final(:,j)+n1(index_now)+i*k_final(:,j));
    t2_check=2*(n_final(:,j)+i*k_final(:,j))./(n_final(:,j)+n2+i*k_final(:,j));
    t2_r_check=2*(n2)./(n_final(:,j)+n2+i*k_final(:,j));
    r2_check=((n_final(:,j)+i*k_final(:,j))-n2)./((n_final(:,j)+i*k_final(:,j))+n2);   
    d_n_check=exp(i*2*pi.*Frequency/C.*(n_final(:,j)).*Thickness);   %注意! -1*n!
    d_check=exp(i*2*pi.*Frequency/C.*(n_final(:,j)+i*k_final(:,j)).*Thickness);   %注意! -1*n!
    %d0_check=exp(i*2*pi.*Frequency/C.*(Distance_fit(j)-Thickness));   %注意! -1*n!

    Spectrum_Upper_Temp_Check=((r1_check)./r_BK7);                       %神說是cosine
    Spectrum_Lower_Temp_Check=((t1_check.*t1_r_check.*r2_check.*(d_check.^2))./r_BK7); 
    T_abs_check(:,j)=abs(t_sub.*t1_check.*t2_check.*(d_check)./(1-r2_check.*r1_r_check.*(d_check).^2)).^2;
    %T_abs=abs((2.*n2./(n2+n1)).*(2*n1./(n+n1)).*(2*n./(n+n2)).*exp(i*2*pi.*Frequency/C.*n.*Thickness)./(1-((n-n1)./(n+n1)).*((n-n2)./(n+n2)).*exp(i*4*pi.*Frequency/C.*n.*Thickness))).^2;
    %T_abs_check(:,j)=abs(t_sub.*t1_check.*t2_check.*(d_check)).^2;
    QQ=abs(r2_check.*r1_r_check.*(d_check).^2);
    end
end
    


%dlmwrite('Wavelength_micron.txt',Wavelength_micron,'delimiter','\t','newline','pc','precision','%.12f');
plot(Wavelength_micron(Wavelength_micron<0.65),n_final(Wavelength_micron<0.65,:));
plot(Wavelength_micron(Wavelength_micron<0.65),k_final(Wavelength_micron<0.65,:));

plot(Wavelength_micron(Wavelength_micron<0.65),T_abs_check(Wavelength_micron<0.65,:));

dlmwrite('n_final.txt',n_final,'delimiter','\t','newline','pc');

dlmwrite('k_final.txt',k_final,'delimiter','\t','newline','pc');

dlmwrite('T_abs_check.txt',T_abs_check,'delimiter','\t','newline','pc');

dlmwrite('Thickness.txt',Thickness,'delimiter','\t','newline','pc');

dlmwrite('Wavelength_micron.txt',Wavelength_micron,'delimiter','\t','newline','pc');

dlmwrite('Frequency.txt',Frequency,'delimiter','\t','newline','pc');
%%

plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),n_final(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,2));

xlabel('Wavelength (micron)');
ylabel('n');

%%
plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),k_final(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1));

xlabel('Wavelength (micron)');
ylabel('k');
%%
plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1));

%plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1)/max(T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,5)));

xlabel('Wavelength (micron)');
ylabel('T ');


%plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1)/max(T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1)),Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),QQ(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index)/max(QQ(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index)));
%%
VIEW_INDEX=1;
%xlabel('Wavelength (micron)');
%ylabel('T (percent)');
plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),Spectrum_Upper_Record(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,VIEW_INDEX),Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),Spectrum_Lower_Record(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,VIEW_INDEX));
%plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),Spectrum_Upper_Record(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1)+Spectrum_Lower_Record(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1));
%plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1)/max(T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,5)));

plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),real(Spectrum_Devided(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,:)),Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),Spectrum_Upper_Record(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,VIEW_INDEX)+Spectrum_Lower_Record(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,VIEW_INDEX),Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),Spectrum_Upper_Record(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,VIEW_INDEX),Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),Spectrum_Lower_Record(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,VIEW_INDEX));
xlabel('Wavelength (micron)');

%plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),d_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1));
%%

end
end