%% Options

clear all
index_calculation=1;

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
Number_of_Loop=5;
Thickness=[-0.2:0.01:0.2]*1E-6+2.52.*1E-6;
Load_Previous_Result=0;
Position_Error=0*1E-6;%[0.019:0.001:0.024]*1E-6;%[0.02:0.0003:0.023]*1E-6;
%%
n_should=1.7;
delta_n=0.2;
Wavelength_Center=600;

pixel_1=700;               % DC cutoff
%pixel_2=1800;%2000               % the end of 2
%pixel_3=1120;               %for saperation

Wavelength_Considered_Min=500;          %nm
Wavelength_Considered_Max=700;

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


cd('D:\120328\Glass_inter\');
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

r_BK7=((n_bk7-1)./(n_bk7+1));
t_AN100=(2*(n_bk7)./(n_bk7+1));
%Spectrum_Original=Spectrum_Reference./(r_BK7).^2;

%% Data Loading

Filter_inter=0;
cd('D:\120328\Glass_inter\');
for j=0:99
    Filter_inter=Filter_inter+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Filter_inter=Filter_inter/100;

Filter_sam=0;
cd('D:\120328\Glass_sam\');
for j=0:99
    Filter_sam=Filter_sam+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Filter_sam=Filter_sam/100;


ref=0;
cd('D:\120328\ref\');
for j=0:99
    ref=ref+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
ref=ref/100;


Glass_inter=0;
cd('D:\120328\Glass_inter\');
for j=0:99
    Glass_inter=Glass_inter+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Glass_inter=Glass_inter/100;


Glass_sam=0;
cd('D:\120328\Glass_sam\');
for j=0:99
    Glass_sam=Glass_sam+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Glass_sam=Glass_sam/100;

Spectrum_Sample_Wavelength=Filter_inter(:,2)-Filter_sam(:,2)-ref(:,2);
Spectrum_Reference_Wavelength=Glass_inter(:,2)-Glass_sam(:,2)-ref(:,2);

Spectrum_Sample_Wavelength=Spectrum_Sample_Wavelength-Spectrum_Sample_Wavelength(1);

Spectrum_Reference_Wavelength=Spectrum_Reference_Wavelength-Spectrum_Reference_Wavelength(1);

plot(Wavelength,Spectrum_Sample_Wavelength,Wavelength,Spectrum_Reference_Wavelength);

plot(Wavelength,Spectrum_Sample_Wavelength,Wavelength,Spectrum_Reference_Wavelength);

%% Spectrum generation


Gauss_Window=2;

Spectrum_Sample_Frequency=(Spectrum_Sample_Wavelength.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectrum_Sample=interp1(Frequency_Old,Spectrum_Sample_Frequency,Frequency); 
Spectrum_Sample(isnan(Spectrum_Sample))=0;
Spectrum_Sample(Frequency<Min_Frequency)=0;
%Signal_Sample=fft(Spectrum_Sample).*gaussmf(Position_micron,[10 27.4]);  
%Signal_Sample(round(length(Signal_Sample)/2)+1:end)=0;
%Spectrum_Sample=(ifft(Signal_Sample));
%Spectrum_Sample=2*Spectrum_Sample(1:N_f);

Spectrum_Reference_Frequency=(Spectrum_Reference_Wavelength.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency); 
Spectrum_Reference(isnan(Spectrum_Reference))=0;
Spectrum_Reference(Frequency<Min_Frequency)=0;
Spectrum_Reference(N_f+1:N_t)=0;

%Signal_Reference=fft(Spectrum_Reference).*gaussmf(Position_micron,[Gauss_Window Spatial_Range(j)]);  
%Signal_Reference(round(length(Signal_Reference)/2)+1:end)=0;
%Spectrum_Reference(:,j)=(ifft(Signal_Reference));
%Spectrum_Reference(:,j)=2*Spectrum_Reference(1:N_f);

plot(Frequency,Spectrum_Sample);


Spatial_Range=10:0.5:50;

for j=1:length(Spatial_Range)
    Spectrum_Sample(N_f+1:N_t)=0;
    Signal_Sample=fft(Spectrum_Sample).*gaussmf(Position_micron,[Gauss_Window Spatial_Range(j)]);  
    Signal_Sample(round(length(Signal_Sample)/2)+1:end)=0;
    Spectrum_Sample_Part=(ifft(Signal_Sample));
    Time_Spectrum_Sample(:,j)=2*Spectrum_Sample_Part(1:N_f);
    
end
imagesc(abs(Time_Spectrum_Sample),'xdata',Spatial_Range,'ydata',Frequency);