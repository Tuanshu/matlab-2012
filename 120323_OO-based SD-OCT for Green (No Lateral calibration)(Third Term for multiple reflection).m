%% Options

clear all
index_calculation=1;

array=1:999;
Gauss_Window=10;
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
Thickness=([-0.1:0.01:0.1]+2.6).*1E-6;
Load_Previous_Result=0;
Position_Error=[0]*1E-6;%[0.019:0.001:0.024]*1E-6;%[0.02:0.0003:0.023]*1E-6;
%%
n_should=1.6;
delta_n=0.2;
Wavelength_Center=600;

pixel_1=700;               % DC cutoff
%pixel_2=1800;%2000               % the end of 2
%pixel_3=1120;               %for saperation

Wavelength_Considered_Min=550;          %nm
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


cd('D:\120313\Filter_inter\');
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
cd('D:\120313\Filter_inter\');
for j=0:99
    Filter_inter=Filter_inter+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Filter_inter=Filter_inter/100;

Filter_sam=0;
cd('D:\120313\Filter_sam\');
for j=0:99
    Filter_sam=Filter_sam+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Filter_sam=Filter_sam/100;


ref=0;
cd('D:\120313\ref\');
for j=0:99
    ref=ref+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
ref=ref/100;


Glass_inter=0;
cd('D:\120313\Glass_inter\');
for j=0:99
    Glass_inter=Glass_inter+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Glass_inter=Glass_inter/100;


Glass_sam=0;
cd('D:\120313\Glass_sam\');
for j=0:99
    Glass_sam=Glass_sam+importdata(sprintf('D%i.txt',j));      % Data_2: the glass data 1
end
Glass_sam=Glass_sam/100;

Spectrum_Sample_Wavelength=Filter_inter(:,2)-Filter_sam(:,2)-ref(:,2);
%Spectrum_Sample_Wavelength=ref(:,2);
Spectrum_Reference_Wavelength=Glass_inter(:,2)-Glass_sam(:,2)-ref(:,2);

Spectrum_Sample_Wavelength=Spectrum_Sample_Wavelength-Spectrum_Sample_Wavelength(1);

Spectrum_Reference_Wavelength=Spectrum_Reference_Wavelength-Spectrum_Reference_Wavelength(1);

plot(Wavelength,Spectrum_Sample_Wavelength,Wavelength,Spectrum_Reference_Wavelength);

%% Spectrum generation

Gauss_Window=10;

Spectrum_Sample_Frequency=(Spectrum_Sample_Wavelength.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectrum_Sample=interp1(Frequency_Old,Spectrum_Sample_Frequency,Frequency); 
Spectrum_Sample(isnan(Spectrum_Sample))=0;
Spectrum_Sample(Frequency<Min_Frequency)=0;
Spectrum_Sample(N_f+1:N_t)=0;
Signal_Sample=fft(Spectrum_Sample);%.*gaussmf(Position_micron,[10 27.4]);  
Signal_Sample(round(length(Signal_Sample)/2)+1:end)=0;
Spectrum_Sample=(ifft(Signal_Sample));
Spectrum_Sample=2*Spectrum_Sample(1:N_f);

Spectrum_Reference_Frequency=(Spectrum_Reference_Wavelength.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency); 
Spectrum_Reference(isnan(Spectrum_Reference))=0;
Spectrum_Reference(Frequency<Min_Frequency)=0;
Spectrum_Reference(N_f+1:N_t)=0;
Signal_Reference=fft(Spectrum_Reference);%.*gaussmf(Position_micron,[Gauss_Window 27.4]);  
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

plot(Position_micron,Signal_Sample,Position_micron,abs(Signal_Sample));
plot(Position_micron,Signal_Sample,Position_micron,abs(Signal_Sample));

xlabel('Optical path Difference (micron)');
ylabel('Interference Signal');


plot(Wavelength_micron(Wavelength_micron<0.80),real(Spectrum_Devided(Wavelength_micron<0.80,:)),Wavelength_micron(Wavelength_micron<0.80),abs(Spectrum_Devided(Wavelength_micron<0.80,:)));

xlabel('Wavelength (micron)');
ylabel('A');

%% to generate the phase spectrum

%% Start the n k fitting    

if index_calculation==1   
    
    n_final(1:length(Frequency),1:length(Thickness))=1.5;
    k_final(1:length(Frequency),1:length(Thickness))=0;
    T_check(1:length(Frequency),1:length(Thickness))=0;  
    Spectrum_Check(1:length(Frequency),1:length(Thickness))=0;
 
    for q=1:length(Thickness)
        Position_Error_Now=Position_Error;
        Thickness_Now=Thickness(q);
        value_temp=10000000000;
        n(1:length(Frequency),1)=1.5;
        k(1:length(Frequency),1)=0;

        %T_max_check(1:length(Spectrum),1:length(n_should))=0;

        value_total_1(1:Number_of_Loop)=0;
        value_total_2(1:Number_of_Loop)=0;

        n_o=(n_should-delta_n):0.0001:(n_should+delta_n);
        k_o=-0.1:0.0002:0.3;
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
    
        n_should_index=find(n_o-n_should>0,1,'first');
    
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
            
                r1=(n1(index_now)-(n_temp+i*k(index_now)))./(n_temp+n1(index_now)+i*k(index_now));
                r1_r=((n_temp+i*k(index_now))-n1(index_now))./(n_temp+n1(index_now)+i*k(index_now));    
                t1=2*(n1(index_now))./(n_temp+n1(index_now)+i*k(index_now));
                t1_r=2*(n_temp+i*k(index_now))./(n_temp+n1(index_now)+i*k(index_now));
                t2=2*(n_temp+i*k(index_now))./(n_temp+n2+i*k(index_now));
                t2_r=2*(n2)./(n_temp+n2+i*k(index_now));
                r2=((n_temp+i*k(index_now))-n2)./((n_temp+i*k(index_now))+n2);   
                d=exp(i*2*pi.*Frequency(index_now)/C.*(n_temp+i*k(index_now)).*Thickness_Now);   %注意! -1*n!
                d0=exp(i*2*pi.*Frequency(index_now)/C.*(Position_Error_Now));   %注意! -1*n!
                Spectrum_Upper_Temp=((r1)./r_BK7(index_now));                       %神說是cosine
                Spectrum_Lower_Temp=((t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now))./(1-r1_r.*(d.^2).*r2);                         %神說是cosine                
                %Spectrum_Multi_Temp=((t1.*t1_r.*r2.*(d.^2).*(d.^2).*r1_r.*r2)./r_BK7(index_now));                         %神說是cosine
                Spectrum_Temp=Spectrum_Upper_Temp+Spectrum_Lower_Temp;%+Spectrum_Multi_Temp;
                
                Merit=(((Spectrum_Temp)-(Spectrum_Devided(index_now).*(d0.^2))).^2); %(abs(abs(Spectrum_Temp)-abs(Spectrum(index_now))).^2)+

                if range_specified == 1
                    if (index_now == Frequency_Center_Index) || (index_now == (Frequency_Center_Index -1))
                        range_upper=min(n_should_index+delta_n_number,length(Merit));
                        range_lower=max(n_should_index-delta_n_number,1);
                    elseif (index_now > Frequency_Center_Index) || (index_now < (Frequency_Center_Index -1))
                        range_upper=min(index_old+round(delta_n_number/20),length(Merit));
                        range_lower=max(index_old-round(delta_n_number/20),1);
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
                d=exp(i*2*pi.*Frequency(index_now)/C.*(n(index_now)+i*k_temp).*Thickness_Now);   %注意! -1*n!
                %d0=exp(i*2*pi.*Frequency(index_now)/C.*(Distance_fit(j)-Thickness));   %注意! -1*n!
                d0=exp(i*2*pi.*Frequency(index_now)/C.*(Position_Error_Now));   %注意! -1*n!
                Spectrum_Upper_Temp=((r1)./r_BK7(index_now));                       %神說是cosine
                Spectrum_Lower_Temp=((t1.*t1_r.*r2.*(d.^2))./r_BK7(index_now))./(1-r1_r.*(d.^2).*r2);                         %神說是cosine
                %Spectrum_Multi_Temp=((t1.*t1_r.*r2.*(d.^2).*(d.^2).*r1_r.*r2)./r_BK7(index_now));                         %神說是cosine
                Spectrum_Temp=Spectrum_Upper_Temp+Spectrum_Lower_Temp;%+Spectrum_Multi_Temp;
                Merit=(((Spectrum_Temp)-(Spectrum_Devided(index_now).*(d0.^2))).^2); %(abs(abs(Spectrum_Temp)-abs(Spectrum(index_now))).^2)+


                [value index]=min(Merit);
                value_total_2(p)=value_total_2(p)+abs(value);
                k(index_now)=k_temp(index);
            end

        end


        n_final(:,q)=n;
        k_final(:,q)=k;    
    end

    
end    

%%
for p=1:size(k_final,2)
    for q=1:size(k_final,1)
        if q>1
            if k_final(q,p)>0.08
                k_final(q,p)=k_final(q-1,p);
            end
        end
    end
end
                

%% Fitting

if index_calculation == 1
%% Checking
T_abs_check(1:N_f,1:length(Thickness))=0;
for j=1:length(Thickness)
    if index_calculation == 1        
    t_sub=2.*n2./(n2+n1);
    r1_check=(n1-(n_final(:,j)+i*k_final(:,j)))./(n_final(:,j)+n1+i*k_final(:,j));
    r1_r_check=((n_final(:,j)+i*k_final(:,j))-n1)./(n_final(:,j)+n1+i*k_final(:,j));    
    t1_check=2*(n1)./(n_final(:,j)+n1+i*k_final(:,j));
    t1_r_check=2*(n_final(:,j)+i*k_final(:,j))./(n_final(:,j)+n1+i*k_final(:,j));
    t2_check=2*(n_final(:,j)+i*k_final(:,j))./(n_final(:,j)+n2+i*k_final(:,j));
    t2_r_check=2*(n2)./(n_final(:,j)+n2+i*k_final(:,j));
    r2_check=((n_final(:,j)+i*k_final(:,j))-n2)./((n_final(:,j)+i*k_final(:,j))+n2);   
    d_n_check=exp(i*2*pi.*Frequency/C.*(n_final(:,j)).*Thickness_Now);   %注意! -1*n!
    d_check=exp(i*2*pi.*Frequency/C.*(n_final(:,j)+i*k_final(:,j)).*Thickness_Now);   %注意! -1*n!
    %d0_check=exp(i*2*pi.*Frequency/C.*(Distance_fit(j)-Thickness));   %注意! -1*n!

    Spectrum_Upper_Temp_Check=((r1_check)./r_BK7);                       %神說是cosine
    Spectrum_Lower_Temp_Check=((t1_check.*t1_r_check.*r2_check.*(d_check.^2))./r_BK7); 
    T_abs_check(:,j)=abs(t_sub.*t1_check.*t2_check.*(d_check)./(1-r2_check.*r1_r_check.*(d_check).^2)).^2;
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
end
%%

plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),n_final(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1));

xlabel('Wavelength (micron)');
ylabel('n');


plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),k_final(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1));

xlabel('Wavelength (micron)');
ylabel('k');
%%
plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,8));

%plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1)/max(T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,5)));

xlabel('Wavelength (micron)');
ylabel('T (percent)');


%plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1)/max(T_abs_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1)),Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),QQ(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index)/max(QQ(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index)));

%xlabel('Wavelength (micron)');
%ylabel('T (percent)');


%plot(Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),d_check(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,1));


dlmwrite('QQ.txt',QQ,'delimiter','\t','newline','pc');
