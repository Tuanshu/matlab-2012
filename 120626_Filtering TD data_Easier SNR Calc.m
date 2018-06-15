clear all

%% Setting
SR=40000;
Stage_Speed=2;  %mm/s

Stage_Speed_MTS=(Stage_Speed*1E-3);

C=3E8;

TH_first_peak=0.06;
TH_second_peak=0.01;
Delay_first_peak=80000; %the smaller one, pixel
Min_thickness=12000; %pixel

Max_thickness=16000; %pixel

cd('D:\120716_CCT_4V_PD gain 60dB\');

qqq=1;
xxx=1;
yyy=1;
for j=2:27

%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata(sprintf('Data%d.txt',j));

N_t=length(Data);
N_f=N_t;

Time_Stage=1/SR:1/SR:(1/SR)*N_t;

Time=Time_Stage*Stage_Speed_MTS/C;
Position_micron=Time*C*1E6;

dTime=(Time(2)-Time(1))*2;      %*2 becuase of round trip

Frequency_Max=1/dTime;

Frequency=Frequency_Max/N_f:Frequency_Max/N_f:Frequency_Max;

Wavelength_micron=(C./Frequency)*1E6;

Spectrum=fft(Data,N_f);


Window1=(gaussmf(Frequency,[0.2E14 1.5E14]));
Window1(Frequency>1.5E14)=1;
Window2=(gaussmf(Frequency,[0.2E14 5.5E14]));
Window2(Frequency<5.5E14)=1;
Window=(Window1.*Window2)';
Spectrum=Window.*Spectrum;
Spectrum((round(length(Spectrum)/2)+1):end)=0;

plot(Frequency,Spectrum);

Data_New=ifft(Spectrum);
Data_New=Data_New(1:N_t);

Max_Wavelength_micron=1;

Min_Wavelength_micron=0.6;

Max_Wavelength_micron_index=find(Wavelength_micron<Max_Wavelength_micron,1,'first');

Min_Wavelength_micron_index=find(Wavelength_micron<Min_Wavelength_micron,1,'first');

Data_New(1:40)=Data_New(40);
Data_New((length(Data_New)-40):end)=Data_New((length(Data_New)-40));

clear left_index_array max_index_array max_index_array_second max_value_array_raw max_value_array max_index_array_raw max_index_array_second_raw max_value_array_second_raw max_value_array max_index_array_second max_index_array max_index_array_second_even max_index_array_even max_value_array_even max_index_array_second_odd max_index_array_odd max_value_array_odd
index_array_first_peak=find(Data_New>TH_first_peak);
cont=1;
current_index=index_array_first_peak(1);
p=1;
    while(cont==1)
        left_index_array(p)=current_index;
        p=p+1;
        current_index=index_array_first_peak(find(index_array_first_peak>(current_index+Delay_first_peak),1,'first'));
        if current_index<length(Data_New)
            cont=1;
        else
            cont=0;
        end
    end
left_index_array(length(left_index_array)+1)=length(Data_New);
cont=1;
q=1;
    while(cont==1)
        [max_value_array_raw(q) max_index]=max(abs(Data_New((left_index_array(q)+1):left_index_array(q+1))));
        max_index_array_raw(q)=max_index+left_index_array(q);
        q=q+1;
        if q<length(left_index_array)
            cont=1;
        else
            cont=0;
        end
    end
    x=1;
    y=1;
    ww=1;
    for w=1:length(max_index_array_raw)
        Data_Temp=Data_New;
        if (max_index_array_raw(w)-Max_thickness)>0
            Data_Temp(1:max_index_array_raw(w)-Max_thickness)=0;
        end
        if max_index_array_raw(w)+Max_thickness<length(Data_Temp)
            Data_Temp((max_index_array_raw(w)+Max_thickness):end)=0;
        end
        if ((max_index_array_raw(w)-Min_thickness)>0)&&((max_index_array_raw(w)+Min_thickness)<length(Data_Temp))
            Data_Temp((max_index_array_raw(w)-Min_thickness):(max_index_array_raw(w)+Min_thickness))=0;
        elseif ((max_index_array_raw(w)-Min_thickness)>0)
            Data_Temp((max_index_array_raw(w)-Min_thickness):end)=0;
        elseif ((max_index_array_raw(w)+Min_thickness)<length(Data_Temp))
            Data_Temp(1:(max_index_array_raw(w)+Min_thickness))=0;
        end
        [max_value_array_second_raw(w)  max_index_array_second_raw(w)]=max(abs(Data_Temp));
        if max_value_array_second_raw(w)>TH_second_peak
            max_value_array(ww)=max_value_array_raw(w);
            max_index_array(ww)=max_index_array_raw(w);
            max_index_array_second(ww)=max_index_array_second_raw(w);

            if max_index_array_second(ww)>max_index_array_raw(w)
                max_index_array_second_even(x)=max_index_array_second_raw(w);
                max_index_array_even(x)=max_index_array_raw(w);
                max_value_array_even(x)=max_value_array_raw(w);
                x=x+1;
            elseif max_index_array_second(ww)<max_index_array_raw(w)
                max_index_array_second_odd(y)=max_index_array_second_raw(w);
                max_index_array_odd(y)=max_index_array_raw(w);
                max_value_array_odd(y)=max_value_array_raw(w);
                y=y+1;
            end
            ww=ww+1;
        end
    end
    thickness_array(qqq:(qqq+length(max_index_array)-1))=abs(Position_micron(max_index_array_second)-Position_micron(max_index_array));
    max_array(qqq:(qqq+length(max_index_array)-1))=max_value_array;
    qqq=qqq+length(max_index_array);
    thickness_array_even(xxx:(xxx+length(max_index_array_even)-1))=abs(Position_micron(max_index_array_second_even)-Position_micron(max_index_array_even));
    max_array_even(xxx:(xxx+length(max_index_array_even)-1))=max_value_array_even;
    xxx=xxx+length(max_index_array_even);
    
    thickness_array_odd(yyy:(yyy+length(max_index_array_odd)-1))=abs(Position_micron(max_index_array_second_odd)-Position_micron(max_index_array_odd));
    max_array_odd(yyy:(yyy+length(max_index_array_odd)-1))=max_value_array_odd;
    yyy=yyy+length(max_index_array_odd);
end

xlabel('Wavelength (micron)');
ylabel('Interference Power (a.u.)');

plot(Wavelength_micron(Max_Wavelength_micron_index:Min_Wavelength_micron_index),Spectrum(Max_Wavelength_micron_index:Min_Wavelength_micron_index));

[C ind]=max(abs(Data_New));
plot(Position_micron-Position_micron(ind),Data_New/max(abs(Data_New)),Position_micron-Position_micron(ind),abs(Data_New)/max(abs(Data_New)));
xlabel('OPD (micron)');
ylabel('Interference Signal');
plot(thickness_array);
xlabel('Nth Scan');
ylabel('Measured Thickness (OPD, micron)');
Standard_Deviation=std(thickness_array);
Standard_Deviation_even=std(thickness_array_even);
Standard_Deviation_odd=std(thickness_array_odd);


max_array_fitting=polyfit([1:length(max_array)],max_array,1);
max_array_even_fitting=polyfit([1:length(max_array_even)],max_array_even,1);
max_array_odd_fitting=polyfit([1:length(max_array_odd)],max_array_odd,1);
%max_array_even_fitting=polyfit([1:length(max_array_even)],max_array_even,1);
%max_array_odd_fitting=polyfit([1:length(max_array_odd)],max_array_odd,1);

max_array_filtered=max_array-(max_array_fitting(1).*[1:length(max_array)]+max_array_fitting(2));

max_array_even_filtered=max_array_even-(max_array_even_fitting(1).*[1:length(max_array_even)]+max_array_even_fitting(2));
max_array_odd_filtered=max_array_odd-(max_array_odd_fitting(1).*[1:length(max_array_odd)]+max_array_odd_fitting(2));

%max_array_even_filtered=max_array_even-(max_array_even_fitting(1).*[1:length(max_array_even)]+max_array_even_fitting(2));

%max_array_odd_filtered=max_array_odd-(max_array_odd_fitting(1).*[1:length(max_array_odd)]+max_array_odd_fitting(2));

plot(max_array_filtered);
xlabel('Nth Scan');
ylabel('Interference Signal Power Deviation (a.u.)');

plot(thickness_array);
xlabel('Nth Scan');
ylabel('Measured Thickness (OPD, micron)');


%SD=std(abs(Data_New(1.26E5:1.36E5)));

SD=std(max_array_odd_filtered);
SNR=10*log10((max(abs(max_array))^2)/(SD^2));
%SNR_even=10*log10((max(abs(max_array_even))^2)/(SD^2));
%SNR_odd=10*log10((max(abs(max_array_odd))^2)/(SD^2));
%plot(Wavelength_micron,Spectrum);
%xlabel('Wavelength (micron)');
%label('Interference Spectrum');
%cd('D:\120222\');

%dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');
%plot(Position_micron,Signal_Bscan_Envelope);