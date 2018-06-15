clear all

%%% To find the phase from min voltage to max voltage ONLY

Path='D:\Users\TuanShu\121219_PZT Calibration_renamed\';

Parameter='Frequency';

%Frequency=[0.001 0.002];
Frequency=[0.005];
Sampling_rate=2500; %Hz
smooth_window=10000;
Range=8;        %V
q=1;
min_voltage=(10-Range)/2;
max_voltage=min_voltage+Range;
dvoltage_min=0.001;
dvoltage_max=0.1;
Wavelength_Laser=0.6328;     %micron

dt=1/Sampling_rate;
cd(Path);
clf
axes('FontSize',10);
hold all
for p=1:length(Frequency)
    Data=importdata(sprintf('121219_f%g_V%d.txt',Frequency(p),Range(q)));                   %Asuume: line 1: time, line 2: interference signal, line 3: voltage 121219_f0.001_V2to8.txt©Ç©Ç


    Voltage=Data(:,1)*10;

    Signal=Data(:,2);

    Voltage_min_index_temp=find(Voltage<(min_voltage+dvoltage_min),1,'first');
    Voltage_next_max_index=Voltage_min_index_temp+find(Voltage((Voltage_min_index_temp+1):end)>(max_voltage-dvoltage_max),1,'first');
    [value Voltage_min_index]=min(Voltage((Voltage_min_index_temp+1):Voltage_next_max_index));
    Voltage_min_index=Voltage_min_index_temp+Voltage_min_index;
    
    Voltage=Voltage(Voltage_min_index:Voltage_next_max_index);

    Signal=Signal(Voltage_min_index:Voltage_next_max_index);


    Time_index=1:length(Signal);

    Spectrum=ifft(Signal);

    Lower_Band=20;

    Upper_Band=fix(length(Spectrum)/2);

    Spectrum(1:Lower_Band)=0;
    Spectrum(Upper_Band:end)=0;

    Signal_New=fft(Spectrum,[],1);
    Phase_original=angle(Signal_New);
    Phase=unwrap(Phase_original);

    Position=-1*Wavelength_Laser*Phase/(2*pi)/2;
    Position=Position-Position(find(Voltage>5,1,'first'));
    
    Position=smooth(Position,smooth_window);
    Voltage=smooth(Voltage,smooth_window);
    
    %Velocity=diff(Position)/dt;
    %Velocity(length(Position))=Velocity(length(Position)-1);
    %plot1=subplot(2,1,1);
    %plot(Voltage,Position);
    %xlabel('Voltage (V)');
    %ylabel('PZT Position (micron)');
    %plot2=subplot(2,1,2);

    Position_relative=Position-(Voltage-5)*10;
    plot(Voltage,Position_relative);

    disp(p);
end
hold off
xlabel('Voltage (V)','FontSize',12);
ylabel('PZT Position (micron)','FontSize',12);
%legend('Frequency = 0.001 Hz','Frequency = 0.002 Hz','Frequency = 0.005 Hz','Frequency = 0.01 Hz','Frequency = 0.02 Hz','Frequency = 0.05 Hz');
%%% Output: voltage vs. position
%%% OCT data: must know the voltage for each frame (put DAQ read in the
%%% same loop of frame grabber?)
ylabel('Signal','FontSize',12);
