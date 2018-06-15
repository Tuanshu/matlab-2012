clear all

%%% To find the phase from min voltage to max voltage ONLY

Path='D:\Users\TuanShu\121219_PZT Calibration_renamed\';

Frequency=[0.001 0.002 0.005];

Wavelength_Laser=0.6328;     %micron

Sampling_rate=1;             %Hz 其實好像用不到

cd(Path);
newplot
hold all
for p=1:length(Frequency)
    Data=importdata(sprintf('121219_f%g_V2.txt',Frequency(p)));                   %Asuume: line 1: time, line 2: interference signal, line 3: voltage 121219_f0.001_V2to8.txt怪怪


    Voltage=Data(:,1)*10;

    Signal=Data(:,2);

    [Voltage_min_value Voltage_min_index]=min(Voltage);

    [Voltage_next_max_value Voltage_next_max_index]=max(Voltage((Voltage_min_index+1):end));

    Voltage_next_max_index=Voltage_min_index+Voltage_next_max_index;

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
    plot(Voltage,Position);
    disp(p);
end
hold off
%%% Output: voltage vs. position
%%% OCT data: must know the voltage for each frame (put DAQ read in the
%%% same loop of frame grabber?)
