clear all;
close all;
clc;

ref_sub=1;

total_OPD=156/4;        %micron ask ¾G¤D¹Å
axial_resolution=1.5;   %micron

lateral_step=1;    %micron
length_lateral=2250;

window_ref=10;

data_Bscan=importdata('D:\120111\3mm 1 (not enough range)(the 0.4 one).txt');
data_Bscan=data_Bscan(:,1:round(length_lateral/lateral_step));


%Is=sum(Is_o,2)/size(Is_o,2);

axial_position=[0:total_OPD/(size(data_Bscan,1)-1):total_OPD]';  

lateral_position=[0:lateral_step:lateral_step*(size(data_Bscan,2)-1)]';  

short_pass=2.5;   %in terms of pixel

long_pass=750;

%% in the frequency domain: 2/2250 for each pixel

long_pass_f=round(1/long_pass/(1/length_lateral));
short_pass_f=round(1/short_pass/(1/length_lateral));

%% filtering and manual hilbert

data_Bscan_f=fft(data_Bscan,[],1);
start_index_of_spectrum=50;

end_index_of_spectrum=250;

data_Bscan_f(1:start_index_of_spectrum,:)=0;
data_Bscan_f(end_index_of_spectrum:end,:)=0;
data_Bscan_env=abs(ifft(data_Bscan_f,[],1));


plot(axial_position,data_Bscan_env(:,1));

data_Bscan_env(1:50,:)=0;
data_Bscan_env((size(data_Bscan_env,1)-49:end),:)=0;

%% Finding the inerface
[value_max index_max]=max(data_Bscan_env,[],1);
profile=axial_position(index_max);
profile=profile';
%% Another mehod
for j=1:size(data_Bscan_env,2)
    index_left=find(data_Bscan_env(:,j)>value_max(j)/2,1,'first');
    index_right=find(data_Bscan_env(:,j)>value_max(j)/2,1,'last');
    FWHM(j)=axial_position(index_right)-axial_position(index_left);
    profile_another(j)=(axial_position(index_right)+axial_position(index_left))/2;
end


%% Ref sub

if ref_sub==1
    ref_Bscan=importdata('D:\120111\3mm ref.txt');  
    ref_Bscan=ref_Bscan(:,1:round(length_lateral/lateral_step)); 
    
    
    ref_Bscan_f=fft(ref_Bscan,[],1);
    
    ref_Bscan_f(1:start_index_of_spectrum,:)=0;
    ref_Bscan_f(end_index_of_spectrum:end,:)=0;
    ref_Bscan_env=abs(ifft(ref_Bscan_f,[],1));
    
    ref_Bscan_env(1:50,:)=0;
    ref_Bscan_env(((size(ref_Bscan_env,1)-49):end),:)=0;
    
    [value_max_ref index_max_ref]=max(ref_Bscan_env,[],1);
    profile_ref=axial_position(index_max_ref);
    
    for j=1:size(ref_Bscan_env,2)
        index_left_ref=find(ref_Bscan_env(:,j)>value_max_ref(j)/2,1,'first');
        index_right_ref=find(ref_Bscan_env(:,j)>value_max_ref(j)/2,1,'last');
        FWHM_ref(j)=axial_position(index_right_ref)-axial_position(index_left_ref);
        profile_ref_another(j)=(axial_position(index_right_ref)+axial_position(index_left_ref))/2;
    end
end

profile=profile_another;

if ref_sub==1
    plot(lateral_position,profile,lateral_position,profile_another);
    profile_ref=profile_ref_another;
    profile_ref_sub=profile-profile_ref;
else
    profile_ref_sub=profile;
end

fitting=polyfit([1:length(profile_ref_sub)],profile_ref_sub,1);

fitted_curve=fitting(1).*[1:length(profile_ref_sub)]+fitting(2);

plot(lateral_position,fitted_curve,lateral_position,profile_ref_sub);

profile_sub=profile_ref_sub-fitted_curve;

profile_f=fft(profile_sub);

profile_f_filtered=profile_f;

profile_f_filtered(1:long_pass_f)=0;

profile_f_filtered((length(profile_f_filtered)-long_pass_f+1):end)=0;

profile_f_filtered(short_pass_f:(length(profile_f_filtered)-short_pass_f+1))=0;

profile_new=real(ifft(profile_f_filtered));


plot(lateral_position,profile_sub,lateral_position,profile_new);

profile_new=profile_new-mean(profile_new);

Ra=sum(abs(profile_new))/size(profile_new,2);
Rq=(sum((profile_new).^2)/size(profile_new,2)).^0.5;
Ry=max(profile_new)-min(profile_new);




imagesc(10*log10(data_Bscan_env(end:-1:1,:)),'xdata',lateral_position,'ydata',axial_position);

plot(lateral_position,profile_sub,lateral_position,profile_new);
%% Ra, Ry, Rz, Rq calculation

%dlmwrite('D:\position.txt',position,'delimiter','\t','newline','pc');
