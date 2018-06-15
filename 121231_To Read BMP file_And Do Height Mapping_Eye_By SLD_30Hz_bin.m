clear all

%Too bad axial resolution
fitting=0;
new_bit_depth=8;
cd('D:\Users\TuanShu\121231_Green Color Filter\3\Sam\');
set_peak_index=500;
Averaging_Factor=10;

Binning_Factor=1;

Stage_speed=1.5;  %micron/sec
Sampling_rate=28;   %Hz

Frame_Axial_Spacing=Stage_speed/Sampling_rate*5/3;  %micron

Objective_Focal_Length=4/0.85;   %cm
Porjection_Lens_Focal_Length=50;    %cm
Pixel_Size=7.4;       %micron

%Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;
Lateral_Spacing=600/(538-330);

image_index=0:30000;
ROI=[1 61;1 80];      %up, down, left, right
%ROI=[43 342;363 662];      %up, down, left, right
%ROI=[1 384;257 768];      %up, down, left, right
%%ROI=[151 234;457 568];      %up, down, left, right
Image_Stack(1:ceil(length(image_index)/Averaging_Factor),1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
Low_Pass=100;   %pixel
High_Pass=400;
Image_Temp_2(1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
for p=1:length(image_index)
    Image_Temp=imread(sprintf('%06i.tif',image_index(p)));   
    
    [m,n]=size(Image_Temp); %M is the original matrix

    Image_Temp=sum(reshape(Image_Temp,Binning_Factor,[]),1);
    Image_Temp=reshape(Image_Temp,m/Binning_Factor,[]).';

    Image_Temp=sum(reshape(Image_Temp,Binning_Factor,[]) ,1);
    Image_Temp=reshape(Image_Temp,n/Binning_Factor,[]).';
    if rem(p,Averaging_Factor)==0
        Image_Temp_2=Image_Temp(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2))+Image_Temp_2;
        Image_Stack(ceil(p/Averaging_Factor),:,:)=Image_Temp_2;
        Image_Temp_2=0;
    else
        Image_Temp_2=Image_Temp(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2))+Image_Temp_2;
    end
    %Image_Stack(ceil(p/Averaging_Factor),:,:)=Image_Temp;
    %Image=imread(sprintf('1.jpg',image_index(p)));    
    disp(p);
end
Image_Stack(size(Image_Stack,1),:,:)=Image_Stack(size(Image_Stack,1)-1,:,:);

Start_DC=Image_Stack(1,:,:);
End_DC=Image_Stack(end,:,:);
DC_volume(1:size(Image_Stack,1),1:size(Image_Stack,2),1:size(Image_Stack,3))=0;
for q=1:size(Image_Stack,1)
    DC_volume(q,:,:)=Start_DC+(End_DC-Start_DC)*q/size(Image_Stack,1);
end
Image_Stack=Image_Stack-DC_volume;
[Max_value Max_index]=max(Image_Stack(:,30,40));

Image_Stack=circshift(Image_Stack,set_peak_index-Max_index);

FFT_Stack=fft(Image_Stack,[],1);

FFT_Stack(1:Low_Pass,:,:)=0;
FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
FFT_Stack(High_Pass:end,:,:)=0;

%clear Image_Stack Image Image_Temp
%FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
Image_Stack_New=(ifft(FFT_Stack,[],1));
%plot(real(FFT_Stack(:,500,400)));
New_Max=max(Image_Stack_New(:));

Image_Stack_New_Normalized=Image_Stack_New./max(Image_Stack_New(:));



%

Sam=Image_Stack_New_Normalized(:,30,40);
Ref=Image_Stack_New_Normalized(:,30,40);
Sam_real=real(Sam);
Ref_real=real(Ref);
plot(1:length(Sam),Sam,1:length(Ref),Ref);
dlmwrite('Sam.txt',Sam,'delimiter','\t','newline','pc');

%dlmwrite('Ref.txt',Ref,'delimiter','\t','newline','pc');


dlmwrite('Sam_real.txt',Sam_real,'delimiter','\t','newline','pc');

%dlmwrite('Ref_real.txt',Ref_real,'delimiter','\t','newline','pc');
%for rr=1:size(Image_Stack_New_Normalized,3)
%    Temp_Bscan(:,:)=Image_Stack_New_Normalized(:,:,rr);
%    imwrite(Temp_Bscan,sprintf('Bscan_%i.png',rr),'Bitdepth',new_bit_depth);
%    disp(rr);
%end

for q=1:size(Image_Stack,1)
    Temp_Slides(:,:)=Image_Stack_New_Normalized(q,:,:);
    imwrite(Temp_Slides,sprintf('Enface_Binned_by_%i_%i.png',Binning_Factor,image_index(q)),'Bitdepth',new_bit_depth);
    disp(q);
end
%Temp_Sides_Read=imread(sprintf('%i_New.png',image_index(q)));








