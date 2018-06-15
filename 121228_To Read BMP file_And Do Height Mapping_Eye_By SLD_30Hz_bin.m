clear all

%Too bad axial resolution
fitting=0;
new_bit_depth=8;
cd('D:\Users\TuanShu\121227_GlassLower_50_50_f0.005\');

Binning_Factor=1;

Stage_speed=1.5;  %micron/sec
Sampling_rate=28;   %Hz

Frame_Axial_Spacing=Stage_speed/Sampling_rate*5/3;  %micron

Objective_Focal_Length=4/0.85;   %cm
Porjection_Lens_Focal_Length=50;    %cm
Pixel_Size=7.4;       %micron

%Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;
Lateral_Spacing=600/(538-330);

image_index=0:65000;
ROI=[1 50;1 50];      %up, down, left, right
%ROI=[43 342;363 662];      %up, down, left, right
%ROI=[1 384;257 768];      %up, down, left, right
%%ROI=[151 234;457 568];      %up, down, left, right
Image_Stack(length(image_index),1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
Low_Pass=30;   %pixel
High_Pass=200;
for p=1:length(image_index)
    Image_Temp=imread(sprintf('%06i.tif',image_index(p)));   
    
    [m,n]=size(Image_Temp); %M is the original matrix

    Image_Temp=sum(reshape(Image_Temp,Binning_Factor,[]),1);
    Image_Temp=reshape(Image_Temp,m/Binning_Factor,[]).';

    Image_Temp=sum( reshape(Image_Temp,Binning_Factor,[]) ,1);
    Image_Temp=reshape(Image_Temp,n/Binning_Factor,[]).';
    
    %Image=imread(sprintf('1.jpg',image_index(p)));    
    Image_Stack(p,:,:)=Image_Temp(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2));
    disp(p);
end


FFT_Stack=fft(Image_Stack,[],1);

FFT_Stack(1:Low_Pass,:,:)=0;
FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
%FFT_Stack(High_Pass:end,:,:)=0;

%clear Image_Stack Image Image_Temp
%FFT_Stack(round(size(FFT_Stack,1)/2):end,:,:)=0;
Image_Stack_New=abs(ifft(FFT_Stack,[],1));
%plot(real(FFT_Stack(:,500,400)));
New_Max=max(Image_Stack_New(:));

Image_Stack_New_Normalized=Image_Stack_New./max(Image_Stack_New(:));


%for q=1:length(image_index)
%    Temp_Slides(:,:)=Image_Stack_New_Normalized(q,:,:);
%    imwrite(Temp_Slides,sprintf('Enface_Binned_by_%i_%i.png',Binning_Factor,image_index(q)),'Bitdepth',new_bit_depth);
%    disp(q);
%end


%for rr=1:size(Image_Stack_New_Normalized,3)
%    Temp_Bscan(:,:)=Image_Stack_New_Normalized(:,:,rr);
%    imwrite(Temp_Bscan,sprintf('Bscan_%i.png',rr),'Bitdepth',new_bit_depth);
%    disp(rr);
%end

%Temp_Sides_Read=imread(sprintf('%i_New.png',image_index(q)));








