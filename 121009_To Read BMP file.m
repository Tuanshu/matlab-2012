clear all

cd('D:\Users\TuanShu\121008_Test5\');

Stage_speed=1;  %micron/sec
Sampling_rate=60;   %Hz

Frame_Axial_Spacing=Stage_speed/Sampling_rate;  %micron

Objective_Focal_Length=4;   %cm
Porjection_Lens_Focal_Length=50;    %cm
Pixel_Size=5;       %micron

Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;

image_index=1:700;
ROI=[1 480;1 640];      %up, down, left, right
Image_Stack(length(image_index),1:(ROI(1,2)-ROI(1,1)+1),1:(ROI(2,2)-ROI(2,1)+1))=0;
Low_Pass=60;   %pixel
High_Pass=120;
for p=1:length(image_index)
    Image=imread(sprintf('%i.png',image_index(p)));    
    %Image=imread(sprintf('1.jpg',image_index(p)));    
    Image_Stack(p,:,:)=Image(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),1);
    disp(p);
end

FFT_Stack=fft(Image_Stack,[],1);
FFT_Stack(1:Low_Pass,:,:)=0;
FFT_Stack(round(length(FFT_Stack)/2):end,:,:)=0;
Image_Stack_New=abs(ifft(FFT_Stack,[],1));
 plot(real(Image_Stack_New(:,300,150)));

[max_value max_index]=max(Image_Stack_New,[],1);
 
plot(real(FFT_Stack(:,300,150)));

Height(:,:)=Frame_Axial_Spacing*max_index(1,:,:);

imagesc(Height,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');


plot(0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),Height(100,:));
xlabel('(Micron)');
ylabel('(Micron)');