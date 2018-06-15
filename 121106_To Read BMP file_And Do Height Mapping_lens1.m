clear all

cd('D:\Users\TuanShu\121106_Lens1\');

Stage_speed=1;  %micron/sec
Sampling_rate=33;   %Hz

Frame_Axial_Spacing=Stage_speed/Sampling_rate*5/3;  %micron

Objective_Focal_Length=4/0.85;   %cm
Porjection_Lens_Focal_Length=50;    %cm
Pixel_Size=7.4;       %micron

%Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;
Lateral_Spacing=600/(538-330);

image_index=1:900;
ROI=[1 768;1 1024];      %up, down, left, right
Image_Stack(length(image_index),1:(ROI(1,2)-ROI(1,1)+1),1:(ROI(2,2)-ROI(2,1)+1))=0;
Low_Pass=20;   %pixel
High_Pass=400;
for p=1:length(image_index)
    Image=imread(sprintf('%i.png',image_index(p)));    
    %Image=imread(sprintf('1.jpg',image_index(p)));    
    Image_Stack(p,:,:)=Image(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),1);
    disp(p);
end

FFT_Stack=fft(Image_Stack,[],1);
FFT_Stack(1:Low_Pass,:,:)=0;
FFT_Stack(round(length(FFT_Stack)/2):end,:,:)=0;
Image_Stack_New=real(ifft(FFT_Stack,[],1));
%plot(real(FFT_Stack(:,300,400)));

[max_value max_index]=max(Image_Stack_New,[],1);

MaxValue(:,:)=max_value(1,:,:);

%plot(real(FFT_Stack(:,300,400)));

Height(:,:)=Frame_Axial_Spacing*max_index(1,:,:);
Height=max(max(Height))-Height;
[maxvalue1D maxindex1D]=max(Height(:));
Xgrid(1:size(Height,1),1:size(Height,2))=0;
Ygrid(1:size(Height,1),1:size(Height,2))=0;

for p=1:size(Height,2)
    Xgrid(:,p)=(p-1)*Lateral_Spacing;
end
for q=1:size(Height,1)
    Ygrid(q,:)=(q-1)*Lateral_Spacing;
end

Xmax=Xgrid(maxindex1D);
Ymax=Ygrid(maxindex1D);
FitSurface= fittype( @(c,r,a, b, x, y) c+(r^2-(x-a).^2-(y-b).^2).^0.5, 'independent', {'x', 'y'},'dependent', 'z');    %(x-A)^2+(y-B)^2+(z-C)^2=R^2
                                                                                                                        %z=C+(R^2-(x-A)^.2-(y-B)^.2).^0.5
                                                                                                                        
Starting_X=1650;                                                                                                                                                                                                                              
Starting_Y=1250;
R_Weight=500;
R_Weight_Inner=0;
Weight_Func=MaxValue/max(max(MaxValue));
Weight_Func(:,:)=1;
Weight_Func((Xgrid-Starting_X).^2+(Ygrid-Starting_Y).^2>R_Weight^2)=0;

Weight_Func((Xgrid-Starting_X).^2+(Ygrid-Starting_Y).^2<R_Weight_Inner^2)=0;

FitPara=fit([Xgrid(:),Ygrid(:)],Height(:),FitSurface,'Weight',Weight_Func(:),'StartPoint',[-15000,16000,Starting_X,Starting_Y]);

CheckCurve=FitPara.c+((FitPara.r)^2-(Xgrid-FitPara.a).^2-(Ygrid-FitPara.b).^2).^0.5;

imagesc(MaxValue/max(max(MaxValue)),'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
axis equal


imagesc(Weight_Func,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
axis equal

imagesc(CheckCurve,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
axis equal

plot((1:(size(Image_Stack_New,1)))*Frame_Axial_Spacing,real(Image_Stack_New(:,240,360)));
xlabel('Axial Position (Micron)');
ylabel('Amplitude (a.u.)');


imagesc(Height,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
axis equal

Height2=Height;
Height2((Xgrid-Starting_X).^2+(Ygrid-Starting_Y).^2>R_Weight^2)=0;
plot(0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),Height2(500,:),0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),CheckCurve(500,:));
xlabel('(Micron)');
ylabel('(Micron)');
legend('Measured height','Fitting curve')

Fitted_Curvature=FitPara.r %(micron)