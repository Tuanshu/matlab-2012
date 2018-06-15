clear all

%Too bad axial resolution
fitting=0;
new_bit_depth=8;
cd('D:\Users\TuanShu\121126_LCD2\');

Stage_speed=1;  %micron/sec
Sampling_rate=33;   %Hz

Frame_Axial_Spacing=Stage_speed/Sampling_rate*5/3;  %micron

Objective_Focal_Length=4/0.85;   %cm
Porjection_Lens_Focal_Length=50;    %cm
Pixel_Size=7.4;       %micron

%Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;
Lateral_Spacing=600/(538-330);

image_index=80:280;
ROI=[300 700;300 700];      %up, down, left, right
Image_Stack(length(image_index),1:(ROI(1,2)-ROI(1,1)+1),1:(ROI(2,2)-ROI(2,1)+1))=0;



upper_range_of_1st_frame_1st_line=500;

lower_range_of_1st_frame_1st_line=750;

select_height=600;

window_size=150; %2*window_size
window_size_1st_line=150; %2*window_size

upper_range=upper_range_of_1st_frame_1st_line;
lower_range=lower_range_of_1st_frame_1st_line;
upper_range_of_1st_line=upper_range_of_1st_frame_1st_line;
lower_range_of_1st_line=lower_range_of_1st_frame_1st_line;


for p=1:length(image_index)
    Image=imread(sprintf('Bscan_%i.png',image_index(p)));    
    
    for q=1:size(Image,2)
        if q==1
            [max_value max_index]=max(Image(upper_range_of_1st_line:lower_range_of_1st_line,q));
            real_max_index=(max_index+upper_range_of_1st_line-1);
        else
            [max_value max_index]=max(Image(upper_range:lower_range,q));
            real_max_index=(max_index+upper_range-1);
        end
        Image(:,q)=circshift(Image(:,q),select_height-real_max_index);

        if q==1
            upper_range_of_1st_line=real_max_index-window_size_1st_line;
            lower_range_of_1st_line=real_max_index+window_size_1st_line;
        end
        upper_range=max(real_max_index-window_size,1);
        lower_range=min(real_max_index+window_size,size(Image,1));
    end
    imwrite(Image,sprintf('Bscan_%i_Aligned.png',p),'Bitdepth',new_bit_depth);
    fprintf('%d\n',p);
end

