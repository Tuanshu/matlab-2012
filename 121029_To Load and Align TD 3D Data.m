clear all

cd('D:\Users\TuanShu\121025\3D\sample04_1');

file_name='sample04';
file_index=[1:200];


upper_range_of_1st_frame_1st_line=850;

lower_range_of_1st_frame_1st_line=900;

select_height=900;

window_size=20; %2*window_size
window_size_1st_line=20; %2*window_size

upper_range=upper_range_of_1st_frame_1st_line;
lower_range=lower_range_of_1st_frame_1st_line;
upper_range_of_1st_line=upper_range_of_1st_frame_1st_line;
lower_range_of_1st_line=lower_range_of_1st_frame_1st_line;

for p=1:length(file_index)
    Frame(:,:,p)=importdata(sprintf('%s_%i.txt',file_name,file_index(p)));
    for q=1:size(Frame,2)
        if q==1
            [max_value max_index]=max(Frame(upper_range_of_1st_line:lower_range_of_1st_line,q,p));
            real_max_index=(max_index+upper_range_of_1st_line-1);
        else
            [max_value max_index]=max(Frame(upper_range:lower_range,q,p));
            real_max_index=(max_index+upper_range-1);
        end
        Frame(:,q,p)=circshift(Frame(:,q,p),select_height-real_max_index);

        if q==1
            upper_range_of_1st_line=real_max_index-window_size_1st_line;
            lower_range_of_1st_line=real_max_index+window_size_1st_line;
        end
        upper_range=real_max_index-window_size;
        lower_range=real_max_index+window_size;
    end
    new_file_name=sprintf('%s_%i_Aligned.txt',file_name,file_index(p));
    dlmwrite(new_file_name,Frame(:,:,p),'delimiter','\t','newline','pc');
    fprintf('%d\n',p);
end
%%
 imagesc(Frame(:,:,50),'ydata',(1:size(Frame,1))*0.06,'xdata',(1:size(Frame,2))*5); %3V
%% 
 imagesc(Frame(:,:,200)); %3V
