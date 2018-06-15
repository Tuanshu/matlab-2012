clear all

cd('D:\Users\TuanShu\121004\3D\sample04_5');
%axial pixel resolution:0.06 micron
%lateral pixel resolution: 10 micron

file_name='sample04';
file_index=1:100;


upper_range_of_1st_frame_1st_line=880;

lower_range_of_1st_frame_1st_line=940;

select_height=900;

window_size=20; %2*window_size
window_size_1st_line=20; %2*window_size

upper_range=upper_range_of_1st_frame_1st_line;
lower_range=lower_range_of_1st_frame_1st_line;
upper_range_of_1st_line=upper_range_of_1st_frame_1st_line;
lower_range_of_1st_line=lower_range_of_1st_frame_1st_line;
Height=800;
Enface(1:400,1:100)=0;
for p=1:length(file_index)
    Frame=importdata(sprintf('%s_%i_Aligned.txt',file_name,file_index(p)));
    Enface(:,p)=Frame(800,:);
end
 imagesc(Enface);