clear all

cd('D:\Users\TuanShu\121004\3D\sample04_5');
%axial pixel resolution:0.06 micron
%lateral pixel resolution: 10 micron

file_name='sample04';
file_index=80;


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
    Frame=importdata(sprintf('%s_%i_Aligned.txt',file_name,file_index));

 imagesc(Frame(end:-1:1,:),'ydata',(1:size(Frame,1))*0.06-40,'xdata',(1:size(Frame,2))*10-500);
 xlabel('Lateral Position (micron)');
 ylabel('Optical Path Difference (micron)');
%axis equal