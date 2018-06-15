clear all

cd('D:\Users\TuanShu\121004\3D\sample04_5');

file_name='sample04';


upper_range_of_1st_frame_1st_line=880;

lower_range_of_1st_frame_1st_line=940;



    Frame=importdata(sprintf('%s_%i_Aligned.txt',file_name,1));
 imagesc(Frame);