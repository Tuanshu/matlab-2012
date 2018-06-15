clear all;

Data3D(1:30,1:30,1:30)=0;   %axial, lateral, lateral

%1: Expansion
%0: Contraction

Array=[1 1 1 0 0 0 1 1 1 0 0 0];

Threshold=1;
%% Thresholding


Data3D = max(Data3D, Threshold);
Data3D (Data3D==Threshold) = 0;
Data3D (Data3D~=0) = 1;

%%

for j=1:length(array)
    
    DData3D_1(2:size(DData3D_1,1),:,:)=diff(Data3D,1,1);
    DData3D_2(:,2:size(DData3D_2,2),:)=diff(Data3D,1,2);
    DData3D_3(:,:,2:size(DData3D_3,3))=diff(Data3D,1,3);
    DData3D=DData3D_1+DData3D_2+DData3D_3;

    if Array == 1
        
        
    
    
    if Array == 1
        index_1=find(Data3D==1);
        [index_1_x index_1_y index_1_z]=ind2sub(index_1);
        for p=1:length(index_1)
            if Data3D()
        
while (p<size(Data3D,1) && p<size(Data3D,1) 