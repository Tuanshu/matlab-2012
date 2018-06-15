clear all

%Too bad axial resolution
cd('D:\Users\TuanShu\121126_Eye\');
Start_index=200;
End_index=3500;

for p=1:fix((End_index-Start_index)/2)
    imwrite(fix(imread(sprintf('New_%i.png',Start_index+(p-1)*2))+imread(sprintf('New_%i.png',Start_index+(p-1)*2+1)))/2,sprintf('Reducedby2_%i.png',p),'Bitdepth',8);
    disp(p);
end

