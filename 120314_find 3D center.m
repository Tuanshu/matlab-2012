clear all


cd('D:\');
Data=importdata('2D.txt');      % Data_2: the glass data 1
%Data(1:100,1:50)=0;
Distance_lateral_pixel=2;

%Data(45,20)=20;

%Data(1,1)=10;

%Data(99,40)=10;

Direction=1;    %Direction: 1 plus; 2 minus

R_start=100;
R_grid_rough=1;

x_start=40;
y_start=90;

x_delta=0.1;
y_delta=0.1;

X_grid_1D=(1:size(Data,2))*Distance_lateral_pixel;
Y_grid_1D=(1:size(Data,1))'*Distance_lateral_pixel;

Empty_grid_X(1:size(Data,1))=1;
Empty_grid_X=Empty_grid_X';
Empty_grid_Y(1:size(Data,2))=1;

X_grid=Empty_grid_X*X_grid_1D;

Y_grid=Y_grid_1D*Empty_grid_Y;

[value_1 index_1]=max(Data);
[value_2 index_2]=max(value_1);

max_index_x=index_2;
max_index_y=index_1(index_2);
max_cord_x=X_grid(max_index_y,max_index_x);
max_cord_y=Y_grid(max_index_y,max_index_x);
max_cord_z=Data(max_index_y,max_index_x);

loop=0;
R_current=R_start;
Stopping=0;
while Stopping < 3
    Merit=((((X_grid-max_cord_x).^2)+((Y_grid-max_cord_y).^2)+((Data-(max_cord_z-R_current)).^2)).^0.5-R_current).^2;
    Merit_sum=sum(sum(Merit));
    if loop == 0
        Merit_sum_previous=Merit_sum+1;
    end
    if (Direction == 1)
        if (Merit_sum<Merit_sum_previous)
            R_current=R_current+R_grid_rough;
            if Stopping > 0
                Stopping=Stopping-1;
            end
        else
            R_current=R_current-R_grid_rough;
            Direction=-1;
            Stopping=Stopping+2;
        end
    elseif (Direction == -1)
        if (Merit_sum<Merit_sum_previous)
            R_current=R_current-R_grid_rough;
            if Stopping > 0
                Stopping=Stopping-1;
            end
        else
            R_current=R_current+R_grid_rough;
            Direction=1;
            Stopping=Stopping+2;
        end
    end
    Merit_sum_previous=Merit_sum;
    loop=loop+1;
end

%%
x_guess=max_cord_x;
y_guess=max_cord_y;
z_guess=max_cord_z-R_current;
R_guess=R_current;