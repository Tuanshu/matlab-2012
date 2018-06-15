clear all

Lateral_res=30; %micron
Radius_of_Curvature=7.8;    %mm
Radius_of_Curvature_Reference=7;%mm
Lateral_Position=0:0.1:5;   %mm
Slope=-(Lateral_Position./(Radius_of_Curvature^2-Lateral_Position.^2).^0.5-Lateral_Position./(Radius_of_Curvature_Reference^2-Lateral_Position.^2).^0.5);
Corresponding_Axial_res=abs(Lateral_res*Slope);
plot(Lateral_Position,Corresponding_Axial_res);

