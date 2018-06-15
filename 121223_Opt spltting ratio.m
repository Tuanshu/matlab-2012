clear all
E_ASE=1;
E_eye=1;
E_ASE_tune=E_eye:0.001:10;
r=0:0.001:1;
S_inter=(r.^2).*((1-r).^2)*E_ASE^2;
r_max=E_eye/E_ASE;
Array(1:length(r))=0;
Array(find(r>r_max,1,'first'))=max(S_inter);
plot(r,S_inter);
xlabel('reflectivity');
ylabel('S_i_n_t_e_r');

S_inter_toASE=((E_eye./E_ASE_tune).^2).*((1-(E_eye./E_ASE_tune)).^2).*E_ASE_tune.^2;

plot(E_ASE_tune,S_inter_toASE);

xlabel('E_A_S_E (unit: E_e_y_e)');
ylabel('S_i_n_t_e_r');

[value index]=max(S_inter_toASE);
E_ASE_tune(index)