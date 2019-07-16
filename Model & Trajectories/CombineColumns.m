%% Last Edit 26/04/2017

disp('Generating the linear combinations between the columns of the symbolic dymanic model ...');
Wc1 = Wi(:,IC1);
Wc1 = Wc1';

Wc2 = Wi(:,IC2);
Wc2 = Wc2';

Wc = Wc1 - M*Wc2;
Wc = Wc';
disp('Done.');
%% Write on a .m file the obtained simbolic matrix
disp('Write in Wc.m the result obtained in order to reduce computational cost for the trajectory optimization ...');
FP = fopen('Wc.m','w');
fprintf(FP,'Wci = [');
for i = 1:size(Wc,1)
    for j = 1:size(Wc,2)
        fprintf(FP,char(Wc(i,j)));
        if (j==size(Wc,2)) && (i ==size(Wc,1))
            break;
        end
        if j == size(Wc,2)
            fprintf(FP,';...\n');
        else
            fprintf(FP,',');
        end
    end
end    
fprintf(FP,'];');
fclose(FP);
disp('Done.');

clear Wc Wc1