b1 = 0.1434; b2 = 0.1391; f1 = 0.3302; f2 = 0.3576;
ml1 = 20; ml2 = 10; mm1 = 1; mm2 = 1;
I1 = 0.4667; I2 = 0.2333; Im2 = 21.18; Im1 = 12.10;
mr1 = 0.25*(ml1+1);mr2 = 0.25*ml2; 
gammaR_ref = [b1;...
                                   b2;...
                                   f1;...
                                   f2;...
                                  Im2;...
          ml1/2 + ml2/2 + mm2/2 + mr1;...
 I1 + Im1/100 - ml1/4 - ml2/4 - mm2/4;...
                          ml2/2 + mr2;...
                           I2 - ml2/4];
%{
[b1;...
                                    b2;...
                                    f1;...
                                    f2;...
                                   Im2;...
                  I2 + ml1 + mm2 + mr1;...
 I1/16 - I2/16 + Im1 - ml1/16 - mm2/16;...
                              ml2 - I2;...
                              I2 + mr2];
%}