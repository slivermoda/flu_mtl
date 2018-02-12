function Display_RES(RES, Method)
fprintf( strcat('RMSE_', Method, ' = %f\n'), RES);
% fprintf( strcat('Acc_', Method, ' = %f +- %f\n'), mean(RES.Acc), std(RES.Acc) );
% fprintf( strcat('Sen_', Method, ' = %f +- %f\n'), mean(RES.Sen), std(RES.Sen) );
% fprintf( strcat('Spe_', Method, ' = %f +- %f\n'), mean(RES.Spe), std(RES.Spe) );


