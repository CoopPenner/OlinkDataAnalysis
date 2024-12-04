function [cvOlink] = olinkCVGen(data2use,samp1, samp2)
%simple function to calculate coefficient of variation in service of qc




cvOlink=nan(1,size(data2use,2));
    for dd=1: size(data2use,2 )
        val1= data2use(samp1, dd);
        val2= data2use(samp2, dd);

    cvOlink(dd)= std([val1,val2])/  mean([val1,val2]) ;
    end
