%This code will be used to output newly processed proteiomics datasets! We
%will go from scratch for everything


addpath('/Users/pennerc/Documents/MATLAB/OlinkDataAnalysis')

%the first dataset was run on 3 different plates in two different batches,
%it seems like there was not a bridging sample between plates only between
%batches... I guess let's see? 

origDataRaw= readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/olink data 1/20210120_olinkdata.csv');
origData2use=origDataRaw(:,12:end);
origNames=origData2use.Properties.VariableNames;
origplate=origDataRaw.PlateID;
origBatch=origDataRaw.Batch;
origDaysPost=origDataRaw.DaysBetween;



intraSamplePlate1_1a=173; % index for RBM301-01 pooled plasma
intraSamplePlate1_1b=174; % index for RBM301-01 

intraSamplePlate1_2a = 87; %index for 143922
intraSamplePlate1_2b = 88; %index for 143923

inter_samplePlate1_Plate3_1a = 175; %index for RBM302-01 
inter_samplePlate1_Plate3_1b  = 176; %index for RBM302-02 


inter_samplePlate1_Plate3_2a = 77; %index for 142965
inter_samplePlate1_Plate3_2b = 78; %index for 142966

intraSamplePlate3a=182;
intraSamplePlate3b=183;

intraSamplePlate2_1a=179;
intraSamplePlate2_1b=180;

allSamp=origDataRaw.SampleID;


%% removing proteins with very high intrasample variation

%iterate through each protein 

intraSamplePlate1a=173; % index for RBM301-01 pooled plasma
intraSamplePlate1b=174; % index for RBM301-01 


[cvOlinkintraPlate1_1] = olinkCVGen(origData2use,intraSamplePlate1a, intraSamplePlate1b);
[cvOlinkintraPlate1_2] = olinkCVGen(origData2use,intraSamplePlate1_2a, intraSamplePlate1_2a);


[cvOlinkintraPlate2_1] = olinkCVGen(origData2use,intraSamplePlate2_1a, intraSamplePlate2_1b);

sum(cvOlinkintraPlate2_1> .2 & cvOlinkintraPlate1_2>.2)

[cvOlinkintraPlate3_1] = olinkCVGen(origData2use,intraSamplePlate2a, intraSamplePlate2b);





%% reorienting data correctly


%sigh so annoying they saved everything as NA instead of nan..........


%ok first orig dataset
allDatCount=nan(1,size(origData2use,2)); %so in the original dataset there are a ton of proteins with exactly 32-33 usable samples.... 
% not sure why, but I'm removing these from the analysis
for dd=1:size(origData2use, 2)
    count=0;
    for tt=1:size(origData2use,1)
        if iscell(origData2use{tt,dd})
            if strcmp(origData2use{tt,dd}, 'NA')
        count=count+1; %iterate through every cell in the table and mark how many nans there are
            end
        end
        allDatCount(dd)=count;
    end
end

origData2use(:,allDatCount>100)=[]; %for now I simply delete these proteins
origRemNames=origNames(allDatCount>150);
origNames(allDatCount>150)=[]; %this is our updated name list
origData2use=origData2use{:,:}; %now that it is a single data type we can easily create a table of values


[coeff,scoreTot,latent,tsquared,explained] = pca(origData2use(~isnan(origDaysPost),: )  );



figure
hold on
plot(  scoreTot( contains(origplate(~isnan(origDaysPost) ),'Plate1')  ,1)   , scoreTot( contains(origplate(~isnan(origDaysPost)),'Plate1') ,2)  ,'.' ,'color', 'r', 'MarkerSize',25)
plot(  scoreTot( contains(origplate(~isnan(origDaysPost)),'Plate2')  ,1)   , scoreTot( contains(origplate(~isnan(origDaysPost)),'Plate2') ,2)  ,'.' ,'color', 'g', 'MarkerSize',25)
 plot(  scoreTot( contains(origplate(~isnan(origDaysPost)),'Plate3')  ,1)   , scoreTot( contains(origplate(~isnan(origDaysPost)),'Plate3') ,2)  ,'.' ,'color', 'k', 'MarkerSize',25)

xlabel(['PCA1 ', num2str(explained(1)), ' variance explained' ] )
ylabel(['PCA2 ', num2str(explained(2)), ' variance explained' ] )

legend({'Plate1', 'Plate2', 'Plate3'})

figure
hold on
plot(  scoreTot( contains(origBatch,'Batch1')  ,1)   , scoreTot( contains(origBatch,'Batch1') ,2)  ,'.' ,'color', 'c', 'MarkerSize',25)
plot(  scoreTot( contains(origBatch,'Batch2')  ,1)   , scoreTot( contains(origBatch,'Batch2') ,2)  ,'.' ,'color', 'm', 'MarkerSize',25)

xlabel(['PCA1 ', num2str(explained(1)), ' variance explained' ] )
ylabel(['PCA2 ', num2str(explained(2)), ' variance explained' ] )



%% doing the same for replication dataset


repDataRaw= readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/olink data 2/20220726_olinkreplication.csv');
repData2use=repDataRaw(:,12:end);
repNames=repData2use.Properties.VariableNames;
repDaysPost=repDataRaw.DaysBetween;
repID=repDataRaw


boxKey= readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/KP_olinkreplicationFINAL Blinded Box Key 7.5.2022.xlsx');