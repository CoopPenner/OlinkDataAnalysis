%olink proteiomics
%% reading in data
dataRead=true;
%read in metaData
if dataRead

%reading in all of our data
DRSScoresOrig=readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/olink data 1/01_pd_drs.csv'); %DRS scores from original cohort (this includes PD AZ and HC)
DRSScoresRep=readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/olink data 2/olinkdata2_drs.csv'); %DRS scores from the rep, this should just be PD

origData=readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/olink data 1/20220517_olinkoriginalqualitycontrolvalues.csv');
repData=readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/olink data 2/20220729_olink2_QC.csv');
end

% I refer to orig and rep for this first part despite the fact that we are
% actually running them together as a function of improper matching... 

for condAtPlay=1:2 %sigh this is a bit sloppy I know

            if condAtPlay==1 %first orig and then rep
        DRS2use=DRSScoresOrig;
        data2use=origData;
            elseif condAtPlay==2
        DRS2use=DRSScoresRep;
        data2use=repData;
            end
    

  sampleDate=data2use.SampleDate; % the date a given sample was taken

    DRSTotAge=DRS2use.DRSTotalAge; %this will be our predictor variable (age adjusted total score)
    ID=DRS2use.INDDID; %INDDID from the DRS scores
    countID=data2use.INDDID; %INDDID from count data
    testDate=DRS2use.TestDate; %all of our test dates
    % Globaldx=DRSScores.GlobalDx;
    % Parkinsondx=DRSScores.PDCDx;

    %let's make sure our datapoint is represented in our count data
     
    
    noDate=isnat(testDate); ID(noDate)=[]; testDate(noDate)=[]; DRSTotAge(noDate)=[]; %when there is no date we remove that point
    
    [~,EachID]=findgroups(countID); %every participant in a given count dataset

[uniqueIDX] = unique(countID,'first');
idxToRemoveTot = find(not(ismember(1:numel(uniqueIDX),i))  ); %in some cases multiple vials were drawn on the same day and processed...

    
    slopeCollect=nan(1,length(countID)); %this will be the slope of a line fit to DRS score
    year2DecCollect=nan(1,length(countID)); %this will be the number of years before they initiate a two year decline in score ( what we term clinically meaningful)
    startScore=nan(1,length(countID)); % I'll probably exclude folks who had frank dementia at baseline
    sexCollect=nan(1,length(countID)); %collecting biological sex 
    ageCollect=nan(1,length(countID));
        
            for tt= 1:length(countID) %iterating through each patient that has a testDate



sampleDate2use= sampleDate(tt); %the date a sample was taken

if size(sampleDate2use,1)>1
    sampleDate2use=sampleDate2use(1,1); %sometimes this is duplicated, unsure why
end

             relDates=testDate(   ID==countID(tt)    ); %all the test dates
             relScores=DRSTotAge(ID==countID(tt)); %all the scores
             %sigh these are sometimes out of order  
            [testYear,testMonth]=ymd(relDates);
            [sampleYear,sampleMonth]=ymd(sampleDate2use);
if condAtPlay==2
    sampleYear=sampleYear+2000;
end

        
        pureTime=datenum(relDates); %convert it to pure time to sort
        pureTimeSample=datenum(sampleDate2use);
        [sorted,n]=sort(pureTime); %sorting on absolute value of time
        testYear=testYear(n); testMonth=testMonth(n); relScores=relScores(n); %then just outputting the reorder


        
        %now make every test date relative to the sample date

            testYear=testYear-sampleYear;
            testYear=testYear*12;
            testMonth=testMonth-sampleMonth;
            
            relativeTime=testYear+testMonth;

            preSamp=relativeTime<0;
            relativeTime(preSamp)=[]; relScores(preSamp)=[]; %remove any test dates that came before the sample was taken

        if ~isempty(relScores)

         [uniqueA, i,j] = unique(relativeTime,'first');
         idxToRemove = find(not(ismember(1:numel(relativeTime),i)))   ;
         relativeTime(idxToRemove)=[]; relScores(idxToRemove)=[]; %sometimes for some reason there are duplicate times
        relativeTime((isnan(relScores)))=[];
        relScores(isnan(relScores))=[]; 
        relativeTime=relativeTime/12; %everything contextualized in years, since that is the standard by which slopes are presented


        scoreVec=diff(relScores);  sigChange=find( (relScores-relScores(1) ) <=-2 ,1); %here we identify when scores sink two points below first measure.      
         P=polyfit(relativeTime,relScores,1);
        slopeCollect(tt)=P(1);
     
if isnan(P(1)) && condAtPlay==2
d=1
end

        if ~isempty(sigChange)      
    year2DecCollect(tt)= relativeTime(sigChange)  ;
        else
     year2DecCollect(tt)= nan ;
        end
    startScore(tt)=relScores(1); %collect starting point
    ageCollect(tt)= data2use.AGE__YRS(tt); %collect age at first check

    if strcmp(data2use.Sex(tt),'Male')
    sexCollect(tt)= 1; %coding male as 1
    elseif strcmp(data2use.Sex(tt),'Female')
    sexCollect(tt)= 2; %coding female as 2
    end

        end
            end

            if condAtPlay==1
        origSlope=slopeCollect;
        origID=countID;
        origDec=year2DecCollect;
        origDement=startScore;
        idxToRemoveOrig=idxToRemoveTot;
        origAge=ageCollect;
        origSex=sexCollect;
            elseif condAtPlay==2
        repSlope=slopeCollect;
        repID=countID;
        repDec=year2DecCollect;
        repDement=startScore;
        idxToRemoveRep=idxToRemoveTot;
        repAge=ageCollect;
        repSex=sexCollect;
            end
end



%% Plotting some basic clinical metrics

allDem=[origDement,repDement];
allAge=[origAge,repAge];
allSex=[origSex, repSex];


figure
allDec=[origDec, repDec];
histogram(allDec(allDem~=1),'BinWidth',.2)

figure
subplot(1,3,1)
allSlope=[origSlope, repSlope];
histogram(allSlope,'BinWidth',.2, 'FaceColor','b')
xlabel('DRS slope')


subplot(1,3,2)
histogram(allDem,'BinWidth',1, 'FaceColor','c')
vline(6,'r--')
vline(12,'m--')
xlabel('First DRS score after blood draw')


subplot(1,3,3)
histogram(allDec,'BinWidth',.3, 'FaceColor','m')
xlabel('Years from blood draw to DRS loss of 2')


figure

scatter(allDec,allSlope )


figure

scatter(allDec(allDem>6),allSlope(allDem>6))




%% reorienting data correctly

origData2use=origData(:,12:end);
origNames=origData2use.Properties.VariableNames;

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



% now rep dataset
%ok thank god this isn't some weird mixed table

repData2use=repData(:,11:end);
repNames=repData2use.Properties.VariableNames;
repData2use=repData2use{:,:};


%% unfortunately these datasets don't have the same proteins so we have to match them and only take proteins that are present in both

totalVarNames=[];
totalProtTable=[];
removeNames=[];

for dd=1:length(origNames)
    
    if(  sum((strcmp(repNames, origNames{dd}) ) )>0 )
    
        totalVarNames=[totalVarNames, {origNames{dd}}];
        totalProtTable=[totalProtTable, [origData2use(:,dd); repData2use(:, ( strcmp(repNames, origNames{dd})  ) )   ]  ]; %these are straightforward perfect matches

    elseif contains('_',origNames{dd}) % a lot of the ones that are missed have an underscore and a number or something
    origNames{dd}
    
        end
            elseif sum(contains(origRemNames, origNames{dd}))>0 % if this was a protein we know was explicitly removed in the previous step
                removeNames=[removeNames, {origNames{dd}}];
            
            end

end



repNames(contains(repNames, 'GZMA'))

totalSlope=[origSlope,repSlope];
totalDec=[origDec, repDec];
totalDement=[origDement,repDement];

%% normalizing and combining datasets 

[bridgeID,bridgePos]= intersect(repData.OlinkID, origData.OlinkID);



%checking to make sure that bridges are appropriate tho obvi too late to
%change anything lol





%first uploading the normalized run from Kristen


normRun=readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/2_SECOND TRY/220925_normalized_run1.csv');


%just checking how these samples went through


[origIDPres, origLoc]= intersect(normRun.SampleID, origData.OLINKID);








%% first let's actually experiment with random forest methods

%this code was generated with assistance from chatgpt

% Assume X is your proteomics data matrix and Y is a vector of labels (0 for controls, 1 for patients)
rfModel = TreeBagger(100, X, Y, 'OOBPrediction', 'on', 'Method', 'classification');

% Calculate feature importance
importanceScores = rfModel.OOBPermutedPredictorDeltaError;

% Sort and display the most important features
[sortedScores, featureIdx] = sort(importanceScores, 'descend');
topFeatures = featureIdx(1:10); % Top 10 most important proteins
disp('Top 10 important proteins:');
disp(topFeatures);







%% ok now to distribute our samples evenly across DRS slope, run set, age, sex, and potentially time to decrement. 
% Additionally I will remove patients who are demented at baseline at this
% stage.
 
totalAge=[origAge,repAge];
totalSex=[origSex, repSex];

totalrunNum=[ ones(1,length(origSlope)), 2* ones(1,length(repSlope))]; %we will balance using original run metric as well
totalrunNum(isnan(totalAge))=nan; %nanning out samples without a slope to ensure I don't get confused

pt2remove=  totalDement <=6 | (totalDec>=4 & totalSlope<-.6 );







%% initial PCA to evaluate data quality and matching... 


[coeff,scoreTot,latent,tsquared,explained] = pca(totalProtTable);


figure
plot(1:length(explained),explained)

ylabel('explained variance')
xlabel('Principal Component')


slope2use=totalSlope;
dement2use=totalDement<=6;

slopeCuf=-.5;



figure
hold on
plot(  scoreTot( totalSlope>slopeCuf  ,1)   , scoreTot( totalSlope>slopeCuf ,2)  ,'.' ,'color', 'r', 'MarkerSize',25)
plot(  scoreTot( totalSlope<=slopeCuf & totalDec<3 ,1)   , scoreTot( totalSlope<=slopeCuf & totalDec<3 ,2)  ,'.' ,'color', 'c', 'MarkerSize',25)
xlabel(['PCA1 ', num2str(explained(1)), ' variance explained' ] )
ylabel(['PCA2 ', num2str(explained(2)), ' variance explained' ] )

legend({'Slow progressors', 'Fast Progressors'})



figure
hold on
plot(  scoreTot( 1:length(origSlope) ,1)   , scoreTot( 1:length(origSlope) ,2)  ,'.' ,'color', 'r', 'MarkerSize',25)
plot(  scoreTot( length(origSlope)+1:end ,1)   , scoreTot( length(origSlope)+1:end ,2)  ,'.' ,'color', 'c', 'MarkerSize',25)
xlabel(['PCA1 ', num2str(explained(1)), ' variance explained' ] )
ylabel(['PCA2 ', num2str(explained(2)), ' variance explained' ] )

legend({'Original Data', 'Replication Data'})






%% looking at PCA1 loadings

pca1Loading=(coeff(:,1));

[~,posLoadIdx]=maxk(pca1Loading,20); %top 20 biggest contributors positive
[~,negLoadIdx]=mink(pca1Loading,20); %top 20 biggest contributors Neg 


figure
hold on
scatter(1:20,pca1Loading(posLoadIdx))
xVals=1.1:1:length(posLoadIdx)+1; loadings=pca1Loading(posLoadIdx); offset=.2;
text( xVals, loadings, (totalVarNames(posLoadIdx)  ))




pca2Loading=(coeff(:,2));

[~,posLoadIdx]=maxk(pca2Loading,20); %top 20 biggest contributors positive
[~,negLoadIdx]=mink(pca2Loading,20); %top 20 biggest contributors Neg 


figure
hold on
scatter(1:20,pca2Loading(posLoadIdx))
xVals=1.1:1:length(posLoadIdx)+1; loadings=pca2Loading(posLoadIdx); offset=.2;
text( xVals, loadings, (totalVarNames(posLoadIdx)  ))





hold on

scatter([40:-1:21],abs(pca1Loading(negLoadIdx)))
xVals=40.1:-1:21; loadings= abs(pca1Loading(negLoadIdx));
text( xVals, loadings, (totalVarNames(negLoadIdx)  ))

xlim([0,42])

legend({'20 Highest Ranking Positive Loading Coefficients'...
   '20 Highest Ranking Negative Loading Coefficients' })

ylabel('Absolute Value of Loading Coefficient')
a=gca; a.XTickLabel=[];


%% let's see if we can clean up this separation across PCA2 based on replicate number



[coeff,scoreOrig,latent,tsquared,explained] = pca(totalProtTable);


pca2Loading=(coeff(:,2));



[~,posLoadIdx]=maxk(pca2Loading,20); %top 20 biggest contributors positive
[~,negLoadIdx]=mink(pca2Loading,20); %top 20 biggest contributors Neg 

totProtAmmend= totalProtTable; 
totProtAmmend(:,[posLoadIdx(1:10)', negLoadIdx(1:10)' ] )=[];

[coeff,scoreOrig,latent,tsquared,explained] = pca(totProtAmmend);


figure
hold on
title('Proteiomics PCA after removing top 5 negative and positive loading variables')
plot(  scoreOrig( 1:length(origSlope) ,1)   , scoreOrig( 1:length(origSlope) ,2)  ,'.' ,'color', 'r', 'MarkerSize',25)
plot(  scoreOrig( length(origSlope)+1:end ,1)   , scoreOrig( length(origSlope)+1:end ,2)  ,'.' ,'color', 'c', 'MarkerSize',25)
xlabel(['PCA1 ', num2str(explained(1)), ' variance explained' ] )
ylabel(['PCA2 ', num2str(explained(2)), ' variance explained' ] )

legend({'Original Data', 'Replication Data'})




plot(scoreOrig(slope2use>slopeCuf & slope2use<0 &dement2use==0 & dec2use< 3  ,1),scoreOrig(slope2use>slopeCuf &slope2use<0 &dement2use==0 & dec2use< 3  ,2),'.' ,'color', 'c', 'MarkerSize',25)





dec2use=totalDec;
dement2use=totalDement;
slope2use=totalSlope;

slopeCuf=.1;


figure
hold on
plot(scoreOrig(dec2use<slopeCuf &dement2use==0  ,1),scoreOrig( dec2use<slopeCuf &dement2use==0,2)  ,'.' ,'color', 'r', 'MarkerSize',25)
plot(scoreOrig(dec2use>slopeCuf &slope2use<0.1 &dement2use==0 ,1),scoreOrig(dec2use>slopeCuf &slope2use<0.1 &dement2use==0 ,2),'.' ,'color', 'c', 'MarkerSize',25)



% figure
% hold on
% plot(scoreOrig(origDec==0&origSlope<0  ,1),scoreOrig( origDec==0&origSlope<0,2)  ,'.' ,'color', 'r', 'MarkerSize',25)
% plot(scoreOrig(origDec>=1 &origSlope<0 ,1),scoreOrig(origDec>=1 &origSlope<0 ,2),'.' ,'color', 'c', 'MarkerSize',25)
% 





%% lets look at our loadings
pca2Loading=(coeff(:,2));

[~,posLoadIdx]=maxk(pca2Loading,10); %top 20 biggest contributors positive
[~,negLoadIdx]=mink(pca2Loading,10); %top 20 biggest contributors Neg 









figure
hold on
scatter(1:20,pca1Loading(posLoadIdx))
xVals=1.1:1:length(posLoadIdx)+1; loadings=pca1Loading(posLoadIdx); offset=.2;
text( xVals, loadings, (origNames(posLoadIdx)  ))


hold on

scatter([40:-1:21],abs(pca1Loading(negLoadIdx)))
xVals=40.1:-1:21; loadings= abs(pca1Loading(negLoadIdx));
text( xVals, loadings, (origNames(negLoadIdx)  ))

xlim([0,42])

legend({'20 Highest Ranking Positive Loading Coefficients'...
   '20 Highest Ranking Negative Loading Coefficients' })

ylabel('Absolute Value of Loading Coefficient')
a=gca; a.XTickLabel=[];

%% simple initial eval



dec2use=totalDec;
dement2use=totalDement;
slope2use=totalSlope;
slopeCuf=-.75;
decCut=3;

% HitIdx=find(contains(VarNames, KristenHits));

HitIdx=I;



    for dd=1:length(HitIdx)
    varAtPlay=totalVarNames(HitIdx(dd));
    
    proteinCount=(totalProtTable(:,HitIdx(dd) ));
    





figure
hold on
     b=bar(1,   nanmean(     proteinCount(slope2use<=slopeCuf & ~isnan(slope2use) &dement2use==0 & dec2use<decCut)  )  )   ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel('Normalized Protein Level')
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(proteinCount(slope2use>slopeCuf  & ~isnan(slope2use) &dement2use==0 & dec2use<decCut )  ) ) ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


hold on
a=scatter(rand(1, sum(slope2use<=slopeCuf & ~isnan(slope2use) & dement2use==0& dec2use<decCut  )  )+.5, proteinCount(slope2use<=slopeCuf & ~isnan(slope2use)& dement2use==0& dec2use<decCut ), 'Marker', 'o' ) ;
a.CData=[.8 .2 .5]; 
a=scatter(rand(1, sum(slope2use>slopeCuf & ~isnan(slope2use)&dement2use==0& dec2use<decCut )  )+2.5, proteinCount(slope2use>slopeCuf & ~isnan(slope2use)&dement2use==0& dec2use<decCut ), 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];


    p=anovan( proteinCount(  ~isnan(slope2use))  , {slope2use( ~isnan(slope2use)  ) > slopeCuf & dec2use( ~isnan(totalSlope)  ) <decCut & dement2use(~isnan(slope2use))==0  }, 'display','off');


title([varAtPlay, ' p=', num2str(p)]) 

legend({'Fast Progressors (<-.75) ', 'Slow Progressors (>-.75)'})

    end









 %% removing patients listed in both datasets   
    % doubleIdx=false(1,length(origSlope));
    % 
    % 
    % 
    % 
    % 
    % 
    % for dd=1:length(origSlope) 
    %     if sum(slopeCollect==origSlope(dd))>=1
    %         doubleIdx(dd)=true;
    %     end
    % end

    
origID(origSlope==0)=[];
origSlope(origSlope==0)=[];



% 
% 
%     figure
%     hold on
% 
%     histogram(repSlope, 'BinWidth',.01)
%     histogram(origSlope, 'BinWidth',.01)
% 
% 
% figure
% hold on
%      b=bar(1,nanmean(origSlope))  ;
% b.FaceColor = 'flat';
% b.FaceAlpha=.3;
% b.BarWidth=1.5;
% b.CData(1,:) = [.8 .2 .5]; 
% ylabel('Age at Onset (years')
% a=gca; a.XTickLabel=[];
% hold on
% 
% b=bar(3,nanmean(repSlope))  ;
% b.FaceColor = 'flat';
% b.FaceAlpha=.3;
% b.BarWidth=1.5;
% b.CData(1,:) = [0 0.7 .25];
% 
% 
% 
% 
% hold on
% a=scatter(rand(1, length(origSlope)  )+.5, origSlope, 'Marker', 'o' );
% a.CData=[.8 .2 .5]; 
% b=scatter(rand(1, length(repSlope) )+2.5, repSlope, 'Marker', 'o' );
% b.CData(1,:) = [0 0.7 .25];
% 
% 
% slopeMark=zeros(1,length([repSlope,origSlope])); slopeMark(1:length(origSlope))=1;
% anovan([origSlope,repSlope], {slopeMark}  )
% 
% 
% 
% figure
% hold on
%      b=bar(1,sum(origSlope<-.1)/length(origSlope)  )  ;
%      b=bar(3,sum(repSlope<-.1)/ length(repSlope)  )  ;
% 
% 
% 
% [n,p]=prop_test( [sum(origSlope<-.1), sum(repSlope<-.1)], [length(origSlope), length(repSlope)], 'false')




%% rerunning 







VarNames=origData.Properties.VariableNames;
countIDAtPlay=origData.INDDID  ;
metaDataIDAtPlay=origID;
slopeAtPlay=origSlope;


slopeIdxOrig=nan(1,length(countIDAtPlay));

   

for tt=1:length(countIDAtPlay)
                if sum(metaDataIDAtPlay==countIDAtPlay(tt)  )   >0
                slopeIdxOrig(tt)= slopeAtPlay((metaDataIDAtPlay==countIDAtPlay(tt)   )  ) ;
                
                end             
end



VarNames=repData.Properties.VariableNames;
countIDAtPlay=repData.INDDID  ;
metaDataIDAtPlay=repID;
slopeAtPlay=repSlope;


slopeIdxRep=nan(1,height(repData));
   

for tt=1:height(countIDAtPlay)
                if sum(metaDataIDAtPlay==countIDAtPlay(tt))   >0
                slopeIdxRep(tt)= slopeAtPlay((metaDataIDAtPlay==countIDAtPlay(tt)));
                end             
end





    
    

repSlope2use=slopeIdxRep;
origSlope2use=slopeIdxOrig ;


 figure
hold on





     b=bar(1,nanmean(origSlope))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel('Total Age Adjusted DRS slope (points/year)')
a=gca; a.XTickLabel=[];
hold on




     b=bar(3,nanmean(origSlope2use))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.4 .1 .9]; 
ylabel('Total Age Adjusted DRS slope (points/year)')
a=gca; a.XTickLabel=[];
hold on


b=bar(5,nanmean(repSlope))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];




hold on
a=scatter(rand(1, length(origSlope)  )+.5, origSlope, 'Marker', 'o' );
a.CData=[.8 .2 .5]; 


b=scatter(rand(1, length(origSlope2use) )+2.5, origSlope2use, 'Marker', 'o' );
b.CData(1,:) = [.4 .1 .9];


b=scatter(rand(1, length(repSlope2use) )+4.5, repSlope2use, 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];


slopeMark=zeros(1,length([repSlope,origSlope2use])); slopeMark(1:length(origSlope2use))=1;
anovan([origSlope2use,repSlope], {slopeMark}  )


legend({'All INDD patients', 'Initial Cohort', 'Replication Cohort'})



title('Overview of DRS slopes between Olink Cohorts')

lowCut=-.75;

figure
hold on
     b=bar(1,(sum(origSlope<lowCut)/length(origSlope) *100 )  )  ;
     b=bar(3,(sum(origSlope2use<lowCut)/ sum(~isnan(origSlope2use)) *100 )  )   ;
     b=bar(5,(sum(repSlope2use<lowCut)/ sum(~isnan(repSlope2use) ) ) *100 )  ;


ylabel('Percent of Included Participants with a DRS slope <= -.1')
[n,p]=prop_test( [sum(origSlope2use<lowCut), sum(repSlope2use<lowCut)], [sum(~isnan(origSlope2use)), sum(~isnan(repSlope2use))  ], 'false')

a=gca; a.XTickLabel=[];

legend({'All INDD patients', 'Initial Cohort', 'Replication Cohort'})


%%

%% ok let's check out the actual proteiomics and the eleven selected hits specifically


origData=readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/olink data 1/20220517_olinkoriginalqualitycontrolvalues.csv');
repData=readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/olink data 2/20220729_olink2_QC.csv');

KristenHits={'GALNT10','PDGFC', 'SELL', 'PVALB', 'ICAM', 'KYNU' ,'SCF','FCRL5','PTN','VWC2','C2', 'GPNMB'}



VarNames=origData.Properties.VariableNames;
countIDAtPlay=origData.INDDID  ;
metaDataIDAtPlay=origID;
slopeAtPlay=origSlope;


%%
slopeCuf=-.75;

anovaGrab=[];
for dd=12:length(totalVarNames)
    varAtPlay=totalVarNames(dd);
    
    % proteinCount=table2array(origData(:,dd ));
    
    proteinCount=totalProtTable(:,dd);


    try
p=anovan( proteinCount(  ~isnan(totalSlope))  , {totalSlope( ~isnan(totalSlope)  ) > slopeCuf & totalDec( ~isnan(totalSlope)  ) <4 & totalDement(~isnan(totalSlope))==0  }, 'display','off');
anovaGrab(dd)=p;

    catch
        anovaGrab(dd)=nan;
    end

end

anovaGrab(1:11)=nan;



[M,I]=mink(anovaGrab,20)


totalVarNames(I)
anovaGrab(I)





