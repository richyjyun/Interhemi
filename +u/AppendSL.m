clear; close all;

load('E:\U\Code\+u\MetaData.mat');  % Once done, instead of loading just have it as a function input

for i = 1:length(SL)
   SL(i).monkey = 'Ubi';
   D = SL(i).Date;
   S = SL(i).Session_Guger_Train;
   
   if(strcmp(char(S(2)),'NaN'))
       continue;
   end
   
   Session = [char(D),'_',char(S(2))];
   disp(['Session ',Session])
   % u.trainalign3(Session,2);    % Do this separately with a script, let it run overnight. replacing i16 files (train data)
   
   [accdata, trig1, trig2, lefttrials, righttrials, fs] = u.LoadTrain([Session,'.i16']); % load data
   
   SL(i).lefttrials = lefttrials;
   SL(i).righttrials = righttrials;
   SL(i).fs = fs;
   SL(i).accel_raw_r = accdata(:,2);
   SL(i).accel_raw_l = accdata(:,1);
   
   SL(i).trig1 = trig1;
   SL(i).trig2 = trig2;
   
end

save('E:\U\Code\+u\MetaData.mat','SL');