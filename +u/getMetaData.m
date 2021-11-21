function SL = getMetaData(Monkey, overwrite)
% Imports metadata and saves a session list (SL)
% Monkey input is a cell array of monkeys you want to save data of

if(nargin<2)
    overwrite = 0;
end

%% Ubi
if(any(strcmp(Monkey,'Ubi')))
    % Read in data from Google Sheet
    vals = u.GetGoogleSpreadsheet('1JIfI_aEq1QYkWs_dEgET_NyEvcIn_7_TIcHD2w4PSug');
    if(~overwrite)
        load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat');
        oldSL = SL;
    end
    
    % Arrange into a struct. First row is the fields
    SL = struct([]);
    for i = 2:size(vals,1)
        if(~strcmp(vals(i,3),'Ubi') || strcmp(vals(i,3),'1') || (~overwrite && any(strcmp(extractfield(oldSL,'Date'),vals(i,1)))))
            continue;
        end
        if(strcmp(char(vals(i,1)),20170512) || strcmp(char(vals(i,1)),20170516)) % days with 2 files, need to fix first. 
            continue;
        end
        cur = length(SL)+1;
        for j = 1:size(vals,2)
            SL(cur).(char(vals(1,j))) = char(vals(i,j));
        end
    end
    
    % Rearrange the two fields to be a cell array and save other relevant data
    for i = 1:length(SL)
        
        SL(i).Session_Guger_Train = strsplit(SL(i).Session_Guger_Train,'_');
        SL(i).Stim_Trial_Start_Stop = strsplit(SL(i).Stim_Trial_Start_Stop,'_');
        
        D = SL(i).Date;
        S = SL(i).Session_Guger_Train;
        
        if(strcmp(char(S(2)),'NaN'))
            continue;
        end
        
        Session = [char(D),'_',char(S(2))];
        disp(['Session ',Session])
        % u.trainalign3(Session,2);    % Do this separately with a script, let it run overnight. replacing i16 files (train data)
        
        [accdata, trig1, trig2, lefttrials, righttrials, fs, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi([Session,'.i16']); % load data
        
        SL(i).StimHemi = 'L';
        SL(i).lefttrials = lefttrials;
        SL(i).righttrials = righttrials;
        SL(i).fs = fs;
        SL(i).accel_raw_r = double(accdata(:,2));
        SL(i).accel_raw_l = double(accdata(:,1));
        SL(i).lefttrialsuccess = lefttrialsuccess;
        SL(i).righttrialsuccess = righttrialsuccess;
        %         SL(i).cursorintarget = cursorintarget;
        
        SL(i).trig1 = trig1;
        SL(i).trig2 = trig2;
        
        % Has random triggers in the beginning that didn't actually cause
        % stimulation
        if(strcmp(SL(i).Date,'20170127'))
            SL(i).trig1(1:11) = [];
        end
    end
    
    SL = u.AppendThreshold(SL);
%     SL = a.RemoveTwitch(SL);
    SL = a.AppendReactionTimes(SL);
    SL = a.AppendNormalizedDelay(SL);
    
    if(~overwrite)
        SL = [oldSL,SL];
        clear oldSL;
    end
    
    SL = u.SortSL(SL,'Condition');
    
    % Save. v7.3 needed since SL is too large, but turns out v7.3 can't turn
    % off compression so saving / loading takes longer than usual.
    save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat','SL','-v7.3');
end

%% Igor
if(any(strcmp(Monkey,'Igor')))
    cd('F:\Dropbox\repos\abogaard\efetz\DualAccelRt\data')
    load_meta_info;
    
    deleteexp = [20111008;... % these experiments are weird for some reason or other
        20111019;...
        20111024;...
        201110242;...
        20111115;...
        20120214;...
        20120301;...
        20120419;...
        20120419;...
        20120625;...
        20120803;...
        20120807;...
        20120822;...
        20111122;... %gnd
        20111121;... %gnd
        20120112;... %gnd
        20120120;... %gnd, no twitch, bad
        20120828;... %weird accel, no twitch, funky data
        20120829;... %no twitch, weird short session, maybe wrong electrodes
        20120911;... %weird accel, weird data, no twitch
        20120906;... %noisy, bad performance
        20120816;... %bad performance, no twitch
        20120827;... %bad performance
        20120508;... %stim delivered all trials?
        20120512;... %weird triggering, no twitch
        20120926;... % tiny twitch
        ];
    
    skiprow = CheckForDataFiles(meta_info);
    
    cd('F:\Dropbox\repos\abogaard\efetz\DualAccelRt\code')
    load('SL.mat')
    SLIgor = SL;
    clear SL;
    
%     SLIgor(ismember(vertcat(SLIgor(:,1).expid), deleteexp),:) = []; % delete experiments we dont want        
        
    SL = struct([]);
    for i = 1:size(SLIgor,1)
        for j = 1:3
            if j == 1
                SL(i).Date = num2str(SLIgor(i,j).expid);
                ind = find([meta_info{:,2}]==SLIgor(i,j).expid);
                SL(i).Sessions = [char(meta_info{ind,3}),'_',char(meta_info{ind,4}),'_',char(meta_info{ind,5})];
                SL(i).Animal = 'Igor';
                SL(i).Stim_Amp = num2str(meta_info{ind,7});
                SL(i).Stim_Delay = num2str(meta_info{ind,end-1});
                SL(i).Stim_Loc = meta_info{ind,6};
             
                if(any(any(char(meta_info{ind,6}) == 'R')))
                    SL(i).StimHemi = 'R';
                elseif(any(any(char(meta_info{ind,6}) == 'L')))
                    SL(i).StimHemi = 'L';
                else
                    SL(i).StimHemi = 'NaN';
                end
            
                if(strcmp(SL(i).Date,'20120515'))
                    SL(i).StimHemi = 'R';
                elseif(strcmp(SL(i).Date,'20120603'))
                    SL(i).StimHemi = 'L';
                elseif(strcmp(SL(i).Date,'20120517'))
                    SL(i).StimHemi = 'L';
                elseif(strcmp(SL(i).Date,'20120519'))
                    SL(i).StimHemi = 'L';
                elseif(strcmp(SL(i).Date,'20120814'))
                    SL(i).StimHemi = 'L';
                elseif(strcmp(SL(i).Date,'20120830'))
                    SL(i).StimHemi = 'L';
                elseif(strcmp(SL(i).Date,'20120910'))
                    SL(i).StimHemi = 'L';
                elseif(strcmp(SL(i).Date,'20121002'))
                    SL(i).StimHemi = 'L';
                elseif(strcmp(SL(i).Date,'20120209'))
                    SL(i).StimHemi = 'L';
                end
                
                hand = 'X';
                if(any(char(meta_info{ind,9})=='R'))
                    hand = 'R';
                elseif(any(char(meta_info{ind,9})=='L'))
                    hand = 'L';
                end
                
                cond = 'NaN';
                if(~strcmp(SL(i).StimHemi,'NaN'))
                    if(strcmp(SL(i).StimHemi,hand))
                        cond = 'Ipsi';
                    else
                        cond = 'Contra';
                    end
                    if(strmatch('MOVE',char(meta_info{ind,9})))
                        cond = [cond,'_M'];
                    elseif(length(meta_info{ind,9}) > 1)
                        cond = 'NaN';
                    elseif(~any(char(meta_info{ind,9})=='R') && ~any(char(meta_info{ind,9})=='L'))
                        cond = SLIgor(i,j).cond;
                    else
                        cond = [cond,'_',num2str(meta_info{ind,end-1})];
                    end
                else
                    cond = SLIgor(i,j).cond;
                end
                
                SL(i).Condition = cond;
                
                SL(i).Notes = char(meta_info{ind,end});
                SL(i).Bad = [];
                SL(i).lefttrials = SLIgor(i,j).lefttrials;
                SL(i).righttrials = SLIgor(i,j).righttrials;
                SL(i).fs = SLIgor(i,j).fs;
                SL(i).accel_raw_r = double(SLIgor(i,j).accel_raw_r);
                SL(i).accel_raw_l = double(SLIgor(i,j).accel_raw_l);
                SL(i).lefttrialsuccess = SLIgor(i,j).lefttrialsuccess;
                SL(i).righttrialsuccess = SLIgor(i,j).righttrialsuccess;
                SL(i).trig1 = [];
                SL(i).trig2 = SLIgor(i,j).triggers;
            else
%                 if strcmp(SL(i).Date,'20111221')
%                     keyboard;
%                 end
                last = length(SL(i).accel_raw_r)*1000/SL(i).fs;
                leftbadind = SL(i).lefttrials(:,2)>last;
                SL(i).lefttrials(leftbadind,:) = [];
                SL(i).lefttrialsuccess(leftbadind) = [];
                rightbadind = SL(i).righttrials(:,2)>last;
                SL(i).righttrials(rightbadind,:) = [];
                SL(i).righttrialsuccess(rightbadind) = [];
                SL(i).trig2(SL(i).trig2>last) = [];
                SL(i).lefttrials = [SL(i).lefttrials;last+SLIgor(i,j).lefttrials];
                SL(i).righttrials = [SL(i).righttrials;last+SLIgor(i,j).righttrials];
                SL(i).accel_raw_r = [SL(i).accel_raw_r;double(SLIgor(i,j).accel_raw_r)];
                SL(i).accel_raw_l = [SL(i).accel_raw_l;double(SLIgor(i,j).accel_raw_l)];
                SL(i).lefttrialsuccess = [SL(i).lefttrialsuccess;SLIgor(i,j).lefttrialsuccess];
                SL(i).righttrialsuccess = [SL(i).righttrialsuccess;SLIgor(i,j).righttrialsuccess];
                if(j == 2)
                    SL(i).trig1 = last+SLIgor(i,j).triggers;
                end
                SL(i).trig2 = [SL(i).trig2;last+SLIgor(i,j).triggers];
            end
        
        end
    end
    
    clear SLIgor;
    
    SL = u.AppendThreshold(SL);
%     SL = a.RemoveTwitch(SL);
    SL = a.AppendReactionTimes(SL);
    SL = a.AppendNormalizedDelay(SL);
    
    SL = u.SortSL(SL,'Condition');
    
    save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat','SL','-v7.3');
    
%     for(i=1:size(meta_info,1))
%         if(~strcmp(meta_info{i,1},'Igor') || any(i==skiprow) || any(meta_info{i,2} == deleteexp))
%             continue;
%         end
%         
%         ind = length(SL)+1;
%         
%         for j = 1:3 % for pre, cond, post conditions
%             
%             fname = ['../data/' lower(meta_info{i,1}) '/' num2str(meta_info{i,2}) '_' meta_info{i,j+2}];
%             disp(['processing ' fname]);
%             
%             if ~isstr(meta_info{i,j+2}), disp(['skipping ' fname]), continue; end
%             [accdata, triggers, lefttrials, righttrials, fs, lefttrialsuccess, righttrialsuccess] = u.LoadTrainIgor(['../data/' lower(meta_info{i,1}) '/' num2str(meta_info{i,2}) '_' meta_info{i,j+2}]);
%             
%             if isempty(fs), continue; end
%             
%             if(j == 1)
%                 SL(ind).Date = num2str(meta_info{i,2});
%                 SL(ind).Animal = 'Igor';
%                 SL(ind).Stim_Amp = num2str(meta_info{i,7});
%                 SL(ind).Stim_Delay = num2str(meta_info{i,end-1});
%                 
%                 if(any(any(char(meta_info{i,6})=='R')))
%                     SL(i).StimHemi = 'R';
%                 elseif(any(any(char(meta_info{i,6})=='L')))
%                     SL(i).StimHemi = 'L';
%                 else
%                     SL(i).StimHemi = 'NaN';
%                 end
%                 
%                 hand = 'X';
%                 if(any(char(meta_info{i,9})=='R'))
%                     hand = 'R';
%                 elseif(any(char(meta_info{i,9})=='L'))
%                     hand = 'L';
%                 end
%                 
%                 cond = 'NaN';
%                 if(~strcmp(SL(i).StimHemi,'NaN'))
%                     if(strcmp(SL(i).StimHemi,hand))
%                         cond = 'Ipsi';
%                     else
%                         cond = 'Contra';
%                     end
%                     if(strmatch('MOVE',char(meta_info{i,9})))
%                         cond = [cond,'_M'];
%                     else
%                         cond = [cond,'_',num2str(meta_info{i,end-1})];
%                     end
%                 end
%                 
%                 %                 cond = 'NaN';
%                 %                 if(any(char(meta_info{i,9})=='R'))
%                 %                     SL(i).StimHemi = 'R';
%                 %                     if(any(char(meta_info{i,9})=='R'))
%                 %                         cond = 'Ipsi';
%                 %                     elseif(any(char(meta_info{i,9})=='L'))
%                 %                         cond = 'Contra';
%                 %                     end
%                 %                     if(~strcmp(cond,'NaN'))
%                 %                         if(strmatch('MOVE',char(meta_info{i,9})))
%                 %                             cond = [cond,'_M'];
%                 %                         else
%                 %                             cond = [cond,'_',num2str(meta_info{i,end-1})];
%                 %                         end
%                 %                     end
%                 %                 else(any(char(meta_info{i,9})=='L'))
%                 %                     SL(i).StimHemi = 'L';
%                 %
%                 %                     if(any(char(meta_info{i,9})=='L'))
%                 %                         cond = 'Ipsi';
%                 %                     elseif(any(char(meta_info{i,9})=='R'))
%                 %                         cond = 'Contra';
%                 %                     end
%                 %                     if(~strcmp(cond,'NaN'))
%                 %                         if(strmatch('MOVE',char(meta_info{i,9})))
%                 %                             cond = [cond,'_M'];
%                 %                         else
%                 %                             cond = [cond,'_',num2str(meta_info{i,end-1})];
%                 %                         end
%                 %                     end
%                 %                 end
%                 
%                 SL(ind).Condition = cond;
%                 SL(ind).Notes = char(meta_info{i,end});
%                 SL(ind).Bad = [];
%                 SL(ind).lefttrials = lefttrials;
%                 SL(ind).righttrials = righttrials;
%                 SL(ind).fs = fs;
%                 SL(ind).accel_raw_r = double(accdata(:,2));
%                 SL(ind).accel_raw_l = double(accdata(:,1));
%                 SL(ind).lefttrialsuccess = lefttrialsuccess;
%                 SL(ind).righttrialsuccess = righttrialsuccess;
%                 SL(ind).trig1 = [];
%                 SL(ind).trig2 = triggers;
%             else
%                 last = length(SL(ind).accel_raw_r);
%                 SL(ind).lefttrials = [SL(ind).lefttrials;last+lefttrials];
%                 SL(ind).righttrials = [SL(ind).righttrials;last+righttrials];
%                 SL(ind).accel_raw_r = [SL(ind).accel_raw_r;double(accdata(:,2))];
%                 SL(ind).accel_raw_l = [SL(ind).accel_raw_l;double(accdata(:,1))];
%                 SL(ind).lefttrialsuccess = [SL(ind).lefttrialsuccess;lefttrialsuccess];
%                 SL(ind).righttrialsuccess = [SL(ind).righttrialsuccess;righttrialsuccess];
%                 if(j == 2)
%                     SL(ind).trig1 = last+triggers;
%                 end
%                 SL(ind).trig2 = [SL(ind).trig2;last+triggers];
%             end
%             
%             
%         end
%     end

end

%% Kato
if(any(strcmp(Monkey,'Kato')))
    if(~exist('meta_info'))
        cd('F:\Dropbox\repos\abogaard\efetz\DualAccelRt\data')
        load_meta_info;
        skiprow = CheckForDataFiles(meta_info);
    end
    SL = struct([]);
    
    deleteexp = [20150108;%... % this one has no post conditioning
        20141219;... % weird pre average accel
        20141111;... % low trial number
        20141103;... % bad performance
        20150106;... % no stim artifact
        20150109;... % no stim artifact
        20150114;... % no stim artifact
        20141031;... % bad behavior
        20150116;... % no twitch
        ];
    
    for(i=1:size(meta_info,1))
        if(~strcmp(meta_info{i,1},'Kato') || any(i==skiprow))%  || any(meta_info{i,2} == deleteexp))
            continue;
        end
        
        ind = length(SL)+1;
        
        for j = 1:3 % for pre, cond, post conditions
            
            fname = ['../data/' lower(meta_info{i,1}) '/' num2str(meta_info{i,2}) '_' meta_info{i,j+2}];
            disp(['processing ' fname]);
            
            if ~isstr(meta_info{i,j+2}), disp(['skipping ' fname]), continue; end
            [accdata, triggers, lefttrials, righttrials, fs, lefttrialsuccess, righttrialsuccess] = u.LoadTrainKato(['../data/' lower(meta_info{i,1}) '/' num2str(meta_info{i,2}) '_' meta_info{i,j+2} '.i16']);
            
            if isempty(fs), continue; end
            
            if(j == 1)
                SL(ind).Date = num2str(meta_info{i,2});
                SL(ind).Sessions = [char(meta_info{i,j+2}),'_'];
                SL(ind).Animal = 'Kato';
                SL(ind).Stim_Amp = num2str(meta_info{i,7});
                SL(ind).Stim_Delay = num2str(meta_info{i,end-1});
                SL(ind).Stim_Loc = meta_info{i,6};
                SL(ind).StimHemi = 'R';
                
                cond = 'NaN';
                if(any(char(meta_info{i,9})=='R'))
                    cond = 'Ipsi';
                elseif(any(char(meta_info{i,9})=='L'))
                    cond = 'Contra';
                end
                if(~strcmp(cond,'NaN'))
                    if(strmatch('MOVE',char(meta_info{i,9})))
                        cond = [cond,'_M'];
                    else
                        cond = [cond,'_',num2str(meta_info{i,end-1})];
                    end
                end
                
                SL(ind).Condition = cond;
                SL(ind).Notes = char(meta_info{i,end});
                SL(ind).Bad = [];
                SL(ind).lefttrials = lefttrials;
                SL(ind).righttrials = righttrials;
                SL(ind).fs = fs;
                SL(ind).accel_raw_r = double(accdata(:,2));
                SL(ind).accel_raw_l = double(accdata(:,1));
                SL(ind).lefttrialsuccess = lefttrialsuccess;
                SL(ind).righttrialsuccess = righttrialsuccess;
                SL(ind).trig1 = [];
                SL(ind).trig2 = triggers;
            else
                last = length(SL(ind).accel_raw_r);
                SL(ind).Sessions = [SL(ind).Sessions,char(meta_info{i,j+2}),'_'];
                SL(ind).lefttrials = [SL(ind).lefttrials;last+lefttrials];
                SL(ind).righttrials = [SL(ind).righttrials;last+righttrials];
                SL(ind).accel_raw_r = [SL(ind).accel_raw_r;double(accdata(:,2))];
                SL(ind).accel_raw_l = [SL(ind).accel_raw_l;double(accdata(:,1))];
                SL(ind).lefttrialsuccess = [SL(ind).lefttrialsuccess;lefttrialsuccess];
                SL(ind).righttrialsuccess = [SL(ind).righttrialsuccess;righttrialsuccess];
                if(j == 2)
                    if(any(char(meta_info{i,9})=='R'))
                        SL(ind).trig1 = last+righttrials(:,1);
                    elseif(any(char(meta_info{i,9})=='L'))
                        SL(ind).trig1 = last+lefttrials(:,1);
                    end
                end
                SL(ind).trig2 = [SL(ind).trig2;last+triggers];
            end
            
            
        end
        
    end
    
    SL = u.AppendThreshold(SL);
    SL = a.AppendReactionTimes(SL);
    SL = a.AppendNormalizedDelay(SL);
    
    SL = u.SortSL(SL,'Condition');
    
    save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat','SL','-v7.3');
end


end

    function skiptheserows = CheckForDataFiles(meta_info)
        
        skiptheserows = [];
        
        for i = 1:size(meta_info,1)
            
            for j = 1:3
                
                err = 0;
                
                if ~ischar(meta_info{i,j+2}), continue; end
                
                fname = ['../data/' lower(meta_info{i,1}) '/' num2str(meta_info{i,2}) '_' meta_info{i,j+2}];
                
                if strcmp(meta_info{i,1}, 'Kato')
                    if ~exist([fname '.i16'], 'file'), disp(['missing ' fname '.i16']); err = 1; end
                    if ~exist([fname '.bin'], 'file'), disp(['missing ' fname '.bin']); err = 1; end
                    if ~exist([fname '.cfg'], 'file'), disp(['missing ' fname '.cfg']); err = 1; end
                end
                
                if strcmp(meta_info{i,1}, 'Igor')
                    if ~exist([fname '.bin'], 'file'), disp(['missing ' fname '.bin']); err = 1; end
                    if ~exist([fname '.cfg'], 'file'), disp(['missing ' fname '.log']); err = 1; end
                end
                
                if err==1, skiptheserows(end+1) = i; end
                
            end
        end
    end