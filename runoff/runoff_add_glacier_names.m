clear all; close all; clc;

i = 0;
%i = i + 1; inputMatFiles{i} = 'MAR3.9-MIROC5-rcp26-JJA_mean-tidewaterbasins_rignotid.mat';
i = i + 1; inputMatFiles{i} = 'MAR3.9-MIROC5-rcp85-JJA_mean-tidewaterbasins_rignotid.mat';
%i = i + 1; inputMatFiles{i} = 'MAR3.9-NorESM1-rcp85-JJA_mean-tidewaterbasins_rignotid.mat';
%i = i + 1; inputMatFiles{i} = 'MAR3.9-CSIRO-Mk3.6-rcp85-JJA_mean-tidewaterbasins_rignotid.mat';
%i = i + 1; inputMatFiles{i} = 'MAR3.9-HadGEM2-ES-rcp85-JJA_mean-tidewaterbasins_rignotid.mat';
%i = i + 1; inputMatFiles{i} = 'MAR3.9-IPSL-CM5-MR-rcp85-JJA_mean-tidewaterbasins_rignotid.mat';

for iFile = 1:numel(inputMatFiles)
   inputMatFile = inputMatFiles{iFile};
   runoffIn = load(inputMatFile);
   
   glacierNameCsvFile = 'RignotMouginot2012_TableS1_GlacierList.csv';
   fid = fopen(glacierNameCsvFile,'rt');
   fgets(fid);
   glacierNames = textscan(fid,'%d%q%f%f','delimiter',',');
   fclose(fid);
   nnames = length(glacierNames{1});
   
   % Find all basin numbers from the runoffIn struct
   fn = fieldnames(runoffIn.runoff);
   basinNums = [];
   for i = 1:length(fn)
      basinName = fn{i};
      if strfind(basinName, 'basin')
         basinNums(end+1) = str2num(strrep(basinName,'basin',''));
      end
   end
   nbasins = length(basinNums);
   
   % Initialize output struct array
   runoff(nbasins).rignotGlacierID    = [];
   runoff(nbasins).rignotGlacierName  = '';
   runoff(nbasins).glacierX  = [];
   runoff(nbasins).glacierY  = [];
   runoff(nbasins).time = [];
   runoff(nbasins).timeUnits = '';
   runoff(nbasins).runoff = [];
   runoff(nbasins).runoffUnits = '';
   
   % Populate output struct array
   for ibasin = 1:nbasins
      runoff(ibasin).rignotGlacierID = basinNums(ibasin);
      for iname = 1:nnames
         if runoff(ibasin).rignotGlacierID == glacierNames{1}(iname)
            runoff(ibasin).rignotGlacierName = glacierNames{2}{iname};
            [x,y] = polarstereo_fwd(glacierNames{3}(iname), glacierNames{4}(iname), [], [], 70, -45);
            runoff(ibasin).glacierX  = x;
            runoff(ibasin).glacierY  = y;
            break
         end
      end
      runoff(ibasin).time = runoffIn.time.time;
      runoff(ibasin).timeUnits = runoffIn.time.units;
      runoff(ibasin).runoff = eval( sprintf('runoffIn.runoff.basin%d', basinNums(ibasin)) );
      runoff(ibasin).runoffUnits = runoffIn.runoff.units;
   end
   
   % Sort the struct array by rignotGlacierID field
   runoff_sorted(nbasins).rignotGlacierID    = [];
   runoff_sorted(nbasins).rignotGlacierName  = '';
   runoff_sorted(nbasins).glacierX  = [];
   runoff_sorted(nbasins).glacierY  = [];
   runoff_sorted(nbasins).time = [];
   runoff_sorted(nbasins).timeUnits = '';
   runoff_sorted(nbasins).runoff = [];
   runoff_sorted(nbasins).runoffUnits = '';
   
   rignotGlacierIDs = [runoff(:).rignotGlacierID];
   [rignotGlacierIDs_sorted, rignotGlacierIDs_sorted_idx] = sort(rignotGlacierIDs);
   
   i = 1;
   for ibasin = rignotGlacierIDs_sorted_idx
      runoff_sorted(i).rignotGlacierID   = runoff(ibasin).rignotGlacierID;
      runoff_sorted(i).rignotGlacierName = runoff(ibasin).rignotGlacierName;
      runoff_sorted(i).glacierX          = runoff(ibasin).glacierX;
      runoff_sorted(i).glacierY          = runoff(ibasin).glacierY;
      runoff_sorted(i).time              = 1950:2100; %runoff(ibasin).time;
      runoff_sorted(i).timeUnits         = 'year'; %runoff(ibasin).timeUnits;
      runoff_sorted(i).runoff            = runoff(ibasin).runoff;
      runoff_sorted(i).runoffUnits       = runoff(ibasin).runoffUnits;
      i = i + 1;
   end
   
   runoff = runoff_sorted;
   
   save(strrep(inputMatFile, '.mat', '_withGlacierIDs.mat'), 'runoff');

end

