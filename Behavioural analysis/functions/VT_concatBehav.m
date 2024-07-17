function [dat] = VT_concatBehav(sub,sess, exp, dat)

thisPath = ['P' sub '\S' num2str(sess) '\Behaviour\'];
fullPath = ([exp.dataPath, thisPath]);

files = dir(fullPath);
files = files(contains({files.name}, '.mat'));
files = {files.name};
thisBehav = [];
for f = 1:length(files)
    thisFile = files{f};
    cd(fullPath)
    if exist(thisFile, 'file')

        load([fullPath thisFile]);
        thisBlock = str2num(thisFile(end-4)); %Get block code 
        blockInfo = [repmat([str2num(sub), sess, thisBlock],length(Behav),1), [1:length(Behav)]'];
        thisBehav = [blockInfo, Behav];
      
        dat.allBehav = [dat.allBehav; thisBehav];   
        dat.sumBehav = [dat.sumBehav; str2num(sub), sess, thisBlock, mean(thisBehav(:,7))*100];
    end
end

if ~isempty(thisBehav)
disp(['P' sub ' S' num2str(sess), ', Accuracy = ' num2str(nanmean(thisBehav(:,7))*100) '%']);
dat.allBehav = dat.allBehav; 
dat.sumBehav = dat.sumBehav;
end

end
