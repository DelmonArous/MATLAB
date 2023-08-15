clear all;
close all;
fclose('all');
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
clc;

sourcepath = 'C:\Users\delmo\Dropbox\Jobb\Matlab\Olga ANOVA\Cytokines v2';
filelist = getAllFiles(sourcepath);

F = @(s)all(structfun(@(a)isscalar(a)&&isnan(a),s));

for i = 1:length(filelist)
    
    filelist{i}
    
    data = readXLSXdocument(filelist{i});
    [path, name, ext] = fileparts(filelist{i});
    
    y = [];
    X1 = {};
    X2 = [];
%     X3 = {};
    
    for j = 1:size(data,1)
        
        y    = [y; data{j,1}];
        X1{j} = data{j,2};
%         X2    = [X2; data{j,3}];
%         X3{j} = data{j,4};
        
    end
    
    X1 = X1';
    
    %% t-test
%     [h, p] = ttest(y, X1, 'Alpha', 0.05)
    
%     [p, tbl, stats] = anova1([y X1, X2], 'varnames', {'Dose','Time'})
%     [c, m, h, nms] = multcompare(stats)
    
    %% One-way ANOVA w/ Dose as covariate
    [p, tbl, stats] = anova1(y, X1, 'varnames', {'Dose'})
%     [c, m, h, nms] = multcompare(stats)

    %% One-way ANOVA w/ Time as covariate
%     [p, tbl, stats] = anova1(y, X2, 'varnames', {'Time'});

%% Two-way ANOVA
            
%   [p, tbl, stats] = anovan(y, {X1, X2, X3}, 'model', 'linear', ...
%         'varnames', {'Dose', 'Time', 'Gender'})
%     [p, tbl, stats] = anovan(y, {X1, X2, X3}, 'model', 'interaction', ...
%         'varnames', {'Dose', 'Time', 'Gender'});
%   [p, tbl, stats] = anovan(y, {X1, X2, X3}, 'model', 'full', ...
%         'varnames', {'Dose', 'Time', 'Gender'})

end
