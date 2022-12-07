function paths = findMatchingFiles(basepath)
% Syntax:       paths = findMatchingFiles(basepath);
% Description:  Return all files with the given base filename

% Parse base filepath
basepath = regexprep(basepath,'\','/');
[path name ext] = fileparts(basepath);
if ~isempty(path)
    path = [path '/'];
end

% Find matching filenames
list = dir([path name '*' ext]);
names = {list.name};
names = names(:);

%{
% Sort by appended number(s), if any
getNum = @(str) str2double(regexp(str,'\d+','match'));
nums = cell2mat(cellfun(getNum,names,'UniformOutput',false));
[~,idx] = sortrows(nums);
names = names(idx);
%}

% Return complete paths to matching files
fcn = @(name) [path name];
paths = cellfun(fcn,names,'UniformOutput',false);
