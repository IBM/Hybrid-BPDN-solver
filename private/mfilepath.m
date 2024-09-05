function p = mfilepath(filename)
% MFILEPATH   Gets the path of a full filename
%
%    P = MFILEPATH(FILENAME)  is intended for use in conjunction with
%    the MFILENAME('fullpath') built-in function and returns that path
%    part P.
%    This function is useful when reading or writing data or figure
%    files in the directory where the calling function is executed.
%    The path P ends with a FILESEP character.

%    Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%    This file is part of the code accompanying the paper: Joint-sparse
%    recovery from multiple measurements, submitted April 2009.

if exist('fileparts')
   [pathstr,name,ext] = fileparts(filename);
else
   idx = max([1,find(filename == filesep)]);
   pathstr = filename(1:idx-1);
end

p = [pathstr, filesep];
