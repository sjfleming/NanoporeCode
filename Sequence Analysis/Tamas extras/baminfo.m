function info = baminfo(filename, varargin)
%BAMINFO Information about BAM format file.
%
%   INFO = BAMINFO(FILENAME) returns a structure whose fields contain
%   information about a BAM file. FILENAME is a string containing a file
%   name, or a path and a file name, of a BAM file. INFO is a structure
%   with the following fields:
%
%
%       Filename    - name of the file
%       FilePath    - path to the file
%       FileModDate - modification date of the file
%       FileSize    - size of the file in bytes
%
%   Optional information contained in header section of BAM file.  Each of
%   these will be a structure with additional fields based on the
%   information present in the BAM file:
%
%       Header - File format version, sort order and group order
%
%       SequenceDictionary - Sequence Name, sequence length, genome
%       assembly identifier, MD5 checksum of sequence, URI of sequence,
%       species
%
%       ReadGroup - Read group identifier, sample, library, description,
%       platform unit, predicted median insert size, sequencing center,
%       date, platform
%
%       Program - Program name, version, command line
%
%       NOTE: Only the optional fields present in the BAM file will appear
%       in the output structure.
%
%   BAMINFO(...,'NUMOFREADS',T/F) scans the entire BAM file to determine
%   the number of alignment records stored in the file. The count is 
%   returned in the NumReads field. Default is false.
%
%   BAMINFO(...,'SCANDICTIONARY',T/F) scans the entire BAM file to
%   determine the available reference names and a tally for the alignment
%   records stored in the file. The information is returned in the
%   ScannedDictionary and ScannedDictionaryCount fields respectively.
%   Default is false. 
%
%   Example:
%
%       % Get BAM file information.
%       info = baminfo('ex1.bam');
%
%   See also BAMREAD.

%   Copyright 2009-2011 The MathWorks, Inc.

% Read input arguments:
[numReadsFlag, dictionaryFlag] = parse_inputs(varargin{:});

if ischar(filename)
    if ~((exist(filename,'file') || exist(fullfile(pwd,filename),'file')))
        error(message('bioinfo:baminfo:FileNotFound', filename));
    end
    % get the path
    filePath = fileparts(filename);
    % Get filename, filesize and date
    if isempty(filePath)
        filename = which(filename);
        filePath = fileparts(filename);
    end
    fInfo = dir(filename);
    info.Filename = fInfo.name;
    info.FilePath = filePath;
    info.FileSize = uint64(fInfo.bytes);
    info.FileModDate = fInfo.date;
    [headertext, references, lengths] = bioinfoprivate.bamaccessmex('getheader',[info.FilePath filesep info.Filename]);
    if ~isempty(headertext)
        headinfo = saminfo(headertext);
        fn = fields(headinfo);
        for i = 1:numel(fn)
            info.(fn{i}) = headinfo.(fn{i});
        end
    end
    NRefs = numel(references);
    if ~isfield(info,'SequenceDictionary') && NRefs
        % In case there was no header available, the SequenceDictionary is
        % taken from the "Reference" section of the BAM file
        info.SequenceDictionary(NRefs) = struct('SequenceName', [], 'SequenceLength', []);
        for ref_idx = 1:NRefs
            info.SequenceDictionary(ref_idx).SequenceName = references{ref_idx};
            info.SequenceDictionary(ref_idx).SequenceLength = int32(lengths(ref_idx));
        end
    end
else
    error(message('bioinfo:baminfo:InvalidInput'))
end

info.NumReads = uint64([]);
info.ScannedDictionary = cell(0,1);
info.ScannedDictionaryCount = uint64(zeros(0,1));

if numReadsFlag || dictionaryFlag
    [~,counts] = bioinfoprivate.bamaccessmex('baminfo',filename);
    if numReadsFlag
       info.NumReads = uint64(sum([0;counts]));
    end
    if dictionaryFlag && ~isempty(counts)
        if counts(1)==0
            info.ScannedDictionary = {info.SequenceDictionary(counts(2:end)>0).SequenceName}';
            info.ScannedDictionaryCount = uint64(counts(2:end));
        else
            info.ScannedDictionary = [{info.SequenceDictionary(counts(2:end)>0).SequenceName} {'Unmapped'}]';
            info.ScannedDictionaryCount = uint64(counts([2:end 1]));
        end
    end
end

%--------------------------------------------------------------------------
function [numReadsFlag, dictionaryFlag] = parse_inputs(varargin)
% Parse input PV pairs.

% defaults
numReadsFlag = 0;
dictionaryFlag = 0;

if rem(nargin,2) ~= 0
    error(message('bioinfo:baminfo:IncorrectNumberOfArguments', mfilename));
end

okargs = {'numofreads','scandictionary'};

for j=1:2:nargin-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    if k == 1
        numReadsFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    elseif k==2
        dictionaryFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end
