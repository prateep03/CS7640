function paramValue = readParam(options, paramName, paramValue, mandatory)


if nargin < 4,
    mandatory = 0;
end

if isfield(options, paramName)
    paramValue = eval(['options.' paramName ';']);
elseif mandatory
    error(['Param : options.',paramName,' needed but not provided.']);
end