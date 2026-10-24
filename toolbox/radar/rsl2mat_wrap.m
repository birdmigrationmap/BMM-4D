function radar = rsl2mat_wrap(filepath, params)
% Wrapper of rsl2mat adding a few improvements

if exist("rsl2mat", "file") ~= 3
    error("rsl2mat_wrap:MexNotFound", ...
        ["The compiled rsl2mat MEX file is missing. The current wsrlib checkout " ...
        "contains only rsl2mat.m, which is documentation and cannot be called. " ...
        "Rebuild wsrlib with wsrlib_make using a supported Intel MATLAB installation."])
end

% Get the radar name (callid from the filepath
[~,filename]=fileparts(filepath);
callid = char(filename);

% Ensure that filepath is a char
radar = rsl2mat(char(filepath), callid(1:4), params);

% Convert year, month, day etc.. to datetime and remove them
radar.datetime = datetime(radar.year,radar.month,radar.day,radar.hour,radar.minute,radar.second);
radar = rmfield(radar,{'year','month','day','hour','minute','second'});
end
