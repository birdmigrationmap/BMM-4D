function scans_filepath = load_rd(scan_radar, scan_time, varargin)
%LOAD_RD Locate or download raw radar files for requested radar-time pairs.
%   SCANS_FILEPATH = LOAD_RD(SCAN_RADAR, SCAN_TIME, ...) returns a string
%   array with one filepath per requested radar-time combination.
%
%   Options:
%   - 'closest_accepted': maximum time offset in minutes.
%   - 'download': download missing files from NOAA when true.

p = inputParser;
addParameter(p,'closest_accepted',10); %
addParameter(p,'download',false); %
parse(p,varargin{:});

% List of unique date
[scan_day,~,scan_day_id] = unique(dateshift(scan_time,"start","day"));

% Initialize the vector
scans_filepath=strings(numel(scan_radar),numel(scan_time));

for i_s=1:numel(scan_radar)

    for i_d=1:numel(scan_day)

        % get folder location
        flr = "data/raw_scans/nexrad/"+scan_radar(i_s)+"_"+datestr(scan_day(i_d),'yyyymmdd');
        if ~isfolder(flr)
            mkdir(flr);
        end


        % List available file
        fl = dir(flr);
        fl = fl({fl.name}~="..");
        fl = fl({fl.name}~=".");
        fl = fl(~contains({fl.name},"_MDM"));
        fl = fl(~contains({fl.name},"DS_Store"));

        % Get date of the available time.
        fl_date = NaT(numel(fl),1);
        for i_f=1:numel(fl)
            tmp = aws_parse(fl(i_f).name);
            fl_date(i_f) = tmp.datetime;
        end


        id_t = find(scan_day_id==i_d);
        for i_t = 1:numel(id_t)
            % select the closest
            [tmp,id] = min(abs(fl_date-scan_time(id_t(i_t))));
            % check that it is close enough (30min)
            if numel(id)>0 && (minutes(tmp) < p.Results.closest_accepted)
                scans_filepath(i_s,id_t(i_t)) = fl(id).folder+"/"+fl(id).name;
            end
        end

        if p.Results.download && any(scans_filepath(i_s,id_t)=="")

            % Download list of available file (very slow)
            try
                files_day_info = aws_list_datetime(scan_radar(i_s), scan_day(i_d));
            catch err
                warning("load_rd:AwsListFailed", ...
                    "Could not list AWS files for %s on %s: %s", ...
                    scan_radar(i_s), datestr(scan_day(i_d), 'yyyy-mm-dd'), err.message)
                continue
            end
            valid_file = arrayfun(@(x) isscalar(x.datetime) && ~isnat(x.datetime), files_day_info);
            files_day_info = files_day_info(valid_file);

            if numel(files_day_info)==0
                warning("No data for "+scan_radar(i_s))
                continue
            end

            for i_t = 1:numel(id_t)
                % select only the closest
                [closest_dur,id]=min(abs([files_day_info.datetime]-scan_time(id_t(i_t))));

                % Check if the file is ok
                if (minutes(closest_dur) >= p.Results.closest_accepted)
                    warning("No file within %g minutes for %s at %s.", ...
                        p.Results.closest_accepted, scan_radar(i_s), string(scan_time(id_t(i_t))))
                    continue
                end

                tmp = strsplit(files_day_info(id).name,"/");
                if ~exist(flr+"/"+tmp{end},'file')
                    system( "/usr/local/bin/aws s3 cp s3://unidata-nexrad-level2/"+files_day_info(id).name+" "+flr+" --no-sign-request" );

                    if contains(tmp{end},".gz")
                        gunzip(flr+"/"+tmp{end}, flr)
                        delete(flr+"/"+tmp{end})
                    end
                end
            end

            warning("Missing files were downloaded for %s on %s. Re-run load_rd to refresh filepaths.", ...
                scan_radar(i_s), datestr(scan_day(i_d),'yyyy-mm-dd'))
        end
    end
end
