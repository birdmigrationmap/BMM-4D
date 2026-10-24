function sweeps = read_level2_moments(raw_file, moment_names, max_elev)
%READ_LEVEL2_MOMENTS Read chosen moments from a NEXRAD Archive II file.
%   SWEEPS = READ_LEVEL2_MOMENTS(RAW_FILE, MOMENT_NAMES, MAX_ELEV) reads
%   the Message 31 radials of RAW_FILE (e.g. MOMENT_NAMES = ["REF" "CFP"])
%   for elevation cuts below MAX_ELEV degrees. rsl2mat does not expose
%   every moment (for example CFP, the clutter filter power removed), so
%   this reads them directly. The bzip2 records are decompressed by a
%   short Python (standard library) call; the messages are parsed here.
%   SWEEPS(k) has elev (mean of the radials), elevation_number, az [1 n_az],
%   and for each moment a [n_gate n_az] matrix (empty if the cut lacks it) of decoded values with its
%   range (m) [n_gate 1] in <name>_range. Codes: decoded values are
%   (raw - offset) / scale; raw 0 (below threshold) and 1 (range folded)
%   become NaN, and the raw codes are kept in <name>_code (uint8/uint16).
%   Several cuts can share an elevation (split cuts); they stay separate.

decompressed = [tempname '.ar2'];
python_code = "import bz2,struct,sys;f=open(sys.argv[1],'rb');o=open(sys.argv[2],'wb');f.read(24)" + newline + ...
    "while True:" + newline + ...
    " h=f.read(4)" + newline + ...
    " if len(h)<4: break" + newline + ...
    " n=abs(struct.unpack('>i',h)[0]);o.write(bz2.decompress(f.read(n)))";
script_file = [tempname '.py'];
fid = fopen(script_file, "w"); fprintf(fid, "%s", python_code); fclose(fid);
[status, message] = system(sprintf('python3 "%s" "%s" "%s"', script_file, raw_file, decompressed));
delete(script_file)
if status ~= 0
    error("read_level2_moments: decompression failed: %s", message)
end
fid = fopen(decompressed, "r", "ieee-be");
bytes = fread(fid, Inf, "uint8=>uint8");
fclose(fid);
delete(decompressed)

% Walk the messages: 12-byte CTM header, then a 16-byte message header.
% Message 31 has its own length; other messages fill 2432-byte slots.
radials = struct("elevation_number", {}, "elev", {}, "az", {}, "moments", {});
position = 1;
n_bytes = numel(bytes);
while position + 28 <= n_bytes
    header = position + 12;
    message_halfwords = double(typecast(flip(bytes(header:header + 1)), "uint16"));
    message_type = bytes(header + 3);
    if message_type == 31
        body = header + 16;
        azimuth = double(typecast(flip(bytes(body + 12:body + 15)), "single"));
        elevation_number = double(bytes(body + 22));
        elevation = double(typecast(flip(bytes(body + 24:body + 27)), "single"));
        if elevation < max_elev
            block_count = double(typecast(flip(bytes(body + 30:body + 31)), "uint16"));
            moments = struct();
            for i_block = 1:block_count
                offset = double(typecast(flip(bytes(body + 32 + 4 * (i_block - 1):body + 35 + 4 * (i_block - 1))), "uint32"));
                block = body + offset;
                if bytes(block) ~= uint8('D')
                    continue
                end
                name = string(char(bytes(block + 1:block + 3)'));
                name = strtrim(name);
                if ~ismember(name, moment_names)
                    continue
                end
                n_gate = double(typecast(flip(bytes(block + 8:block + 9)), "uint16"));
                first_gate = double(typecast(flip(bytes(block + 10:block + 11)), "uint16"));
                gate_spacing = double(typecast(flip(bytes(block + 12:block + 13)), "uint16"));
                word_size = double(bytes(block + 19));
                scale = double(typecast(flip(bytes(block + 20:block + 23)), "single"));
                value_offset = double(typecast(flip(bytes(block + 24:block + 27)), "single"));
                data_start = block + 28;
                if word_size == 8
                    code = bytes(data_start:data_start + n_gate - 1);
                else
                    code = typecast(flip(bytes(data_start:data_start + 2 * n_gate - 1)), "uint16");
                    code = flip(code);
                end
                moments.(name) = struct("code", code(:), "scale", scale, "offset", value_offset, ...
                    "range", first_gate + gate_spacing * (0:n_gate - 1)');
            end
            radials(end + 1) = struct("elevation_number", elevation_number, "elev", elevation, ... %#ok<AGROW>
                "az", azimuth, "moments", moments);
        end
        position = header + 2 * message_halfwords;
    else
        position = position + 2432;
    end
end

% Group the radials by elevation cut and sort by azimuth
cut_numbers = unique([radials.elevation_number]);
sweeps = cell(1, numel(cut_numbers));
for i_cut = 1:numel(cut_numbers)
    cut = radials([radials.elevation_number] == cut_numbers(i_cut));
    [az, order] = sort([cut.az]);
    cut = cut(order);
    sweep = struct("elevation_number", cut_numbers(i_cut), "elev", mean([cut.elev]), "az", az);
    for name = moment_names
        % Every sweep gets every moment field (empty if absent, e.g. the
        % Doppler cut of a split cut has no CFP)
        sweep.(name) = [];
        sweep.(name + "_code") = [];
        sweep.(name + "_range") = [];
        has_moment = arrayfun(@(r) isfield(r.moments, name), cut);
        if ~any(has_moment)
            continue
        end
        first = cut(find(has_moment, 1)).moments.(name);
        n_gate = numel(first.range);
        code = zeros(n_gate, numel(cut), "like", first.code);
        for i_radial = find(has_moment)
            radial_code = cut(i_radial).moments.(name).code;
            code(1:min(n_gate, numel(radial_code)), i_radial) = radial_code(1:min(n_gate, numel(radial_code)));
        end
        value = (double(code) - first.offset) / first.scale;
        value(code <= 1) = nan;
        sweep.(name) = value;
        sweep.(name + "_code") = code;
        sweep.(name + "_range") = first.range;
    end
    sweeps{i_cut} = sweep;
end
sweeps = [sweeps{:}];
end
