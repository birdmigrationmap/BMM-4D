function vcp = read_vcp(raw_file)
%READ_VCP Scan strategy of a NEXRAD Archive II volume from its header.
%   VCP = READ_VCP(RAW_FILE) reads Message 5 (volume coverage pattern)
%   from the metadata record of RAW_FILE, the first bzip2 record, so only
%   that record is decompressed. VCP has:
%     number   - VCP number (e.g. 212, 215, 35, 31, 32);
%     cuts     - table, one row per elevation cut in scan order:
%                elev (deg), waveform (code: 1 contiguous surveillance,
%                2 contiguous Doppler with ambiguity resolution, 3 without,
%                4 batch, 5 staggered pulse pair), waveform_name, and
%                split (true for an elevation that also has a contiguous
%                surveillance cut: the cleaning reads that cut's reflectivity).
%   Returns number NaN and an empty table when Message 5 is not found.

decompressed = [tempname '.ar2'];
% First record only; older archives without bzip2 records: first 400 kB
python_code = "import bz2,struct,sys;f=open(sys.argv[1],'rb');o=open(sys.argv[2],'wb');f.read(24)" + newline + ...
    "if f.read(8)[4:7]!=b'BZh':" + newline + ...
    " f.seek(24);o.write(f.read(400000));sys.exit(0)" + newline + ...
    "f.seek(24);n=abs(struct.unpack('>i',f.read(4))[0]);o.write(bz2.decompress(f.read(n)))";
script_file = [tempname '.py'];
fid = fopen(script_file, "w"); fprintf(fid, "%s", python_code); fclose(fid);
[status, message] = system(sprintf('python3 "%s" "%s" "%s"', script_file, raw_file, decompressed));
delete(script_file)
if status ~= 0
    error("bmm4d:read_vcp:Decompress", "Decompression failed: %s", message)
end
fid = fopen(decompressed, "r", "ieee-be");
bytes = fread(fid, Inf, "uint8=>uint8");
fclose(fid);
delete(decompressed)

% Metadata messages fill 2432-byte slots: 12-byte CTM header, 16-byte
% message header (type in its 4th byte), then the message body
waveform_names = ["surveillance" "doppler_ar" "doppler" "batch" "staggered"];
vcp = struct("number", NaN, "cuts", table());
uint16_at = @(i) double(typecast(flip(bytes(i:i + 1)), "uint16"));
for position = 1:2432:numel(bytes) - 2432
    header = position + 12;
    if bytes(header + 3) ~= 5
        continue
    end
    body = header + 16;
    vcp.number = uint16_at(body + 4);
    n_cut = uint16_at(body + 6);
    % 22-byte pattern header, then 46 bytes per elevation cut
    first = body + 22 + 46 * (0:n_cut - 1)';
    elev = arrayfun(@(i) uint16_at(i), first) * 180 / 32768;
    waveform = double(bytes(first + 3));
    waveform_name = strings(n_cut, 1);
    known = waveform >= 1 & waveform <= 5;
    waveform_name(known) = waveform_names(waveform(known));
    split = arrayfun(@(e) any(waveform == 1 & abs(elev - e) < 0.05), elev);
    vcp.cuts = table(round(elev, 2), waveform, waveform_name, split, ...
        VariableNames=["elev" "waveform" "waveform_name" "split"]);
    return
end
end
