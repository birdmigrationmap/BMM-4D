function case_name = get_case_name(prefix, scans, i_case, i_sim)
%GET_CASE_NAME Build a deterministic file stem for one synthetic case.

case_stamp = string(datetime(scans.time(i_case), "Format", "yyyyMMdd_HHmmss"));
case_name = prefix + i_case + "_sim_" + i_sim + "_" + string(scans.station(i_case)) + "_" + case_stamp;
end
