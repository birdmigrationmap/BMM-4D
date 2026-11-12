function paths = nr_shared_path(paths)
%NR_SHARED_PATH Paths on the Vogelzug share, as this machine sees them (dev).
%   PATHS = NR_SHARED_PATH(PATHS) maps paths on the research-unit share
%   between the cluster (/mnt/ceph/DOM_Forschung/UNIT_Vogelzug) and the
%   laptop's mount (/Volumes/UNIT_Vogelzug), whichever exists here. The
%   shared night set's manifests store cluster paths. Other paths and ""
%   are returned unchanged.

cluster_root = "/mnt/ceph/DOM_Forschung/UNIT_Vogelzug";
laptop_root = "/Volumes/UNIT_Vogelzug";
paths = string(paths);
if isfolder(cluster_root)
    paths = replace(paths, laptop_root, cluster_root);
elseif isfolder(laptop_root)
    paths = replace(paths, cluster_root, laptop_root);
end
end
