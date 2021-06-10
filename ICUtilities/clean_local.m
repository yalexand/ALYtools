function clean_local()
    myCluster = parcluster('local')
    delete(myCluster.Jobs);
end

