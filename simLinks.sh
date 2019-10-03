## Create simLinks for large folders

#mkdir -p /mnt/storage/cbarros/LandscapesInMotion/analyses

## move folders to storage
rsync -a /home/cbarros/GitHub/LandscapesInMotion/analyses/cache /mnt/storage/cbarros/LandscapesInMotion/analyses
rsync -a /home/cbarros/GitHub/LandscapesInMotion/analyses/DAfires_expAnalyses_files /mnt/storage/cbarros/LandscapesInMotion/analyses
rsync -a /home/cbarros/GitHub/LandscapesInMotion/analyses/DAfires_expAnalyses_cache /mnt/storage/cbarros/LandscapesInMotion/analyses
rsync -a /home/cbarros/GitHub/LandscapesInMotion/analyses/fireDataJoins /mnt/storage/cbarros/LandscapesInMotion/analyses
rsync -a /home/cbarros/GitHub/LandscapesInMotion/analyses/fireDataSummary /mnt/storage/cbarros/LandscapesInMotion/analyses
rsync -a /home/cbarros/GitHub/LandscapesInMotion/analyses/FireEvents /mnt/storage/cbarros/LandscapesInMotion/analyses
rsync -a /home/cbarros/GitHub/LandscapesInMotion/analyses/firesInFMAs /mnt/storage/cbarros/LandscapesInMotion/analyses

## after checking all is good remove the contents of the folders before creating the links

## create simLinks
ln -s /mnt/storage/cbarros/LandscapesInMotion/analyses/cache ~/GitHub/LandscapesInMotion/analyses/cache
ln -s /mnt/storage/cbarros/LandscapesInMotion/analyses/DAfires_expAnalyses_files /home/cbarros/GitHub/LandscapesInMotion/analyses/DAfires_expAnalyses_files
ln -s /mnt/storage/cbarros/LandscapesInMotion/analyses/DAfires_expAnalyses_cache /home/cbarros/GitHub/LandscapesInMotion/analyses/DAfires_expAnalyses_cache
ln -s /mnt/storage/cbarros/LandscapesInMotion/analyses/fireDataJoins /home/cbarros/GitHub/LandscapesInMotion/analyses/fireDataJoins
ln -s /mnt/storage/cbarros/LandscapesInMotion/analyses/fireDataSummary /home/cbarros/GitHub/LandscapesInMotion/analyses/fireDataSummary
ln -s /mnt/storage/cbarros/LandscapesInMotion/analyses/FireEvents /home/cbarros/GitHub/LandscapesInMotion/analyses/FireEvents
ln -s /mnt/storage/cbarros/LandscapesInMotion/analyses/firesInFMAs /home/cbarros/GitHub/LandscapesInMotion/analyses/firesInFMAs

