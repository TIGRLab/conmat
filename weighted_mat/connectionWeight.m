% w = 2 / (n_roi_1 + n_roi_2) * n_fibers

% this function outputs the weighted connectivity matrices in a csv for length_normalized and volume_normalized connection density
% RUN voxelCount.m FIRST!

function output = connectionWeight(dtiList, lengthList)

    for dtiIndex = 1:length(dtiList)
        try
            dtiSplit = strsplit(char(dtiList(dtiIndex)), '_');
            dtiSubjID = strjoin(dtiSplit(1:4), '_');

            dtiFileName = strcat(char(dtiSubjID), '_bedpostX_det_connectivity.csv');
            lengthsFileName = strcat(char(dtiSubjID), '_bedpostX_det_length.csv');
            voxelCountsFileName = strcat(char(dtiSubjID), '_ROIVoxelCount.csv');

            conMat = csvread(dtiFileName, 1, 0);
            lengthsMat = csvread(lengthsFileName, 1, 0);
            voxelCounts = csvread(voxelCountsFileName);
            weightedMat = [];
            weightedLength = [];

            for roi1 = 1:268
                for roi2 = 1:268
                    n_fibers = conMat(roi1, roi2);
                    weightedMat(roi1, roi2) = 2 * n_fibers / (voxelCounts(roi1) + voxelCounts(roi2));

                    if lengthsMat(roi1, roi2) > 0
                        weightedLength(roi1, roi2) = 2 / (voxelCounts(roi1) + voxelCounts(roi2)) / lengthsMat(roi1, roi2);
                    end
                end
            end

            weightedStreamName = strcat(char(dtiSubjID), '_bedpostX_det_connectivity_weighted.csv');
            weightedLengthName = strcat(char(dtiSubjID), '_bedpostX_det_length_weighted.csv');
            conmatPath = strcat('/scratch/lliu/SPINS/pipelines/conmat/', char(dtiSubjID));
            csvwrite(fullfile(conmatPath, weightedStreamName), weightedMat, 1, 0);
         	csvwrite(fullfile(conmatPath, weightedLengthName), weightedLength, 1, 0);

        catch
            continue
        end
    end
end
