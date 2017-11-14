% in bash
% module load freesurfer
% module load SPM

% this counts the number of voxels within each ROI of the 268 Shen atlas

function output = voxelCount(dtiList, fmriList)

    fmriSubjIDList = {};
    for i = 1:length(fmriList)
        try
        fmriSplit = strsplit(char(fmriList(i)), '_');
        fmriSubjID = cellstr(strjoin(fmriSplit(1:4), '_'))

        niiFile = char(gunzip(strcat(fmriSubjID, '_registered_shen.nii.gz')));
        headerInfo = spm_vol(niiFile);
        data = spm_read_vols(headerInfo);
        voxels = [];

        for i = 1:268
            voxels(i) = length(find(data == i));
        end

        csvName = strcat(char(fmriSubjID), '_ROIVoxelCount.csv');
        csvwrite(csvName, voxels);
        catch
            continue
        end
    end

end
