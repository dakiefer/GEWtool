function S = save(obj,file)
    % determine where to save:
    if nargin < 2
        path = fileparts(which('aluminum.json'));
        file = fullfile(path, [obj.name '.json']);
    end
    if exist(file, 'file')==2
        [file,path] = uiputfile([obj.name '.json'], 'File already exists: select where to save.');
        if file==0, warning('Cancelled save operation.'); return; end
        file = fullfile(path, file);
    end
    % extract relevant data:
    mat.name = obj.name;
    mat.rho = obj.rho;
    if strcmp(obj.symmetry, 'isotropic') 
        mat.lambda = obj.lambda;
        mat.mu = obj.mu;
    else
        mat.C = obj.C;
    end
    mat.symmetry = obj.symmetry;
    S = jsonencode(mat);
    try
        fid = fopen(file, 'w'); fprintf(fid, S); fclose(fid); 
    catch e
        fclose(fid); warning('Could not save file'); rethrow(e);
    end
    fprintf('Saved to %s.\n', file);
end
