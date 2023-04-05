function S = save(obj,file)
    % save - Save Material object into a json-formatted text file.
    %
    % Arguments:
    % - obj:   Material object.
    % - file:  (optional) file name as a string or char array. Default: GEWtool's database folder.
    %
    % Return value:
    % - s:     The json-formatted string that has been saved.
    %
    % If "file" exists already, a window will prompt the user where to save instead 
    % (re-select to overwrite). The saved properties are "name", "rho" and
    % either "C" (anisotropic) or "lambda" and "mu" (isotropic).
    % 
    % 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    % determine where to save:
    if nargin < 2
%         path = fileparts(which('aluminum.json'));
        path = '.';
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
    S = jsonencode(mat, "PrettyPrint", true);
    try
        fid = fopen(file, 'w'); fprintf(fid, S); fclose(fid); 
    catch e
        fclose(fid); warning('Could not save file'); rethrow(e);
    end
    fprintf('Saved to %s.\n', file);
end
