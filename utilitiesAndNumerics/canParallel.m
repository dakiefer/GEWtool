function canPar = canParallel()
% CANPARALLEL - Check if a code can be run on a parallel pool.
% Checks if the parallel computing toolbox is installed and whether a license 
% can be checked out.
% 
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

persistent hasPar   % save result for future checks

if isempty(hasPar)
    isInstalled = ~isempty(ver('parallel'));
    hasLicense = license('test','Distrib_Computing_Toolbox');
    hasPar = isInstalled && hasLicense;
end

canPar = hasPar;

end
