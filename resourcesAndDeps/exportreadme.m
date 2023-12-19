function exportreadme(fformat)
    % EXPORTREADME - Create HTML version of readme with exact GitHub formating.
    if nargin < 1
        fformat = "";
    end
    setenv('PATH', getenv('PATH')+":/opt/homebrew/bin:/usr/local/bin")
    cmdHTML = 'grip --export readme.md';
    cmdPDF = "pandoc -s -f gfm -t html5 --metadata-file resourcesAndDeps/readme_format.yaml --metadata title="""" --pdf-engine wkhtmltopdf readme.md -o readme.pdf";
    switch fformat
        case "html"
            system(cmdHTML);
        case "pdf"
            system(cmdPDF); 
        otherwise
            system(cmdHTML);
            system(cmdPDF); 
    end
end
