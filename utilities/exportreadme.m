function exportreadme(fformat)
    % EXPORTREADME - Create HTML version of readme with exact GitHub formating.
    if nargin < 1
        fformat = "";
    end
    setenv('PATH', getenv('PATH')+":/opt/homebrew/bin:/usr/local/bin:/Library/TeX/texbin")
    cmdHTML = 'grip --export readme.md';
    % cmdPDF = "pandoc -s -f gfm -t html5 --metadata-file resourcesAndDeps/readme_format.yaml --metadata title="""" --pdf-engine wkhtmltopdf readme.md -o readme.pdf"; old export
    cmdPDF = "pandoc -s --pdf-engine xelatex -V geometry:""left=2.5cm,right=2.5cm,top=1.5cm,bottom=2cm"" readme.md -o readme.pdf";
    switch fformat
        case "html"
            system(cmdHTML);
        case "pdf"
            disp('Exporting to readme.pdf')
            system(cmdPDF); 
        otherwise
            system(cmdHTML);
            disp('Exporting to readme.pdf')
            system(cmdPDF); 
    end
end
