function pass = userTestConfirmation(fig)
    answer = questdlg('Is the plot good?','confirm test', 'yes', 'no', 'no');
    close(fig);
    pass = false;
    if strcmp(answer, 'yes')
        pass = true; 
    end
end
