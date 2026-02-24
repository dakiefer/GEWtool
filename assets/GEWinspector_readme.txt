How to create the gif: 

Record using shift-command-5. This creates the file "GEWinspector.mov". 

Convert mov to gif using ffmpeg: 
>> ffmpeg -i GEWinspector.mov -pix_fmt rgb24 -r 10 -s 880x428 GEWinspector.gif

options are: 
-r: frame rate. 
-s: size
-pix_fmt: the pixel format
