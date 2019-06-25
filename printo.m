%Script that prints a tiff files from top graph
if ishghandle(gcf)
    print -r300 -dtiff fibrils.tif

end