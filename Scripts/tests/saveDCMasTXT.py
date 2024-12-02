import pydicom as pyd

dfname = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\Imaging Data\20241128_112522_RB_Q1_RB_Q1_1_1\36\pdata\1\dicom\MRIm001.dcm";

dcm = pyd.dcmread(dfname)

with open('dcm.txt', 'w') as f:
    f.write(str(dcm))