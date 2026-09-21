function setAugmentedDicomDictionary()

%dictFile = '/home/ajohanss/prog/workspace/umich20150223/dicom-dict-with-private-tags.txt';

%dictFile = '/home/youdae/dicom-dict-with-private-tags.txt'
dictFile = which('dicom-dict-with-private-tags.txt');

dicomdict('set', dictFile);

reloadDicomDictionary();