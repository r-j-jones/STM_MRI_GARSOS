function files = dirdicom(targetdir)
% dirdicom lists all DICOM files in a directory
%
% Usage: files = dirdicom(targetdir)
%    targetdir - char array with directory name
%    files - same structure as from MATLAB dir command
%
% Note: If the length of a path to a file exceeds 255 characters the path
% should be prefixed by "\\?\". For example: "c:\some\where" is changed to
% "\\?\c:\some\where". This allows dirdicom to list these files as well.
% Not using "\\?\" may cause problems when scanning DICOMDIR directories
% which have long directory and file names.

files = dir(targetdir);
files = files(~strcmp({files(:).name},{'..'})); % isdir does not work for .. if path is prefixed by \\?\
files = files(~strcmp({files(:).name},{'.'})); % isdir does not work for . if path is prefixed by \\?\
files = files(~[files.isdir]);

dicommask = false(size(files));

%pg = Progressor('Scanning files.');

for i = 1:numel(files)
    filename = fullfile(targetdir,files(i).name);
    try
                fid = fopen(filename);
                preamble = fread(fid,132,'uint8');
%                 if numel(preamble) == 132 && all(preamble(1:128) == 0) && all(char(preamble(129:132)') == 'DICM')
%                     dicommask(i) = true;
%                 end
                if (numel(preamble) == 132 && all(char(preamble(129:132)') == 'DICM'))
                    dicommask(i) = true;
                elseif numel(preamble) >= 2 
                    groupNumber = typecast(uint8(preamble(1:2)),'uint16');
                    dicommask(i) = (groupNumber == uint16(2) || groupNumber == uint16(8));
                end
                    
                fclose(fid);
                
%         dicommask(i) = isdicom(filename);
    catch
        fprintf('Skipping %s\n', filename);
    end
%     switch exist(filename,'file')
%         case 0
%             fprintf('Skipping %s\n', filename);
%         case 7
%             %Do nothing for folders. Needed for paths width \\?\ because
%             %isdir above does not work.
%         otherwise
%             fid = fopen(filename);
%             preamble = fread(fid,132,'uint8');
%             if numel(preamble) == 132 && all(preamble(1:128) == 0) && all(char(preamble(129:132)') == 'DICM')
%                 dicommask(i) = true;
%             end
%             fclose(fid);
%     end
    %pg.setProgress(i/numel(files));
end

files = files(dicommask);
