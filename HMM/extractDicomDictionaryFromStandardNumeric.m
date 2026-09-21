function [D61, D71, DLN] = extractDicomDictionaryFromStandardNumeric()

persistent pd;

% xmlFile - PS3.6 in the DICOM standard in xml DocBook format. Usually called "part06.xml".

if isempty(pd)    
    %xmlFile = fullfile(fileparts(mfilename('fullpath')),'part06.xml');
    %xmlFile = fullfile(fileparts(mfilename('fullpath')),'dicom_part06_2017c.xml');
    xmlFile = which('part06.xml');
    %xmlFile = which('dicom_part06_2017c.xml');
    
    xmlDoc = xmlread(xmlFile);
    
    % Debugging of new dictionary
    %documentSubtitle = extractSubtitle(xmlDoc);    
    %fprintf('Loading DICOM data dictionary: "%s"\n', documentSubtitle);
    
    T61 = extractElementListFramTable(xmlDoc, '6-1'); % Extract rows from "Table 6-1. Registry of DICOM Data Elements".
    T71 = extractElementListFramTable(xmlDoc, '7-1'); % Extract rows from "Table 7-1. Registry of DICOM File Meta Elements".
    
    %m61 = cellfun(@isempty,strfind(T61(:,1),'x'));
    m61 = ~contains(T61(:,1),'x');
    %m71 = cellfun(@isempty,strfind(T71(:,1),'x'));
    m71 = ~contains(T71(:,1),'x');
    
    T61 = T61(m61,:);
    T71 = T71(m71,:);
    
    T61(:,1) = cellfun(@dicomHexPairToNumber,T61(:,1),'UniformOutput',false);
    T71(:,1) = cellfun(@dicomHexPairToNumber,T71(:,1),'UniformOutput',false);
    
    T61(:,3) = cellfun(@(x)regexprep(x,'[^a-zA-Z0-9]',''),T61(:,3),'UniformOutput',false);
    T71(:,3) = cellfun(@(x)regexprep(x,'[^a-zA-Z0-9]',''),T71(:,3),'UniformOutput',false);
    
    m61 = cellfun(@(x)~isempty(strtrim(x)),T61(:,3));
    m71 = cellfun(@(x)~isempty(strtrim(x)),T71(:,3));
    
    T61 = T61(m61,:);
    T71 = T71(m71,:);
    
    pd.D61 = containers.Map(T61(:,1),num2cell(T61(:,2:end),2));    
    pd.D71 = containers.Map(T71(:,1),num2cell(T71(:,2:end),2));
    pd.DLN = containers.Map(cat(1, T71(:,3), T61(:,3)), cat(1, T71(:,1), T61(:,1)));
end

D61 = pd.D61; % Dictionary for Table 6-1
D71 = pd.D71; % Dictionary for Table 7-1
DLN = pd.DLN; % Dictionary to Tag(e.g. 1048592) from Keyword(e.g. PatientName)

function documentSubtitle = extractSubtitle(doc)

xPath = javax.xml.xpath.XPathFactory.newInstance().newXPath();

expressionTitle = xPath.compile('/*:book/*:subtitle');
subtitles = expressionTitle.evaluate(doc, javax.xml.xpath.XPathConstants.NODESET);
documentSubtitle = char(subtitles.item(0).getTextContent);

function d = extractElementListFramTable(doc, tableLabel)

xPath = javax.xml.xpath.XPathFactory.newInstance().newXPath();

expressionRows = xPath.compile(['//*:table[@label=''' ,tableLabel, ''']/*:tbody/*:tr']);
expressionColumns = xPath.compile('.//*:td/*:para');
tableRows = expressionRows.evaluate(doc,javax.xml.xpath.XPathConstants.NODESET);

d = cell(tableRows.getLength(), 6);

for iDicomTag = 1:tableRows.getLength()
    p = expressionColumns.evaluate(tableRows.item(iDicomTag - 1),javax.xml.xpath.XPathConstants.NODESET);
    d(iDicomTag,1:p.getLength()) = arrayfun(@(x)strtrim(char(p.item(x - 1).getTextContent)), 1:p.getLength(),'UniformOutput', false);
end

d(~cellfun(@ischar,d)) = {''};