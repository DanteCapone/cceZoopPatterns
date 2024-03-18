function [shippos,temperature,salinity,fluorescence,flowmeter, oxygen] = GetMetData(monthdaystart,monthdayend)

% clearvars
% close all

ncolumns=101;
latpos = 67;
lonpos = 68;
salpos = 28;
temppos = 33;
fluorpos = 44;   %FL
timepos = 1;
oxypos=43; %oxygen
%flowposG =    %FM = USW Flow meter (Gallons per minute)
flowposL = 32;   %FI = USW Flow meter (liters per minute)

folder_met = 'C:\Users\DanteACapone\Desktop\Summer_2022\CCE_REU\Decima_Matthews_Gallego_Capone\data\unmerged_metadata\underway_data\met\';


folder = 'tempmetdata\';

curfolder = cd;
folder = [curfolder,'\',folder];

copyfile([folder_met,'2*.MET'], folder)

files = dir(folder);


% filename = '210715.MET'

counter=1;
for i=1:length(files(:,1))
    if length(files(i,:).name)>3
        if strcmp(files(i,:).name(end-2:end),'MET')
            filename=files(i,:).name;
            if str2num(filename(4:6))>=monthdaystart & str2num(filename(4:6))<=monthdayend
                filename
                fid = fopen([folder,filename]);
                header=textscan(fid,'%s %*[^\n]',3);
                header=textscan(fid,repmat('%s ',[1,ncolumns]),1);
                data=textscan(fid,repmat('%f ',[1,ncolumns]));
                lat = data{1,latpos};
                lon = data{1,lonpos};
                tm = data{1,timepos};
                sal = data{1,salpos};
                temp = data{1,temppos};
                fluor = data{1,fluorpos};
                flowL = data{1,flowposL};
                oxy=data{1,oxypos};
                
                if counter==1
                    shippos = [lon,lat];
                    time=tm;
                    temperature = temp;
                    salinity = sal;
                    fluorescence = fluor;
                    flowmeter = flowL;
                    oxygen = oxy;
                else
                    shippos = [shippos;lon,lat];
                    time = [time;tm];
                    temperature = [temperature;temp];
                    salinity = [salinity;sal];
                    fluorescence = [fluorescence;fluor];
                    flowmeter = [flowmeter;flowL];
                    oxygen = [flowmeter;oxy];
                end
                counter = counter + 1;
                
                fclose(fid);
            end
        end
    end
end

temperature(find(flowmeter<0.5))=NaN;
salinity(find(flowmeter<0.5))=NaN;
fluorescence(find(flowmeter<0.5))=NaN;


% figure
% scatter(lon,lat,50,temp);
% colorbar