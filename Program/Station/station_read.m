%% read the data from NOAA station and return daily-averaged variable
%  the missing data is indicated with NaN

tic;
clear;
path='../data/data_station/NOAA/';
dirlist=dir(strcat(path,'*.txt'));
begin_line=110;
time_begin=datenum('2000-01-01');
time_end=datenum('2021-12-31');
location=zeros(length(dirlist),2);
name=cell(length(dirlist),1);
for i = 1:length(dirlist)    
    standard_data=zeros(time_end-time_begin+1,9); 
    num_data=zeros(time_end-time_begin+1,9);
    file_name=dirlist(i).name;
    disp(strcat('processing file: ',file_name," ",num2str(i)));
    tem_path=strcat(path,file_name);
    file_content=importdata(tem_path,'\t',begin_line-1);
    variable_title=file_content.textdata{begin_line-1,2:end};
    tem_name = strsplit(file_content.textdata{7},{': '});
    name(i)=tem_name(2);
    position = strsplit(file_content.textdata{8},{' ',','});
    location(i,1) = str2double(position{3}(1:end-1));
    location(i,2) = str2double(position{4}(1:end-1));
    if position{3}(end)=='S'
        location(i,1) = -location(i,1);
    end
    if position{4}(end)=='W'
        location(i,2) = 360 - location(i,2);
    end

    tem_data=file_content.data;
    tem_shape=size(tem_data);
    data=nan(tem_shape(1),9);
    data(:,1:tem_shape(2))=tem_data;
    tem_data=data;
    data(isnan(data))=0;
    
    for j = 1:tem_shape(1)
        time_current=file_content.textdata{begin_line+j-1,1};
        index_time=datenum(time_current(1:10)) - time_begin + 1;
        standard_data(index_time,:)=standard_data(index_time,:)+data(j,:);
        num_data(index_time,:) = num_data(index_time,:) + ~isnan(tem_data(j,:));    
    end
    standard_data=standard_data./num_data;
    save(strcat('../data/data_station/analysis/station/',file_name,'.mat'),'standard_data');    
    disp('completed!')

end
save('../data/data_station/analysis/station/location.mat','location');
save('../data/data_station/analysis/station/dirlist.mat','dirlist');
save('../data/data_station/analysis/station/name.mat','name');
toc;