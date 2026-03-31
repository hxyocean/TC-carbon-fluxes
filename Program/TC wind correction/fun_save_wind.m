%% this function is used to store the rebuilt wind data
function fun_save_wind(path_save, tem_current_time, wind_tc, wind_filt,...
    wind_orign, lon, lat)
    
    tem_current_time = datestr(datenum(tem_current_time,'yyyymmdd'), 'yyyymmdd');
    tem_file = fullfile(path_save);

    tem_file_name = ['CCMP_',tem_current_time, '_V2.0_TCWIND_REC.nc'];

    ncid = netcdf.create(fullfile(tem_file, tem_file_name),'NC_CLOBBER');
    dimid_lon = netcdf.defDim(ncid,'lon',length(lon)); 
    dimid_lat = netcdf.defDim(ncid,'lat',length(lat)); 
%     dimid_time = netcdf.defDim(ncid,'time',4); 
 
    m_type = 'NC_FLOAT';
    varid_wind_tc = netcdf.defVar(ncid,'wind_tc',m_type,[dimid_lon dimid_lat]);
    varid_wind_filt = netcdf.defVar(ncid,'wind_filt',m_type,[dimid_lon dimid_lat]);
    varid_wind_orign = netcdf.defVar(ncid,'wind_orign',m_type,[dimid_lon dimid_lat]);


    varid_lon = netcdf.defVar(ncid,'lon',m_type,[dimid_lon]);
    varid_lat = netcdf.defVar(ncid,'lat',m_type,[dimid_lat]);

    netcdf.endDef(ncid)
    netcdf.putVar(ncid,varid_wind_tc,wind_tc);
    netcdf.putVar(ncid,varid_wind_filt,wind_filt);
    netcdf.putVar(ncid,varid_wind_orign,wind_orign);
    netcdf.putVar(ncid,varid_lon,lon);
    netcdf.putVar(ncid,varid_lat,lat);
    netcdf.close(ncid);
end