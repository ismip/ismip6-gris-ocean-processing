
%Need TEOS Toolbox!
%Also, requires interpBedmachineGreenland

disp('   -- Making ISMIP6 Grid');
[X Y]=meshgrid(-720e3:1e3:960e3,-3450e3:1e3:-570e3);

%Create connectivities

	disp('   -- Interpolating BedMachine bed and mask');
	BED = interpBedmachineGreenland(X,Y,'bed');
	M   = interpBedmachineGreenland(X,Y,'mask');
	disp('   -- Updating with RTOPO2');
	pos = find(BED==-9999);
	BED(pos) = interpRTopo2(X(pos),Y(pos));
	pos = find(M==-9999 & BED>0);
	M(pos) = 4;
	pos = find(M==-9999 & BED<=0);
	M(pos) = 0;

	disp('   -- Loading MIROC');
	load ocean_extrap_MIROC5_RCP85.mat

	disp('   -- Define open ocean pixels');
	FAR =(M~=0);
	FRAME = zeros(size(FAR));
	FRAME(1,:) = 1; FRAME(end,:) = 1; FRAME(:,1) = 1; FRAME(:,end) = 1;
	FAR = (1e3*bwdist(FAR) > 200e3 | FRAME);

	%Initialize output
	connectivity = false(size(X,1),size(X,2),numel(z));

	for i=1:numel(z),
		disp(['   -- Depth = ' num2str(z(i)) ' ' num2str(i) '/' num2str(numel(z))]);
		depth = z(i);
		A=(BED<=depth);
		disp('    ... Labelling');
		CC=bwlabel(A);
		disp('    ... Finding unique labels');
		pos=find(FAR & BED<depth);
		list = unique(CC(pos));
		connectivity(:,:,i)=ismember(CC,list);
	end

	save connectivity connectivity

%Compute Deepest index
	disp('   -- Loading MIROC');
	load ocean_extrap_MIROC5_RCP85.mat

	disp('   -- Initializing ID matrix');
	ID=int8(zeros(size(X)));

	for i=1:numel(z),
		disp(['   -- Depth = ' num2str(z(i)) ' ' num2str(i) '/' num2str(numel(z))]);
		pos=find(connectivity(:,:,i));
		ID(pos) = i;
	end

	save ID ID

%Basin Ids
	disp('   -- Initializing Basin ID matrix');
	BasinID=single(NaN(size(X)));

	disp('   -- Define open ocean pixels');
	FAR =(M~=0);
	FRAME = zeros(size(FAR));
	FRAME(1,:) = 1; FRAME(end,:) = 1; FRAME(:,1) = 1; FRAME(:,end) = 1;
	FAR = (150*bwdist(FAR) > 100e3 | (FRAME & M==0) );

	for i=1:numel(basins)
		disp(['   -- Basin ' num2str(i) ]);
		FLAGS=ContourTest(basins(i).X(1:end-1),basins(i).Y(1:end-1),X,Y);
		if i==4,
			A=((BED<0 & M>0) | (M==0 & FLAGS))  & Y<-1.8e6;
		elseif i==7,
			A=((BED<0 & M>0) | (M==0 & FLAGS))  & Y>=-1.8e6;
		else
			A=((BED<0 & M>0) | (M==0 & FLAGS));
		end
		disp('    ... Labelling');
		CC=bwlabel(A);
		disp('    ... Finding unique labels');
		pos=find(FAR & FLAGS);
		list = unique(CC(pos));
		pos=find(ismember(CC,list));
		BasinID(pos) = i;
	end

	%BasinID = fillmissing(BasinID,'nearest');
	BasinID = int8(BasinID);

	save BasinID BasinID

%Create Forcings
	%File Name
	filename = 'ocean_extrap_MIROC5_RCP85.mat';
	filename = 'ocean_extrap_ACCESS_RCP85.mat';
	filename = 'ocean_extrap_CSIRO_RCP85.mat';
	%filename = 'ocean_extrap_HadGEM_RCP85.mat';
	%filename = 'ocean_extrap_IPSLCM_RCP85.mat';
	%filename = 'ocean_extrap_MIROC5_RCP26.mat';
	%filename = 'ocean_extrap_MIROC5_RCP85.mat';
	%filename = 'ocean_extrap_NorESM_RCP85.mat';

	disp(['   -- Loading ' filename]);
	load(filename);

	disp('   -- Initializing Temperature, Salinity and depth matrices');
	temperature =single(NaN([size(X),numel(year)]));
	salinity    =single(NaN([size(X),numel(year)]));

	%Extrapolate
	for i=1:numel(basins)
		disp(['   -- Basin ' num2str(i)]);
		[posi posj]=find(BasinID==i & ID~=0);
		ind = sub2ind(size(BasinID),posi,posj);
		pos=ID(ind);
		for t=1:numel(year),
			ind = sub2ind(size(temperature),posi,posj,t*ones(size(posi)));
			temperature(ind) = T(i,t,pos);
			salinity(ind)    = S(i,t,pos);
		end
	end

	%Get depth from BedMachine
	depth  =interpBedmachineGreenland(X,Y,'bed');
	pos = find(depth==-9999);
	depth(pos) = interpRTopo2(X(pos),Y(pos));
	% since I only gave Mathieu temperatures to 2000 m, cap depth here
	depth(find(depth<-2000)) = -2000;

	disp('Converting T and S to TF');
	% pressure from depth
	% varying latitude from 60 to 80 results in 0.1% change in pressure
	lat0 = 70; % assume constant latitude = 70
	p = gsw_p_from_z(depth,lat0);

	% mask pressure by NaNs from ocean props
	NaNmask = squeeze(salinity(:,:,1))./squeeze(salinity(:,:,1));
	p = p.*NaNmask;

	% set any below 0 values to 0 (though after Mathieu's fix, there are none)
	p(find(p<0)) = 0;

	% absolute salinity from practical salinity
	% assuming constant lat-lon results in error on order of 0.01%
	long0 = -45;
	for i=1:size(salinity,3),
		disp(['Converting practical to absolute salinity: ',num2str(year(i))]);
		SA(:,:,i) = gsw_SA_from_SP(salinity(:,:,i),p,long0,lat0);
	end

	% in-situ temperature from potential temperature
	for i=1:size(salinity,3),
		disp(['Converting potential to in-situ temperature: ',num2str(year(i))]);
		temperature_insitu(:,:,i) = gsw_t_from_pt0(SA(:,:,i),temperature(:,:,i),p);
	end

	% in-situ freezing point
	for i=1:size(salinity,3),
		disp(['Calculating in-situ freezing point: ',num2str(year(i))]);
		temperature_freezing(:,:,i) = gsw_t_freezing(SA(:,:,i),p);
	end

	% thermal forcing
	TF = temperature_insitu - temperature_freezing;

	%Save as netcdf
	outputfile = strrep(filename,'.mat','_forcing.nc');

	% 1. Create netCDF file handle and Attributes
	%Checker: https://urldefense.com/v3/__http://pumatest.nerc.ac.uk/cgi-bin/cf-checker.pl__;!!Mih3wA!Q5BC0_XH4PIGfYmklBueGlV7ZhaZjuN3AeyoY1dSIK1JxkyRPhT3P5S3xLE9uuE$ 
	mode = netcdf.getConstant('NETCDF4');
	mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
	ncid=netcdf.create(outputfile,mode);
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title',['Ocean thermal forcing. Prepared for ISMIP6 (contact: mmorligh@uci.edu) ' date()]);
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Donald Slater, Mathieu Morlighem, Denis Felikson');
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'version',date());
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'nx',size(X,2));
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'ny',size(X,1));
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Projection','Polar Stereographic North (70N, 45W)');
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'proj4','+init=epsg:3413');
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'xmin',min(X(:)));
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'ymax',max(Y(:)));
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'no_data',NaN);
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'spacing',1000.);

	% 2. Define dimensions
	% Define x 
	x_id     = netcdf.defDim(ncid,'x',size(X,2));
	x_var_id = netcdf.defVar(ncid,'x','NC_FLOAT',x_id);
	netcdf.putAtt(ncid,x_var_id,'long_name',    'Cartesian x-coordinate');
	netcdf.putAtt(ncid,x_var_id,'standard_name','projection_x_coordinate');
	netcdf.putAtt(ncid,x_var_id,'axis','x');
	netcdf.putAtt(ncid,x_var_id,'units',        'm');

	% Define y
	y_id     = netcdf.defDim(ncid,'y',size(X,1));
	y_var_id = netcdf.defVar(ncid,'y','NC_FLOAT',y_id);
	netcdf.putAtt(ncid,y_var_id,'long_name',    'Cartesian y-coordinate');
	netcdf.putAtt(ncid,y_var_id,'standard_name','projection_y_coordinate');
	netcdf.putAtt(ncid,y_var_id,'axis','y');
	netcdf.putAtt(ncid,y_var_id,'units',        'm');

	%define time
	time_id     = netcdf.defDim(ncid,'time',numel(year));
	time_var_id = netcdf.defVar(ncid,'time','NC_int',time_id);
	netcdf.putAtt(ncid,time_var_id,'standard_name','time');
	netcdf.putAtt(ncid,time_var_id,'calendar','gregorian');
	netcdf.putAtt(ncid,time_var_id,'axis','T');
	netcdf.putAtt(ncid,time_var_id,'units','days since 1900-1-1 00:00:00');
	netcdf.putAtt(ncid,time_var_id,'bounds','time_bounds');
	two_id = netcdf.defDim(ncid,'two',2);
	timebounds_var_id = netcdf.defVar(ncid,'time_bounds','NC_int',[two_id time_id]);

	% Define temp
	temp_var_id = netcdf.defVar(ncid,'thermal_forcing','NC_FLOAT',[x_id y_id time_id]);
	netcdf.putAtt(ncid,temp_var_id,'long_name','ocean thermal forcing at effective depth');
	netcdf.putAtt(ncid,temp_var_id,'units',    'degree_C');
	netcdf.putAtt(ncid,temp_var_id,'grid_mapping', 'polar_stereographic');

	% Define mapping variable
	ps_var_id = netcdf.defVar(ncid,'polar_stereographic','NC_BYTE',[]);
	netcdf.putAtt(ncid,ps_var_id,'grid_mapping_name','polar_stereographic');
	netcdf.putAtt(ncid,ps_var_id,'latitude_of_projection_origin',90.);
	netcdf.putAtt(ncid,ps_var_id,'standard_parallel',70.);
	netcdf.putAtt(ncid,ps_var_id,'straight_vertical_longitude_from_pole',-45.);
	netcdf.putAtt(ncid,ps_var_id,'semi_major_axis',6378137.);
	netcdf.putAtt(ncid,ps_var_id,'inverse_flattening',298.257223563);
	netcdf.putAtt(ncid,ps_var_id,'false_easting',0.);
	netcdf.putAtt(ncid,ps_var_id,'false_northing',0.);

	% 3. compress data and end definition (level 1 to 10)
	netcdf.defVarDeflate(ncid,temp_var_id,true,true,2);
	netcdf.endDef(ncid);

	% 4. Place data
	year = year(:);
	time = datenum(year,7,1) - datenum(1900,1,1);
	timebounds = [datenum(year,1,1)-datenum(1900,1,1)  datenum(year+1,1,1) - datenum(1900,1,1)];
	netcdf.putVar(ncid,x_var_id,single(X(1,:)));
	netcdf.putVar(ncid,y_var_id,single(Y(:,1)));
	netcdf.putVar(ncid,time_var_id,int32(time));
	netcdf.putVar(ncid,timebounds_var_id,int32(timebounds'));
	netcdf.putVar(ncid,temp_var_id,single(permute(TF,[2,1,3])));

	% 5. Close file and display headers
	netcdf.close(ncid)
	system(['ncdump -h ' outputfile]);
