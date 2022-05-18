/******************************************************************************
	This do file perform the following:
	1. Import raw data files (one season per file) with the best tracks 
	   (position, intensity etc.) of tropical cyclones (TCs) e.g., tropical 
	   storm, typhoon, hurricane from the Northwest Pacific between 1945 to 2020 
	   downloaded from the Joint Tyhoon Warning Center (JTWC)
	2. Calculate distances and bearing from Hong Kong (HK), moving speed, 
	   approaching bearing and lanfall info. etc.
	3. Harmonize, combine and save files to stata readable (dta) files
	
	Need: all raw data files stored in directory 1
	Make: stata dta files saved in directory 2 and a master file containing all
		  all data from 1945-2020
*******************************************************************************/

* Initialising interface
set more off
graph drop _all
cls

quietly {

	* store folder directories
	local folder_raw "directory 1"

	local folder_stata "directory 2"

	* input latitude and longitude (decimal degree) of HK
	local hk_lat = 22 + 18 / 60 + 7 / 3600
	local hk_lon = 114 + 10 / 60 + 27 / 3600

	cd "`folder_raw'"

	* import raw data files strored from directory 1
	filelist, dir("`folder_raw'")
	keep filename

	local f_n = _N

	* loop all files
	forvalues f = 1/`f_n' {

			local fname = filename[`f']
			
			preserve
				import delimited `fname', varnames(nonames) clear
				
				quiet describe
				local var_no = r(k)
				quiet ds
				local vlist "`r(varlist)'"
				
				* harmonise names for each column
				forvalues i = 1/`var_no' {
				
					local var_name `:word `i' of `vlist''
					local check_name "v`i'"
					if "`var_name'" != "`check_name'" {
					
						rename `var_name' `check_name'
					
					}
				
				}
							
				keep v1 v2 v3 v7 v8 v9
				
				* remove duplicate rows
				if _N > 1 {
					
					sort v2 v3 v7 v8 v9
					quietly by v2 v3 v7 v8 v9: gen dup = cond(_N == 1, 0, _n)
					drop if dup > 1
					
				}
				
				egen v_max = max(v9)
				
				* recode latitude and longitude
				generate lat double = real(substr(v7, 1, length(v7) - 1)) / 10
				replace lat = 0 - lat if substr(v7, -1, 1) == "S"			
				generate lon double = real(substr(v8, 1, length(v8) - 1)) / 10
				replace lon = 0 - lon if substr(v8, -1, 1) == "W"
				
				* only include north hemispheric data
				drop if lat < 0
				
				* recode date
				generate year = real(substr(string(v3, "%12.0g"), 1, 4))
				generate mon = real(substr(string(v3, "%12.0g"), 5, 2))
				generate day = real(substr(string(v3, "%12.0g"), 7, 2))		
				generate date = mdy(mon, day, year)
				generate hr = real(substr(string(v3, "%12.0g"), 9, 2))
				format date %td
				
				generate day2event = date - mdy(1, 1, year) + 1
				
				* generate category of intensity of TCs
				generate tc_cat = .
				replace tc_cat = 1 if v9 < 35
				replace tc_cat = 2 if v9 >= 35 & v9 < 50
				replace tc_cat = 3 if v9 >= 50 & v9 < 65
				replace tc_cat = 4 if v9 >= 65 & v9 < 130
				replace tc_cat = 5 if v9 >= 130
				
				egen tc_max = max(tc_cat)
				
				* calculate distance (in km) of TC position from HK using Haversine formula
				generate dlon double = lon - `hk_lon'
				generate dlat double = lat - `hk_lat'
				generate a double = (sin(dlat*_pi/180/2))^2 + cos(`hk_lat'*_pi/180) * cos(lat*_pi/180) * (sin(dlon*_pi/180/2))^2
				generate c double = 2 * asin(min(1, sqrt(a)))
				generate dis_hk double = 6371 * c
				
				* generate categories of distances from Hong Kong
				generate hk_800km = 0
				replace hk_800km = 1 if dis_hk <= 800
				generate hk_500km = 0 
				replace hk_500km = 1 if dis_hk <= 500
				generate hk_200km = 0 
				replace hk_200km = 1 if dis_hk <= 200			
				
				* calculate inital bearing (in decimal degree) of TC position from HK
				generate term1 double = sin(dlon*_pi/180) * cos(lat*_pi/180)
				generate term2 double = cos(`hk_lat'*_pi/180) * sin(lat*_pi/180)
				generate term3 double = sin(`hk_lat'*_pi/180) * cos(lat*_pi/180) * cos(dlon*_pi/180)
				generate bearing_from_hk1 double = atan2(term1, term2 - term3) * 180/_pi
				replace bearing_from_hk1 = bearing_from_hk1 + 360 if bearing_from_hk1 < 0
				
				* calculate final bearing
				replace dlon = -dlon
				replace term1 = sin(dlon*_pi/180) * cos(`hk_lat'*_pi/180)
				replace term2 = cos(lat*_pi/180) * sin(`hk_lat'*_pi/180)
				replace term3 = sin(lat*_pi/180) * cos(`hk_lat'*_pi/180) * cos(dlon*_pi/180)
				generate bearing_from_hk2 double = atan2(term1, term2 - term3) * 180/_pi
				replace bearing_from_hk2 = bearing_from_hk2 + 180 
				replace bearing_from_hk2 = bearing_from_hk2 - 360 if bearing_from_hk2 >= 360
				replace bearing_from_hk2 = bearing_from_hk2 + 360 if bearing_from_hk2 < 0
				
				* average initial and final bearing
				generate bearing_from_hk double = (bearing_from_hk1 + bearing_from_hk2) / 2
				
				* calculate speed (km/h) and approaching direction (decimal degree) TC
				local r_n = _N
				generate speed double = .
				generate approach_bearing = .
				generate approach_bearing1 = .
				generate approach_bearing2 = .
				
				if _N > 1 {
				
					forvalues i = 2/`r_n' {
						
						* calculate time differences between each data point for each TC
						local lat1 = lat[`i'- 1]
						local lat2 = lat[`i']
						local lon1 = lon[`i'- 1]
						local lon2 = lon[`i']
						local dlon = `lon2' - `lon1'
						local dlat = `lat2' - `lat1'
						local t1 = hr[`i'- 1]
						local t2 = hr[`i']
						if `t1' == . | `t2' == . | day[`i'] == . | day[`i' - 1] == .{
						
							local dt = .
								
						}
						else {
						
							local dt = (24 * day[`i'] + `t2') - (24 * day[`i' - 1] + `t1')
						
						}
						
						* calculate speed
						local a = (sin(`dlat'*_pi/180/2))^2 + cos(`lat1'*_pi/180) * cos(`lat2'*_pi/180) * (sin(`dlon'*_pi/180/2))^2
						local c = 2 * asin(min(1, sqrt(`a')))										
						local speed = 6371 * `c' / `dt'
						replace speed = `speed' in `i'
						
						* calculate approach bearing
						*initial bearing
						local term1 = sin(`dlon'*_pi/180) * cos(`lat2'*_pi/180)
						local term2 = cos(`lat1'*_pi/180) * sin(`lat2'*_pi/180)
						local term3 = sin(`lat1'*_pi/180) * cos(`lat2'*_pi/180) * cos(`dlon'*_pi/180)
						local approach_bearing1 = atan2(`term1', `term2' - `term3') * 180/_pi
						replace approach_bearing1 = `approach_bearing1' in `i'
						replace approach_bearing1 = approach_bearing1 + 360 if approach_bearing1 < 0
						
						*final bearing 
						local dlon = -`dlon'
						local term1 = sin(`dlon'*_pi/180) * cos(`lat1'*_pi/180)
						local term2 = cos(`lat2'*_pi/180) * sin(`lat1'*_pi/180)
						local term3 = sin(`lat2'*_pi/180) * cos(`lat1'*_pi/180) * cos(`dlon'*_pi/180)
						local approach_bearing2 = atan2(`term1', `term2' - `term3') * 180/_pi
						local approach_bearing2 = `approach_bearing2' + 180
						
						*make sure bearing within 0 to 360 degrees
						if `approach_bearing2' >= 360 {
						
							local approach_bearing2 = `approach_bearing2' - 360
						
						} 
						else if `approach_bearing2' < 0 {
						
							local approach_bearing2 = `approach_bearing2' + 360
							
						} 
						else {
						
						}
						
						
						replace approach_bearing2 = `approach_bearing2' in `i'			
										
					}
					
					*average initial and final bearing 
					replace approach_bearing = (approach_bearing1 + approach_bearing2) / 2
				
				}
				
				replace bearing_from_hk = . if speed == 0
				replace approach_bearing = . if speed == 0
				
				* format file
				rename v1 basin
				rename v2 storm_id
				rename v9 v_kts
				
				if _N > 1 {
				
					drop v3 v7 v8 dup dlon-c term1-term3 bearing_from_hk1 bearing_from_hk2 approach_bearing1 approach_bearing2	
				
				}
				else {
				
					drop v3 v7 v8 dlon-c term1-term3 bearing_from_hk1 bearing_from_hk2 approach_bearing1 approach_bearing2		
					
				}

				* determine if TC is in land by joining landmass matrix data file (generated from other sources)
				merge m:1 lat lon using "directory 3\is_land.dta", keepusing(is_land)
				
				keep if storm_id != .
				
				drop _merge
				
				* sort data by TC id, date 
				sort date hr
				
				generate row_no = _n
				egen date_min = min(date)
							
				order basin storm_id year mon day hr date day2event lat lon is_land speed approach_bearing dis_hk hk_800km hk_500km hk_200km bearing_from_hk v_kts tc_cat v_max tc_max row_no date_min
				
				* create file names for dta files 
				sort date hr

				local f1 = year[1]
				local f2 = mon[1]
				if `f2' < 10 {
				
					local f2 "0`f2'"
				
				}
				local f3 = day[1]
				if `f3' < 10 {
				
					local f3 "0`f3'"
				
				}
				local f4 = basin[1]
				local f5 = storm_id[1]
				if `f5' < 10 {
				
					local f5 "0`f5'"
				
				}
				
				* save dta files (one file per season)
				noisily save "`folder_stata'\\`f1'_`f2'_`f3'_`f4'_`f5'_`f'.dta", replace
				
			restore
		
	}

}

* combine all dta files into a single master dta file
quietly {

	cd "`folder_stata'"
	clear
	append using `: dir . files "*.dta"'

	*remove and modify unrealistic entries
	drop if day == 0
	replace v_kts = 0 if v_kts < 0
	replace v_max = 0 if v_max < 0

	*generate TC unique identifier
	sort date_min storm_id row_no 
	bysort year: egen storm_total = total(row_no == 1)

	generate tmp = 0
	replace tmp = 1 if row_no == 1
	sort date_min storm_id row_no 
	bysort year (date_min row_no): gen storm_number = sum(tmp)
	replace storm_number = . if storm_number == 0
	drop tmp

	sort date_min storm_id row_no 

	*save master file
	noisily save "`folder_stata'\master_dtas\master.dta", replace

}
