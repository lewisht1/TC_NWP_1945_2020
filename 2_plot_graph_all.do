
/******************************************************************************
	This do file perform the following:
	1. Import the stata dta files converted from shape files of the world map
	2. Plot the world map and add all TC tracks on top 
	
	Need: stata dta files converted from shape files of the world map (from
		  external source), master.dta
	Make: Graph of all TC tracks from 1945-2020 Northwest Pacific Basin
*******************************************************************************/

* Initialising interface
set more off
graph drop _all
clear
cls

quietly {

	* import the dta files converted from shape files of the world map (from external source)
	cd "directory 3"

	use world_shp.dta
	
	*merge shape file with file contained geographic ids
	merge m:1 _ID using world 

	rename _Y lat
	rename _X lon
	
	*remove Antarctica from the map
	drop if COUNTRY == "Antarctica"
	drop rec_header-_merge
	
	*shift the focus of the world map to Pacific Ocean
	replace lon = 360 + lon if lon < 0
	
	*restrict the latitude and longitude of the map 
	keep if inrange(lat, 0, 60)
	keep if inrange(lon, 100, 180)

	set graphics off

	*plot the world map
	scatter lat lon, aspect(`=60/80') xla(100(5)180, grid) yla(0(5)60, ang(h) grid) yti(, orient(horiz)) plotregion(margin(none)) msize(tiny) msymbol(point) play("map.grec")

	* plot the TC tracks on top of the world map
	cd "directory 4"

	use master.dta, clear
	
	*shift the focus of the world map to Pacific Ocean
	replace lon = 360 + lon if lon < 0
	
	*restrict the latitude and longitude of the map 
	keep if inrange(lat, 0, 60)
	keep if inrange(lon, 100, 180)
	
	*combine TCs orginated from different basins and recreate unique TC id
	replace basin = "1" if basin == "WP"
	replace basin = "2" if basin == "CP"
	destring basin, replace

	sort date_min storm_id row_no

	generate year2 = year
	replace year2 = year2 - 1 if storm_number == .

	generate uid = 1000000 * basin + 100 * year2 + storm_id

	*create a leading row header for each TC track (for saving time of plotting)
	quiet levelsof uid, local(levels)
	local n_r = r(r)
	mat a = J(`n_r', 4, .)
	 
	local c = 1
	foreach i of local levels {
		
		mat a[`c', 1] = `i'
		mat a[`c', 2] = real(substr("`i'", 1, 1))
		mat a[`c', 3] = real(substr("`i'", 2, 4))
		mat a[`c', 4] = real(substr("`i'", -2, 2))
		
		local c = `c' + 1
		
	}

	local n1 = _N + 1
	local obs = _N + `n_r'
	set obs `obs'


	forvalues i = `n1'/`obs' {
		
		replace row_no = 0 in `i'
		replace uid = a[`i' - `n1' - 1, 1] in `i'
		
	}
	 
	sort uid row_no 

	* plot the TC tracks on top of the map
	addplot: line lat lon, xla(100(5)180, grid) yla(0(5)60) lwidth(0.001cm) lcolor(red%50) legend(off) cmissing(n)

	set graphics on
	graph display Graph
	
	* save graph
	noisily graph export "directory 3\Graph.pdf", as(pdf) name("Graph") replace

}
