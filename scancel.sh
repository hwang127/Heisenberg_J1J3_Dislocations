#!/bin/csh
foreach id (`seq 37248858 37249199`) 
	scancel ${id}
end

