# This can build the "chimeras"
# run: $ make
# 
# scripts:
#   convertphotsys.py converts the COSMOS galaxy sample into the photometry frame of SDSS quasars
#   build_chimeras.py  
#      loads the SDSS quasars', also converts the WISE bands into IRAC
#      for each COSMOS galaxy and SDSS quasar, builds all possible chimeras
#   resample.py
#      deduplicates based on similar colors and redshifts
  

all: chimeras-grahsp.fits chimeras-cigale.fits chimeras-fullinfo.fits

clean:
	rm -rf input-gal_fulldata.fits input-gal_converted.fits QSO-fluxes.fits chimeras-list.txt chimeras-*.fits  chimeras/*

.PHONY: clean all
.SECONDARY:
.SUFFIXES:

#input-gal_fulldata.fits: gal/input.fits
#	# attach information from COSMOS2015 on the best-fit galaxy model
#	stilts tmatch2 in1=gal/input.fits values1=id in2=gal/COSMOS2015_Laigle+.fits values2=NUMBER matcher=exact suffix1= suffix2=_COSMOS2015 fixcols=dups out=$@

gal/input.fits: 
	$(MAKE) -C gal

input-gal_converted.fits: convertphotsys.py gal/input.fits
	# convert photometry system from COSMOS to that of the quasar data set
	python3 $^ $@

# QSO.fits is the Shen et al. DR7 quasar catalog, downloaded as a fits file from vizier
QSO-orig.fits: 
	wget -O $@ https://cdsarc.cds.unistra.fr/viz-bin/nph-Cat/fits?J/ApJS/194/45/catalog.dat.gz

QSO-fluxes.fits: convert_QSOs.py QSO-orig.fits
	# build all possible chimera combinations
	python3 $^ $@

chimeras-list.txt: build_chimeras.py input-gal_converted.fits QSO-fluxes.fits
	# build all possible chimera combinations
	mkdir -p chimeras/
	python3 $^ $@

chimeras-full.fits: chimeras-list.txt
	# concatenate the individual files with stilts
	stilts tcatn out=$@ ocmd="sort RANDOM" $$( n=0; while read line; do n=$$((n+1)); echo in$${n}=$${line}; done < $^; echo nin=$$n; )

chimeras-deduplicated-ids.fits: dedup_chimeras.py chimeras-full.fits
	python3 $^ $@

chimeras-deduplicated.fits: chimeras-deduplicated-ids.fits
	stilts tmatch2 out=$@ matcher=exact \
		in1=chimeras-full.fits values1=id suffix1= \
		in2=chimeras-deduplicated-ids.fits values2=id suffix2=_dedup \
		ocmd='delcols "id_dedup"'

chimeras-fullQinfo.fits: chimeras-deduplicated.fits QSO-fluxes.fits
	stilts tmatch2 out=$@ matcher=exact find=best1 join=all1 fixcols=all \
		in1=chimeras-deduplicated.fits values1=ID_SDSS suffix1='' \
		in2=QSO-fluxes.fits values2=SDSS suffix2='_QSO' \
		ocmd='delcols "GroupSize GroupID"'

chimeras-fullinfo-base.fits: chimeras-fullQinfo.fits input-gal_converted.fits
	stilts tmatch2 out=$@ matcher=exact find=best1 join=all1 fixcols=all \
		in1=chimeras-fullQinfo.fits values1=ID_COSMOS suffix1= \
		in2=input-gal_converted.fits values2=id suffix2=_GAL \
		ocmd='delcols "GroupSize GroupID"'

chimeras-fullinfo.fits: add_host_contamination.py chimeras-fullinfo-base.fits
	# estimate the quasar host light added to each COSMOS galaxy
	python3 $^ $@

chimeras-grahsp.fits: chimeras-deduplicated.fits
	cp $^ $@

chimeras-cigale.fits: chimeras-grahsp.fits
	topcat -stilts tpipe in=$< out=$@ \
		cmd='colmeta -name "sdss.up" u_sdss' \
		cmd='colmeta -name "sdss.rp" r_sdss' \
		cmd='colmeta -name "sdss.ip" i_sdss' \
		cmd='colmeta -name "sdss.zp" z_sdss' \
		cmd='colmeta -name "sdss.up_err" u_sdss_err' \
		cmd='colmeta -name "sdss.rp_err" r_sdss_err' \
		cmd='colmeta -name "sdss.ip_err" i_sdss_err' \
		cmd='colmeta -name "sdss.zp_err" z_sdss_err' \
		cmd='colmeta -name "2mass.J" J_2mass' \
		cmd='colmeta -name "2mass.H" H_2mass' \
		cmd='colmeta -name "2mass.Ks" Ks_2mass' \
		cmd='colmeta -name "2mass.J_err" J_2mass_err' \
		cmd='colmeta -name "2mass.H_err" H_2mass_err' \
		cmd='colmeta -name "2mass.Ks_err" Ks_2mass_err' \
		cmd='sort redshift'
