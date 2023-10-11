# This can build the "pure galaxy sample" from COSMOS
# run: $ make

all: input.fits input-cigale.fits

# find isolated galaxies with spec-z
COSMOS2015_Laigle+.fits_filtered_gal_zspec_all.fits: convert_gal.py test_mass_onlygal_extragal.fits cosmos_acs_iphot_200709.coords filternames.txt
	python3 convert_gal.py test_mass_onlygal_extragal.fits $@

# find AGN with spec-z
input-COSMOS-AGN-specz.fits: convert_gal.py COSMOS2015_Laigle+.fits chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits cosmos_acs_iphot_200709.coords filternames.txt
	python3 convert_agn.py COSMOS2015_Laigle+.fits chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits $@
input-COSMOS-AGN-photz.fits: convert_gal.py COSMOS2015_Laigle+.fits chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits cosmos_acs_iphot_200709.coords filternames.txt
	ONLY_PHOTOZ=1 python3 convert_agn.py COSMOS2015_Laigle+.fits chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits $@
input-COSMOS-AGN-photzcheck.fits: convert_gal.py COSMOS2015_Laigle+.fits chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits cosmos_acs_iphot_200709.coords filternames.txt
	PHOTOZ_REPLACE_ZSPEC=1 python3 convert_agn.py COSMOS2015_Laigle+.fits chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits $@
input-COSMOS-AGN-specz-blended.fits: convert_gal.py COSMOS2015_Laigle+.fits chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits cosmos_acs_iphot_200709.coords filternames.txt
	UNBLENDED=0 python3 convert_agn.py COSMOS2015_Laigle+.fits chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits $@
input-COSMOS-AGN-specz-shallow.fits: convert_gal.py COSMOS2015_Laigle+.fits chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits cosmos_acs_iphot_200709.coords filternames.txt
	MAKESHALLOW=1 python3 convert_agn.py COSMOS2015_Laigle+.fits chandra_COSMOS_legacy_opt_NIR_counterparts_20160113_4d.fits $@

input-full.fits: COSMOS2015_Laigle+.fits_filtered_gal_zspec_all.fits
	cp $^ $@

cosmos-ids-remove-donley.txt: Donleycrit.py COSMOS2015_Laigle+.fits_filtered_gal_zspec_all.fits
	# identify IR AGN
	python3 $^ $@

input-full1.fits: input-full.fits cosmos-ids-remove-radio.txt
	# remove radio AGN
	stilts tmatch2 in1=input-full.fits values1=id ifmt2=ASCII in2=cosmos-ids-remove-radio.txt values2=col1  matcher=exact join=1not2 find=all out=$@

input-full2.fits: input-full1.fits cosmos-ids-remove-simbad.txt
	# remove AGN known by simbad classification
	stilts tmatch2 in1=input-full1.fits values1=id ifmt2=ASCII in2=cosmos-ids-remove-simbad.txt values2=col1  matcher=exact join=1not2 find=all out=$@

input.fits: input-full2.fits cosmos-ids-remove-donley.txt
	# remove IR AGN
	stilts tmatch2 in1=input-full2.fits values1=id ifmt2=ASCII in2=cosmos-ids-remove-donley.txt values2=col1  matcher=exact join=1not2 find=all out=$@

input-cigale.fits: input.fits
	# convert to CIGALE-compatible input column format
	topcat -stilts tpipe in=$< out=$@ \
		cmd='colmeta -name "subaru.suprime.IB427" IB427_SCam' \
		cmd='colmeta -name "subaru.suprime.IB427_err" IB427_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB464" IB464_SCam' \
		cmd='colmeta -name "subaru.suprime.IB464_err" IB464_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB484" IB484_SCam' \
		cmd='colmeta -name "subaru.suprime.IB484_err" IB484_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB505" IB505_SCam' \
		cmd='colmeta -name "subaru.suprime.IB505_err" IB505_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB527" IB527_SCam' \
		cmd='colmeta -name "subaru.suprime.IB527_err" IB527_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB574" IB574_SCam' \
		cmd='colmeta -name "subaru.suprime.IB574_err" IB574_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB624" IB624_SCam' \
		cmd='colmeta -name "subaru.suprime.IB624_err" IB624_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB679" IB679_SCam' \
		cmd='colmeta -name "subaru.suprime.IB679_err" IB679_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB709" IB709_SCam' \
		cmd='colmeta -name "subaru.suprime.IB709_err" IB709_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB738" IB738_SCam' \
		cmd='colmeta -name "subaru.suprime.IB738_err" IB738_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB767" IB767_SCam' \
		cmd='colmeta -name "subaru.suprime.IB767_err" IB767_SCam_err' \
		cmd='colmeta -name "subaru.suprime.IB827" IB827_SCam' \
		cmd='colmeta -name "subaru.suprime.IB827_err" IB827_SCam_err' \
		cmd='colmeta -name "subaru.suprime.NB711" NB711_SCam' \
		cmd='colmeta -name "subaru.suprime.NB711_err" NB711_SCam_err' \
		cmd='colmeta -name "subaru.suprime.NB816" NB816_SCam' \
		cmd='colmeta -name "subaru.suprime.NB816_err" NB816_SCam_err' \
		cmd='colmeta -name "cfht.wircam.H" CFHT_H' \
		cmd='colmeta -name "cfht.wircam.H_err" CFHT_H_err' \
		cmd='colmeta -name "cfht.wircam.Ks" CFHT_K' \
		cmd='colmeta -name "cfht.wircam.Ks_err" CFHT_K_err' \
		cmd='colmeta -name "subaru.hsc.y" HSC_y' \
		cmd='colmeta -name "subaru.hsc.y_err" HSC_y_err' 

clean:
	rm -rf COSMOS_galaxies_convertedPP.fits cosmos-ids-remove-donley.txt input-full*.fits input.fits

.PHONY: clean all
.SECONDARY:
.SUFFIXES:
