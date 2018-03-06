#!/usr/bin/env bash

DEFAULT_FOLDER='/ChIP-seq/'

project=$1
folder=${2:-$DEFAULT_FOLDER}

ASSETS=('awscli_asset' 'bedtools_asset' 'bioconductor_asset' 'common_asset' 'filter_qc_asset' 'idr2_asset' 'macs2_asset' 'samtools_asset' 'spp_asset' 'ucsc_asset')
dx mkdir -p "$project:$folder/assets/"

for asset in ${ASSETS[@]}; do
	dest="$project:$folder/assets/$asset"
	echo $dest >> $asset.log
	dx build_asset --no-watch --destination "$project:$folder/assets/$asset" "$asset/" 1>> $asset.log 2>&1 &
done
