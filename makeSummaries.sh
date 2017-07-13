#!/bin/bash                                                                                                                                                                                                                              
set -x

wc -l starFusion/*/star-fusion.fusion_candidates.final.abridged.FFPM > counts_per_case.star-fusion.fusion_candidates.txt
wc -l mapsplice/*/fusions_candidates.txt | sed -e s/"'"/p/g -e s/'"'/p/g  > counts_per_case.mapsplice.fusion_candidates.txt
wc -l ericScript/*/*.results.filtered.tsv | sed -e s/"'"/p/g -e s/'"'/p/g  > counts_per_case.eric.results.filtered.txt
wc -l defuse/*/*.defuse*txt | sed -e s/"'"/p/g -e s/'"'/p/g  > counts_per_case.defuse_results.txt

wc -l starFusion/*/star-fusion.oncofuse.output.txt.FFPM | sed -e s/"'"/p/g -e s/'"'/p/g  > counts_per_case.star-fusion.oncofuse.output.txt
wc -l mapsplice/*/mapsplice.oncofuse.output.txt | sed -e s/"'"/p/g -e s/'"'/p/g  > counts_per_case.mapsplice.oncofuse.output.txt
wc -l ericScript/*/eric.oncofuse.output.txt | sed -e s/"'"/p/g -e s/'"'/p/g  > counts_per_case.ericscript.oncofuse.output.txt
wc -l integrate/*/*txt | sed -e s/"'"/p/g -e s/'"'/p/g  > counts_per_case.integrate.oncofuse.output.txt


grep ".*" starFusion/*/star-fusion.fusion_candidates.final.abridged.FFPM | awk '/FusionName/ && FNR > 1 {next} {OFS="\t"; n = split($1, info, "/"); split(info[n], obs, ":"); print info[1], info[2], obs[1], obs[2], $0}' | cut -f 1-4,\
   6- | sed -e s/"'"/p/g -e s/'"'/p/g  -e 's/#//' > merged.star-fusion.fusion_candidates.txt

grep ".*" starFusion/*/star-fusion.oncofuse.output.txt | sed 's-:.//-:-' | awk '/SAMPLE_ID/ && FNR > 1 {next} {OFS="\t"; n = split($1, info, "/"); split(info[3], obs, ":"); print info[1], info[2], obs[1], obs[2], $0}' | cut -f 1-4,6\
   - | sed -e s/"'"/p/g -e s/'"'/p/g  > merged.star-fusion.oncofuse.output.txt

grep ".*" mapsplice/*/fusions_not_well_annotated.txt mapsplice/*/fusions_well_annotated.txt | awk '{OFS="\t"; n = split($1, info, "/"); split(info[n], obs, ":"); print info[1], info[2], obs[1], obs[2], $0}' | cut -f 1-4,6- | tr "~" \
   "\t" | sed -e s/"'"/p/g -e s/'"'/p/g  > merged.mapsplice.fusions_candidates.txt

grep ".*" mapsplice/*/mapsplice.oncofuse.output.txt | sed s-:./-:- | awk '/SAMPLE_ID/ && FNR > 1 {next} {OFS="\t"; n = split($1, info, "/"); split(info[3], obs, ":"); print info[1], info[2], obs[1], obs[2], $0}' | cut -f 1-4,6- | se\
   d -e s/"'"/p/g -e s/'"'/p/g  > merged.mapsplice.oncofuse.output.txt

grep ".*" ericScript/*/*.results.filtered.tsv | awk '/GeneName1/ && FNR > 1 {next} {OFS="\t"; n = split($1, info, "/"); split(info[n], obs, ":"); print info[1], info[2], obs[1], obs[2], $0}' | cut -f 1-4,6- | sed -e s/"'"/p/g -e s/'\
   "'/p/g  > merged.ericscript.results.filtered.txt

grep ".*" ericScript/*/eric.oncofuse.output.txt | sed s-:.//-:- | awk '/SAMPLE_ID/ && FNR > 1 {next} {OFS="\t"; n = split($1, info, "/"); split(info[3], obs, ":"); print info[1], info[2], obs[1], obs[2], $0}' | cut -f 1-4,6- | sed -\
   e s/"'"/p/g -e s/'"'/p/g  > merged.ericscript.oncofuse.output.txt

grep ".*" integrate/*/*oncofuse.output.txt | sed -e "s;integrate_rnaseq/oncofuse/;;" -e "s;.oncofuse.input.txt;;" | awk '/SAMPLE_ID/ && FNR > 1 {next} {OFS="\t"; n = split($1, info, "/"); split(info[3], obs, ":"); print info[1], inf\
   o[2], obs[1], obs[2], $0}' | cut -f 1-4,6- | sed -e s/"'"/p/g -e s/'"'/p/g  > merged.integrate.oncofuse.output.txt


grep .* defuse/*/*defuse.txt | awk ' /cluster_id/ && NFR > 1 {next} {OFS="\t"; n = split($1, info, "/"); split(info[3], obs, ":"); print info[1], info[2], obs[1], obs[2], $0}' | cut -f 1-4,6- | sed -e s/"'"/p/g -e s/'"'/p/g > merged\
   .defuse.results.txt

