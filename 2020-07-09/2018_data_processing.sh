step=$1

#===========================================
googlebucket="fc-f3db5671-96a8-48d4-bd1c-00a7decfa0a5"

b=/ahg/regevdata2/projects/Cell2CellCommunication/bashrcs/cellranger.bashrc
#source ${b}
CODE=/ahg/regevdata/projects/bn10_cancer_variants/bn10_oana_2018-06-22/src/bn10_analysis
OUT=/ahg/regevdata/projects/bn10_cancer_variants/bn10_oana_2018-06-22/results/2018-06-22
original_sheet=${OUT}/data/sample_sheets/sample_sheet.csv
if [ ! -d "${OUT}" ]; 
then 
    mkdir -p ${OUT}
fi
mypython=/ahg/regevdata/users/oursu/software/anaconda2/bin
myR=/ahg/regevdata/users/oursu/software/anaconda2/bin
bcldirs="/ahg/regev_gp_transfers/180515_SL-HXR_0665_AFCHH2H3CCXY /ahg/regev_gp_transfers/180515_SL-HXR_0666_BFCHGY7FCCXY /ahg/regev_gp_transfers/180515_SL-HXS_0556_BFCHH2FNCCXY /ahg/regev_gp_transfers/180515_SL-HXS_0555_AFCHGYTGCCXY" 
#files
transcriptome_ref=/seq/regev_genome_portal/SOFTWARE/10X/refdata-cellranger-1.2.0/refdata-cellranger-GRCh38-1.2.0

#todo:
#-fastQC

#files that will be made
noncloud_sheet=${OUT}/data/sample_sheets/sample_sheet_noncloud.csv

#===========================================                      
                                                                                                                    
if [[ ${step} == "make_noncloud_samplesheet" ]];
then
    zcat -f ${original_sheet} | grep p53 | sed 's/,/\t/g' | awk '{print $1","$2"_"NR","$3}' > ${noncloud_sheet}
    zcat -f ${original_sheet} | grep KRAS | sed 's/,/\t/g' | awk '{print $1","$2"_"NR","$3}' >> ${noncloud_sheet}
    cat ${noncloud_sheet}
fi

if [[ ${step} == "make_cloud_samplesheet" ]];
then
    echo "Sample,Reference,Flowcell,Lane,Index" > ${OUT}/data/sample_sheets/sample_sheet_cloud.mkfastq.csv
    echo "Sample,Reference,Flowcell,Lane,Index" > ${OUT}/data/sample_sheets/sample_sheet_cloud.count.csv
    for bcldir in $(echo ${bcldirs});
    do
	if [[ $(basename ${bcldir}) == "180515_SL-HXR_0665_AFCHH2H3CCXY" ]];
        then
            dataname="HH2H3CCXY"
        fi
        if [[ $(basename ${bcldir}) == "180515_SL-HXR_0666_BFCHGY7FCCXY" ]];
        then
            dataname="HGY7FCCXY"
        fi
        if [[ $(basename ${bcldir}) == "180515_SL-HXS_0555_AFCHGYTGCCXY" ]];
        then
            dataname="HGYTGCCXY"
        fi
        if [[ $(basename ${bcldir}) == "180515_SL-HXS_0556_BFCHH2FNCCXY" ]];
        then
            dataname="HH2FNCCXY"
        fi
	cat ${noncloud_sheet} | sed 's/,/\t/g' | awk -v bclname=$(basename ${bcldir}) '{print $2",GRCh38,gs://fc-630bd06d-7483-4b4d-beb9-4eebf77d8872/"bclname","$1","$3}' >> ${OUT}/data/sample_sheets/sample_sheet_cloud.mkfastq.csv
	cat ${noncloud_sheet} | sed 's/,/\t/g' | awk -v bclname=$(basename ${bcldir}) -v dname=${dataname} '{print $2",GRCh38,gs://fc-630bd06d-7483-4b4d-beb9-4eebf77d8872/mapped/"bclname"_fastqs/fastq_path/"dname","$1","$3}' >> ${OUT}/data/sample_sheets/sample_sheet_cloud.count.csv
	echo ${OUT}/data/sample_sheets/sample_sheet_cloud.count.csv
    done
    #old firecloud account
    mkdir -p ${OUT}/data/sample_sheets/count_by_sample
    for sample in $(zcat -f ${OUT}/data/sample_sheets/sample_sheet_cloud.count.csv | grep -v Lane | sed 's/,/\t/g' | cut -f1 | sort | uniq);
    do
	echo "Sample,Reference,Flowcell,Lane,Index" > ${OUT}/data/sample_sheets/count_by_sample/sample_sheet_cloud.count.${sample}.csv
	zcat -f ${OUT}/data/sample_sheets/sample_sheet_cloud.count.csv | grep "${sample}," >> ${OUT}/data/sample_sheets/count_by_sample/sample_sheet_cloud.count.${sample}.csv
    done
    #firecloud for this project
    mkdir -p ${OUT}/data/sample_sheets/count_by_sample_${googlebucket}
    for sample in $(zcat -f ${OUT}/data/sample_sheets/sample_sheet_cloud.count.csv | grep -v Lane | sed 's/,/\t/g' | cut -f1 | sort | uniq);
    do
        echo "Sample,Reference,Flowcell,Lane,Index" > ${OUT}/data/sample_sheets/count_by_sample_${googlebucket}/sample_sheet_cloud.count.${sample}.csv
        zcat -f ${OUT}/data/sample_sheets/sample_sheet_cloud.count.csv | grep "${sample}," | sed 's/,/\t/g' | sed 's/fastq_path\//\t/g' | awk -v gb=${googlebucket} '{print $1","$2",gs://"gb"/"$4","$5","$6}' >> ${OUT}/data/sample_sheets/count_by_sample_${googlebucket}/sample_sheet_cloud.count.${sample}.csv
    done
    ls ${OUT}/data/sample_sheets/count_by_sample_${googlebucket}/sample_sheet_cloud.count.${sample}.csv
fi


#actually did this in the cloud
if [[ ${step} == "mkfastq" ]];
then
    for bcldir in $(echo ${bcldirs});
    do
	echo ${bcldir}
	outdir=${OUT}/data/cellranger/mkfastq_local/$(basename ${bcldir})
	mkdir -p ${outdir}
	s=${outdir}/bcl2fastq___$(basename ${bcldir}).sh
	sample_sheet_csv=${noncloud_sheet}
	cmd="cellranger mkfastq --run=${bcldir} --output-dir=${outdir} --csv=${sample_sheet_csv}"
	echo "source ${b}" > ${s}
	echo "cd ${outdir}" >> ${s}
	echo ${cmd} >> ${s}
	chmod 755 ${s}
	echo "qsub -l h_vmem=20G -l h_rt=30:00:00 -l s_rt=30:00:00 -o ${s}.o -e ${s}.e ${s}"
    done
fi


if [[ ${step} == "count_cloud_upload" ]];
then
    for id in A549_KRAS_1 A549_KRAS_2;
    do
	for bcldir in $(echo ${bcldirs});                                                                                                   
	do
            if [[ $(basename ${bcldir}) == "180515_SL-HXR_0665_AFCHH2H3CCXY" ]];
            then
		dataname="HH2H3CCXY"
            fi
            if [[ $(basename ${bcldir}) == "180515_SL-HXR_0666_BFCHGY7FCCXY" ]];
            then
		dataname="HGY7FCCXY"
            fi
            if [[ $(basename ${bcldir}) == "180515_SL-HXS_0555_AFCHGYTGCCXY" ]];
            then
		dataname="HGYTGCCXY"
            fi
            if [[ $(basename ${bcldir}) == "180515_SL-HXS_0556_BFCHH2FNCCXY" ]];
            then
		dataname="HH2FNCCXY"
            fi
            fastq_dir=${OUT}/data/cellranger/mkfastq/$(basename ${bcldir})_fastqs/fastq_path/${dataname}/
            #ls -lh ${fastq_dir}/${id}
	    gsutil -m cp -r ${fastq_dir}/${id} gs://${googlebucket}/${dataname}
        done
	gsutil -m cp ${OUT}/data/sample_sheets/count_by_sample_${googlebucket}/sample_sheet_cloud.count.${id}.csv gs://${googlebucket}
    done
    
fi



if [[ ${step} == "count" ]];
then
    for id in A549_KRAS_1;
    do
	countdir=${OUT}/data/cellranger/count_local/${id}
        mkdir -p ${countdir}
        s=${countdir}/cellranger_count___${id}.sh
	#mkdir -p ${OUT}/data/cellranger/count_local/${id}_fastqs
	fastq_dirs=""
	for bcldir in $(echo ${bcldirs});
        do
            if [[ $(basename ${bcldir}) == "180515_SL-HXR_0665_AFCHH2H3CCXY" ]];
            then
                dataname="HH2H3CCXY"
            fi
            if [[ $(basename ${bcldir}) == "180515_SL-HXR_0666_BFCHGY7FCCXY" ]];
            then
                dataname="HGY7FCCXY"
            fi
            if [[ $(basename ${bcldir}) == "180515_SL-HXS_0555_AFCHGYTGCCXY" ]];
            then
                dataname="HGYTGCCXY"
            fi
            if [[ $(basename ${bcldir}) == "180515_SL-HXS_0556_BFCHH2FNCCXY" ]];
            then
                dataname="HH2FNCCXY"
            fi
	    #mkdir -p ${OUT}/data/cellranger/count_local/${id}_fastqs/${dataname}
	    #cp -r ${OUT}/data/cellranger/mkfastq/*${dataname}_fastqs/fastq_path/${dataname}/${id} ${OUT}/data/cellranger/count_local/${id}_fastqs/${dataname}/
	    fastq_dirs=${OUT}/data/cellranger/mkfastq/$(basename ${bcldir})_fastqs/fastq_path/${dataname}/${id}","${fastq_dirs}
	done
	fastq_dirs_2=$(echo ${fastq_dirs} | sed 's/,$//g')
	echo ${fastq_dirs_2}
	cmd="cellranger count --id=${id} --sample=${id} --transcriptome=${transcriptome_ref} --force-cells=6000 --fastqs=${fastq_dirs_2}"
	echo "source ${b}" > ${s}
	echo "cd ${countdir}" >> ${s}
	echo "${cmd}" >> ${s}
	echo ${s}
	chmod 755 ${s}
	#10 hours is fine! actually, some needed 20h.
	qsub -pe smp 12 -binding linear:12 -l h_vmem=12G -R y -l h_rt=20:00:00 -l s_rt=20:00:00 -e ${s}.e -o ${s}.o ${s}
    done
fi

if [[ ${step} == "aggregate" ]];
then
    for experiment in KRAS;
    do
	aggr_name=${experiment}_channels_1to32
        #make csv
	mkdir -p ${OUT}/data/cellranger/aggr
	csv_aggregation=${OUT}/data/cellranger/aggr/aggr.${aggr_name}.csv
	rm ${csv_aggregation}
	echo "library_id,molecule_h5" > ${csv_aggregation}
	for channel in {1..32};
	do
	    echo "A549_${experiment}_${channel},${OUT}/data/cellranger/count_local/A549_${experiment}_${channel}/A549_${experiment}_${channel}/outs/molecule_info.h5" >> ${csv_aggregation}
	done
	mkdir -p ${OUT}/data/cellranger/aggr/${aggr_name}
	outpath=${OUT}/data/cellranger/aggr/${aggr_name}/${aggr_name}
	s=${outpath}.sh
	cmd="cellranger aggr --id ${aggr_name} --csv ${csv_aggregation} --normalize=mapped --nosecondary"
	echo "source ${b}" > ${s}
	echo "cd ${OUT}/data/cellranger/aggr/${aggr_name}" >> ${s}
	echo "${cmd}" >> ${s}
	echo ${s}
	chmod 755 ${s}
	echo "qsub -pe smp 12 -binding linear:12 -l h_vmem=12G -R y -l h_rt=30:00:00 -l s_rt=30:00:00 -e ${s}.e -o ${s}.o ${s}"
    done
fi

#=============================================
# counting barcodes
#=============================================

if [[ ${step} == "countBC_bowtie_build" ]];
then
    mkdir -p ${OUT}/data/variant_barcodes
    barcode_fasta=${OUT}/data/variant_barcodes/variant_barcodes.fa
    /ahg/regevdata/users/oursu/software/anaconda2/bin/python ${CODE}/pipeline/process_variant_barcodes.py --variant_barcodes ${OUT}/data/metadata/TP53.txt,${OUT}/data/metadata/KRAS.txt --out ${barcode_fasta}
    #add the puto gene to the barcode fasta
    puro_seq=$(echo "atgaccgagtacaagcccacggtgcgcctcgccacccgcgatgatgtccctagggccgtacgcaccctcgccgccgcgttcgccgactaccccgccacgcgccacaccgtcgatccggaccgccacatcgagcgggtcaccgagctgcaagaactcttcctcacgcgcgtcgggctcgacatcggcaaggtgtgggtcgctgatgacggcgccgccgtggcggtctggaccacgccggagagcgtcgaagcgggggcggtgttcgccgagatcggcccgcgcatggccgagttgagcggttcccggctggccgcgcagcaacagatggaaggcctcctggcgccgcaccggcccaaggagcccgcgtggttcctggccaccgtcggcgtttcgcccgaccaccagggcaagggtctgggcagcgccgtcgtgctccccggagtggaagctgccgagcgcgccggggtgcctgccttcctggagacctccgcgccccgcaacctccccttctacgagcggctcggcttcaccgtcaccgccgacgtcgaggtgcccgaaggaccgcgcacctggtgcatgacccgcaagcccggtgccagttaaggatccgttaacgtcgagggacctaataacttcgtatagcatacattatacgaagttatacatgtttaagggttccggttccactaggtacaattcgatatcaagcttatcgataatcaacctctggattacaaaatttgtgaaagattgactggtattcttaactatgttgctccttttacgctatgtggatacgctgctttaatgcctttgtatcatgctattgcttcccgtatggctttcattttctcctccttgtataaatcctggttgctgtctctttatgaggagttgtggcccgttgtcaggcaacgtggcgtggtgtgcactgtgtttgctgacgcaacccccactggttggggcattgccaccacctgtcagctcctttccgggactttcgctttccccctccctattgccacggcggaactcatcgccgcctgccttgcccgctgctggacaggggctcggctgttgggcactgacaattccgtggtgttgtcggggaaatcatcgtcctttccttggctgctcgcctgtgttgccacctggattctgcgcgggacgtccttctgctacgtcccttcggccctcaatccagcggaccttccttcccgcggcctgctgccggctctgcggcctcttccgcgtcttcgccttcgccctcagacgagtcggatctccctttgggccgcctccccgcatcgataccgtcgacctcgatcgagacctagaaaaacatggagcaatcacaagtagcaatacagcagctaccaatgctgattgtgcctggctagaagcacaagaggaggaggaggtgggttttccagtcacacctcaggtacctttaagaccaatgacttacaaggcagctgtagatcttagccactttttaaaagaaaaggggggactggaagggctaattcactcccaacgaagacaagatatccttgatctgtggatctaccacacacaaggctacttccctgattggcagaactacacaccagggccagggatcagatatccactgacctttggatggtgctacaagctagtaccagttgagcaagagaaggtagaagaagccaatgaaggagagaacacccgcttgttacaccctgtgagcctgcatgggatggatgacccggagagagaagtattagagtggaggtttgacagccgcctagcatttcatcacatggcccgagagctgcatccggactgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagca" | tr '[:lower:]' '[:upper:]') #this starts at CDS_4 and ends at the end of LTR3 in the plasmid
    #p53_wt_construct=$(echo "gcttaagcggtcgacggatcgggagatctcccgatcccctatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagtatctgctccctgcttgtgtgttggaggtcgctgagtagtgcgcgagcaaaatttaagctacaacaaggcaaggcttgaccgacaattgcatgaagaatctgcttagggttaggcgttttgcgctgcttcgcgatgtacgggccagatatacgcgttgacattgattattgactagttattaatagtaatcaattacggggtcattagttcatagcccatatatggagttccgcgttacataacttacggtaaatggcccgcctggctgaccgcccaacgacccccgcccattgacgtcaataatgacgtatgttcccatagtaacgccaatagggactttccattgacgtcaatgggtggagtatttacggtaaactgcccacttggcagtacatcaagtgtatcatatgccaagtacgccccctattgacgtcaatgacggtaaatggcccgcctggcattatgcccagtacatgaccttatgggactttcctacttggcagtacatctacgtattagtcatcgctattaccatggtgatgcggttttggcagtacatcaatgggcgtggatagcggtttgactcacggggatttccaagtctccaccccattgacgtcaatgggagtttgttttggcaccaaaatcaacgggactttccaaaatgtcgtaacaactccgccccattgacgcaaatgggcggtaggcgtgtacggtgggaggtctatataagcagcgcgttttgcctgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagtggcgcccgaacagggacttgaaagcgaaagggaaaccagaggagctctctcgacgcaggactcggcttgctgaagcgcgcacggcaagaggcgaggggcggcgactggtgagtacgccaaaaattttgactagcggaggctagaaggagagagatgggtgcgagagcgtcagtattaagcgggggagaattagatcgcgatgggaaaaaattcggttaaggccagggggaaagaaaaaatataaattaaaacatatagtatgggcaagcagggagctagaacgattcgcagttaatcctggcctgttagaaacatcagaaggctgtagacaaatactgggacagctacaaccatcccttcagacaggatcagaagaacttagatcattatataatacagtagcaaccctctattgtgtgcatcaaaggatagagataaaagacaccaaggaagctttagacaagatagaggaagagcaaaacaaaagtaagaccaccgcacagcaagcggccggccgcgctgatcttcagacctggaggaggagatatgagggacaattggagaagtgaattatataaatataaagtagtaaaaattgaaccattaggagtagcacccaccaaggcaaagagaagagtggtgcagagagaaaaaagagcagtgggaataggagctttgttccttgggttcttgggagcagcaggaagcactatgggcgcagcgtcaatgacgctgacggtacaggccagacaattattgtctggtatagtgcagcagcagaacaatttgctgagggctattgaggcgcaacagcatctgttgcaactcacagtctggggcatcaagcagctccaggcaagaatcctggctgtggaaagatacctaaaggatcaacagctcctggggatttggggttgctctggaaaactcatttgcaccactgctgtgccttggaatgctagttggagtaataaatctctggaacagatttggaatcacacgacctggatggagtgggacagagaaattaacaattacacaagcttaatacactccttaattgaagaatcgcaaaaccagcaagaaaagaatgaacaagaattattggaattagataaatgggcaagtttgtggaattggtttaacataacaaattggctgtggtatataaaattattcataatgatagtaggaggcttggtaggtttaagaatagtttttgctgtactttctatagtgaatagagttaggcagggatattcaccattatcgtttcagacccacctcccaaccccgaggggacccgacaggcccgaaggaatagaagaagaaggtggagagagagacagagacagatccattcgattagtgaacggatcggcactgcgtgcgccaattctgcagacaaatggcagtattcatccacaattttaaaagaaaaggggggattggggggtacagtgcaggggaaagaatagtagacataatagcaacagacatacaaactaaagaattacaaaaacaaattacaaaaattcaaaattttcgggtttattacagggacagcagagatccagtttggttagtaccgggcccgctctagaccatagagcccaccgcatccccagcatgcctgctattgtcttcccaatcctcccccttgctgtcctgccccaccccaccccccagaatagaatgacacctactcagacaatgcgatgcaatttcctcattttattaggaaaggacagtgggagtggcaccttccagggtcaaggaaggcacgggggaggggcaaacaacagatggctggcaactagaaggcacagtcgaggctgatcagcgggtttaaacgggccctctaggttaattaacaccacGCTAGCTTCTTCTTCTTCCGCAGCCCGGGCATTTTGTCTCAGTCTGAGTCAGGCCCTTCTGTCTTGAACATGAGTTTTTTATGGCGGGAGGTAGACTGACCCTTTTTGGACTTCAGGTGGCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGAGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGAAGGGTGAAATATTCTCCATCCAGTGGTTTCTTCTTTGGCTGGGGAGAGGAGCTGGTGTTGTTGGGCAGTGCTCGCTTAGTGCTCCCTGGGGGCAGCTCGTGGTGAGGCTCCCCTTTCTTGCGGAGATTCTCTTCCTCTGTGCGCCGGTCTCTCCCAGGACAGGCACAAACACGCACCTCAAAGCTGTTCCGTCCCAGTAGATTACCACTGGAGTCTTCCAGTGTGATGATGGTGAGGATGGGCCTCCGGTTCATGCCGCCCATGCAGGAACTGTTACACATGTAGTTGTAGTGGATGGTGGTACAGTCAGAGCCAACCTCAGGCGGCTCATAGGGCACCACCACACTATGTCGAAAAGTGTTTCTGTCATCCAAATACTCCACACGCAAATTTCCTTCCACTCGGATAAGATGCTGAGGAGGGGCCAGACCATCGCTATCTGAGCAGCGCTCATGGTGGGGGCAGCGCCTCACAACCTCCGTCATGTGCTGTGACTGCTTGTAGATGGCCATGGCGCGGACGCGGGTGCCGGGCGGGGGTGTGGAATCAACCCACAGCTGCACAGGGCAGGTCTTGGCCAGTTGGCAAAACATCTTGTTGAGGGCAGGGGAGTACGTGCAAGTCACAGACTTGGCTGTCCCAGAATGCAAGAAGCCCAGACGGAAACCGTAGCTGCCCTGGTAGGTTTTCTGGGAAGGGACAGAAGATGACAGGGGCCAGGAGGGGGCTGGTGCAGGGGCCGCCGGTGTAGGAGCTGCTGGTGCAGGGGCCACGGGGGGAGCAGCCTCTGGCATTCTGGGAGCTTCATCTGGACCTGGGTCTTCAGTGAACCATTGTTCAATATCGTCCGGGGACAGCATCAAATCATCCATTGCTTGGGACGGCAAGGGGGACAGAACGTTGTTTTCAGGAAGTAGTTTCCATAGGTCTGAAAATGTTTCCTGACTCAGAGGGGGCTCGACGCTAGGATCTGACTGCGGCTCCTCCATGGTGGCgctagctcacgacacctgaaatggaagaaaaaaactttgaaccactgtctgaggcttgagaatgaaccaagatccaaactcaaaaagggcaaattccaaggagaattacatcaagtgccaagctggcctaacttcagtctccacccactcagtgtggggaaactccatcgcataaaacccctccccccaacctaaagacgacgtactccaaaagcgcgagaactaatcgaggtgcctggacggcgcccggtactccgtggagtcacatgaagcgacggctgaggacggaaaggcccttttcctttgtgtgggtgactcacccgcccgctctcccgagcgccgcgtcctccattttgagctccctgcagcagggccgggaagcggccatctttccgctcacgcaactggtgccgaccgggccagccttgccgcccagggcggggcgatacacggcggcgcgaggccaggcaccagagcaggccggccagcttgagactacccccgtccgattctcggtggccgcgctcgcaggccccgcctcgccgaacatgtgcgctgggacgcacgggccccgtcgccgcccgcggccccaaaaaccgaaataccagtgtgcagatcttggcccgcatttacaagactatcttgccagaaaaaaagcgtcgcagcaggtcatcaaaaattttaaatggctagagacttatcgaaagcagcgagacaggcgcgaaggtgccaccagattcgcacgcggcggccccagcgcccaggccaggcctcaactcaagcacgaggcgaaggggctccttaagcgcaaggcctcgaactctcccacccacttccaacccgaagctcgggatcaagaatcacgtactgcagccaggtggaagtaattcaaggcacgcaagggccataacccgtaaagaggccaggcccgcgggaaccacacacggcacttacctgtgttctggcggcaaacccgttgcgaaaaagaacgttcacggcgactactgcacttatatacggttctcccccaccctcgggaaaaaggcggagccagtacacgacatcactttcccagtttaccccgcgccaccttctctaggcaccggttcaattgccgacccctccccccaacttctcggggactgtgggcgatgtgcgctctgcccactgacgggcacatgcatggcggtaatacggttatccacgcggccgcgttaatggaacaggaactaaatttacccccgggtaggggaggcgcttttcccaaggcagtctggagcatgcgctttagcagccccgctgggcacttggcgctacacaagtggcctctggcctcgcacacattccacatccaccggtaggcgccaaccggctccgttctttggtggccccttcgcgccaccttctactcctcccctagtcaggaagttcccccccgccccgcagctcgcgtcgtgcaggacgtgacaaatggaagtagcacgtctcactagtctcgtgcagatggacagcaccgctgagcaatggaagcgggtaggcctttggggcagcggccaatagcagctttgctccttcgctttctgggctcagaggctgggaaggggtgggtccgggggcgggctcaggggcgggctcaggggcggggcgggcgcccgaaggtcctccggaggcccggcattctgcacgcttcaaaagcgcacgtctgccgcgctgttctcctcttcctcatctccgggcctttcgctcgaggccaccatgaccgagtacaagcccacggtgcgcctcgccacccgcgatgatgtccctagggccgtacgcaccctcgccgccgcgttcgccgactaccccgccacgcgccacaccgtcgatccggaccgccacatcgagcgggtcaccgagctgcaagaactcttcctcacgcgcgtcgggctcgacatcggcaaggtgtgggtcgctgatgacggcgccgccgtggcggtctggaccacgccggagagcgtcgaagcgggggcggtgttcgccgagatcggcccgcgcatggccgagttgagcggttcccggctggccgcgcagcaacagatggaaggcctcctggcgccgcaccggcccaaggagcccgcgtggttcctggccaccgtcggcgtttcgcccgaccaccagggcaagggtctgggcagcgccgtcgtgctccccggagtggaagctgccgagcgcgccggggtgcctgccttcctggagacctccgcgccccgcaacctccccttctacgagcggctcggcttcaccgtcaccgccgacgtcgaggtgcccgaaggaccgcgcacctggtgcatgacccgcaagcccggtgccagttaaggatccgttaacgtcgagggacctaataacttcgtatagcatacattatacgaagttatacatgtttaagggttccggttccactaggtacaattcgatatcaagcttatcgataatcaacctctggattacaaaatttgtgaaagattgactggtattcttaactatgttgctccttttacgctatgtggatacgctgctttaatgcctttgtatcatgctattgcttcccgtatggctttcattttctcctccttgtataaatcctggttgctgtctctttatgaggagttgtggcccgttgtcaggcaacgtggcgtggtgtgcactgtgtttgctgacgcaacccccactggttggggcattgccaccacctgtcagctcctttccgggactttcgctttccccctccctattgccacggcggaactcatcgccgcctgccttgcccgctgctggacaggggctcggctgttgggcactgacaattccgtggtgttgtcggggaaatcatcgtcctttccttggctgctcgcctgtgttgccacctggattctgcgcgggacgtccttctgctacgtcccttcggccctcaatccagcggaccttccttcccgcggcctgctgccggctctgcggcctcttccgcgtcttcgccttcgccctcagacgagtcggatctccctttgggccgcctccccgcatcgataccgtcgacctcgatcgagacctagaaaaacatggagcaatcacaagtagcaatacagcagctaccaatgctgattgtgcctggctagaagcacaagaggaggaggaggtgggttttccagtcacacctcaggtacctttaagaccaatgacttacaaggcagctgtagatcttagccactttttaaaagaaaaggggggactggaagggctaattcactcccaacgaagacaagatatccttgatctgtggatctaccacacacaaggctacttccctgattggcagaactacacaccagggccagggatcagatatccactgacctttggatggtgctacaagctagtaccagttgagcaagagaaggtagaagaagccaatgaaggagagaacacccgcttgttacaccctgtgagcctgcatgggatggatgacccggagagagaagtattagagtggaggtttgacagccgcctagcatttcatcacatggcccgagagctgcatccggactgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagcatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagggcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttgccattgctacaggcatcgtggtgtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataataccgcgccacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccacctgac" | tr '[:lower:]' '[:upper:]')
    echo ">Puro_res" >> ${barcode_fasta}
    echo ${puro_seq} >> ${barcode_fasta}
    #echo ">whole_construct" >> ${barcode_fasta}
    #echo ${p53_wt_construct} >> ${barcode_fasta}
    echo ${barcode_fasta}
    s=${OUT}/data/variant_barcodes/variant_barcodes.index.sh
    echo "${mypython}/bowtie-build -f ${barcode_fasta} ${barcode_fasta}.index" > ${s}
    echo ${s}
    chmod 755 ${s}
fi

if [[ ${step} == "countBC_align" ]];
then
    index=${OUT}/data/variant_barcodes/variant_barcodes.fa.index
    for bcldir in $(echo ${bcldirs});
    do
	dataname=$(basename ${bcldir} | sed 's/AFC/\t/g' | sed 's/BFC/\t/g' | cut -f2)
	echo "============ ${dataname}"
	for sample in $(zcat -f ${OUT}/data/sample_sheets/sample_sheet_cloud.csv | grep -v Lane | sed 's/,/\t/g' | cut -f1 | sort | uniq);
	do
	    for mm in 1 2;
	    do
		outdir=${OUT}/data/count_barcodes/mismatches${mm}/$(basename ${bcldir})
		mkdir -p ${outdir}
		general=$(ls ${OUT}/data/cellranger/mkfastq/$(basename ${bcldir})_fastqs/fastq_path/${dataname}/${sample}/${sample}_*_I1_*.fastq.gz | sed 's/I1/SUBSTITUTE/g') #.head2.5M.gz
		echo ${general}
		r1=$(echo ${general} | sed 's/SUBSTITUTE/R1/g')
		r2=$(echo ${general} | sed 's/SUBSTITUTE/R2/g')
		i1=$(echo ${general} | sed 's/SUBSTITUTE/I1/g')
		outpath=${outdir}/${sample}.$(basename ${bcldir}).variantBarcodeCounts.gz
		s=${outpath}.sh
		cmd="${mypython}/bowtie --verbose -n ${mm} -v ${mm} -l 98 -a ${index} ${r2} | gzip > ${outpath}"  
		echo "#!/usr/bin/env bash" > ${s}
		echo "${cmd}" >> ${s}
		cmd2="${mypython}/python ${CODE}/pipeline/add_cellbc_and_UMI_to_mapped_barcodes.py --mapped ${outpath} --R1 ${r1} --out ${outpath}.with_cellbc_and_UMI.gz"
		echo "${cmd2}" >> ${s}
		#echo "${myR}/Rscript ${CODE}/pipeline/plot_variant_barcode_QC.R ${outpath}.with_cellbc_and_UMI.gz" >> ${s}
		chmod 755 ${s}
		echo ${s}
		qsub -o ${s}.o -e ${s}.e ${s}
	    done
	done
    done
fi

if [[ ${step} == "countBC_merge" ]];
then
    for mm in 0;
    do
	for sampletype in p53 KRAS;
	do
	    combined_count=${OUT}/data/count_barcodes/mismatches${mm}/summary/mismatches${mm}.${sampletype}.vbc.counts
	    combined_cells=${OUT}/data/count_barcodes/mismatches${mm}/summary/mismatches${mm}.${sampletype}.cells
	    mkdir -p ${OUT}/data/count_barcodes/mismatches${mm}/summary
	    rm ${combined_count} ${combined_count}.gz ${combined_cells} ${combined_cells}.gz
	    for bcldir in $(echo ${bcldirs}); 
	    do
		for sample in $(zcat -f ${OUT}/data/sample_sheets/sample_sheet_cloud.csv | grep -v Lane | grep ${sampletype} | sed 's/,/\t/g' | cut -f1 | sort | uniq);
		do
		    if [[ ${sample} == *${sampletype}* ]];
		    then
			echo ${sample}
			vbc=${OUT}/data/count_barcodes/mismatches${mm}/$(basename ${bcldir})/${sample}.$(basename ${bcldir}).variantBarcodeCounts.gz.with_cellbc_and_UMI.gz
			zcat -f ${vbc} | awk -v thissample=${sample} -v thisflowcell=$(basename ${bcldir}) '{print thissample"\t"thisflowcell"\t"$0}' >> ${combined_count}
			cells=${vbc}.cells.gz
			zcat -f ${cells} | awk -v thissample=${sample} -v thisflowcell=$(basename ${bcldir}) '{print thissample"\t"thisflowcell"\t"$0}' >> ${combined_cells}
		    fi
		done
	    done
	    #count all cells ::: awk '{print $1"with"$3"\t"$4}' | awk '{a[$1] += $2} END{for (i in a) print i, a[i]}' | gzip
	    zcat -f ${combined_count} | gzip > ${combined_count}.gz
	    zcat -f ${combined_cells} | awk '{print $1"with"$3"\t"$4}' | awk '{a[$1] += $2} END{for (i in a) print i"\t"a[i]}' | gzip > ${combined_cells}.gz
	    rm ${combined_count} ${combined_cells}
	    echo ${combined_count}.gz
	    ${myR}/Rscript ${CODE}/pipeline/plot_variant_barcode_QC.R ${combined_count}.gz
	done
    done
fi
