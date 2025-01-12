params.input = null
params.output_dir = 'starsolo_out'
params.star_index = '/home/gjsx/append-ssd/work/star_hs_v46'
params.cb_whitelist = null
// default cb & umi len of star is the same as 10x 3' v2 & 5'
params.cb_len = null
params.umi_len = null
params.solo_feat = "Gene GeneFull_Ex50pAS Velocyto"
params.multimap = "Unique EM"
params.add_cr_conf = "--outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4"
params.add_5p_conf = " "
params.debug = false

process starAlign {
  publishDir "${params.output_dir}/${sample_id}", mode: 'move'
  conda '/home/gjsx/micromamba/envs/star'
  tag "Aligning ${sample_id}"
  cpus  16
  memory 100.GB
  maxForks 3
  
  input:
  tuple val(sample_id), val(readgroup), path(reads1), path(reads2)
  path star_index
  val cb_whitelist
  val cb_len
  val umi_len
  val solo_feat
  
  output:
  path "Aligned.sortedByCoord.out.bam"
  path "Solo.out"
  
  script:
  """
  ulimit -HSn 16000
  
  STAR --genomeDir $star_index --runThreadN ${task.cpus} \
--outSAMattributes NH CB UB \
\
--soloCellFilter EmptyDrops_CR ${params.add_cr_conf} ${params.add_5p_conf} \
--soloMultiMappers ${params.multimap} \
\
--soloType CB_UMI_Simple --readFilesCommand zcat --soloBarcodeReadLength 0 \
--soloCBlen $cb_len --soloUMIstart ${cb_len+1} --soloUMIlen $umi_len \
--soloCBwhitelist $cb_whitelist \
--soloFeatures $solo_feat --outSAMtype BAM SortedByCoordinate \
--outSAMattrRGline 'ID:${readgroup}' 'SM:${readgroup}' 'PL:illumina' \
--readFilesIn ${reads2} ${reads1}

cp -vr Solo.out/GeneFull_Ex50pAS/raw Solo.out/GeneFull_Ex50pAS/raw_full
mv -v Solo.out/GeneFull_Ex50pAS/raw_full/UniqueAndMult-EM.mtx \
Solo.out/GeneFull_Ex50pAS/raw_full/matrix.mtx

STAR --runMode soloCellFiltering \
Solo.out/GeneFull_Ex50pAS/raw_full Solo.out/GeneFull_Ex50pAS/filtered_full/ \
--soloCellFilter EmptyDrops_CR

pigz -vk Solo.out/*/filtered*/*
  """
}

process sleepEcho {
  tag "hearing ${sample_id}"
  debug true
  maxForks 7
  input:
  tuple val(sample_id), val(readgroup), path(reads1), path(reads2)
  path star_index
  val cb_whitelist
  
  output:
  val sample_id
  
  script:
  """
  echo "sample id: " ${sample_id}
  echo "index: " ${star_index}
  echo "whitelist: " ${cb_whitelist}
  echo "read group:" ${readgroup}
  echo "read1: " ${reads1}
  echo "read2: " ${reads2}
  sleep 7s
  """
}

workflow {
  log.info """\
    BILL RNAVAR   P I P E L I N E
    ===================================
    cb_whitelist  : ${params.cb_whitelist}
    cb_len        : ${params.cb_len}
    umi_len       : ${params.umi_len}
    debug         : ${params.debug}
    """
    .stripIndent()

  Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .set{ ch_input }
    
  ch_index = Channel.fromPath(params.star_index)
  
  ch_whitelist = params.cb_whitelist ? file(params.cb_whitelist) : "None"
  //params.add_5p_conf = params.add_5p_conf ? params.add_5p_conf : " "
  
  if (params.debug) {
    sleepEcho(ch_input, ch_index.first(), ch_whitelist)
  } else {
    starAlign(ch_input, ch_index.first(), ch_whitelist, params.cb_len,
    params.umi_len, params.solo_feat)
  }
  
  
  
}