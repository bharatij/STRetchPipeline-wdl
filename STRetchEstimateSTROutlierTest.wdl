#Workflow : STRetchPipeline-Step3: STRetch Estimate allele size and One vs all samples outlier test on large number of WGS samples
workflow STRetchOutlier {
String pipeline_version = "1.0"
meta {
	author: "Bharati Jadhav"
    	email: "bharati.jadhav@mssm.edu"
    	description: "Run STRetch estimate allele size and one vs all outlier test"	
  }
	
  File sample_map
  String output_prefix  #"example: Test.STRetch.Outlier" 
  Array[String] chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
  String gotc_docker = "bharatij/stretch_pipeline:stretch" #updated estimatesSTR.py to add chromosome option and work with thausands of samples
    
  Int preemptible_tries = 3
    
  #get list of medain cov, locus, str counts
  call PartitionSampleNameMap {
    input:
        sample_map = sample_map,
        preemptible_tries = preemptible_tries
  }
  
  Array[File] inputSamples = read_lines(PartitionSampleNameMap.samplenames)
  Array[File] fMedCov = read_lines(PartitionSampleNameMap.medcov)
  Array[File] fLocusCt = read_lines(PartitionSampleNameMap.locusct)
  Array[File] fStrCt =   read_lines(PartitionSampleNameMap.strct)
    
  #estimate STR and Outlier on each chromosome in parallel
  scatter (chrom in chroms) {   
    call RunEstimateSTR { 
	input: 
	  fMedCov =fMedCov,
	  fLocusCt =fLocusCt,
	  fStrCt =fStrCt,
          chrom=chrom,
	  output_prefix = output_prefix,
	  docker_image = gotc_docker,
	  preemptible_tries = preemptible_tries
    }
    
  }
    
  Array[File] partStrs = RunEstimateSTR.strtsvs
  call CombineChunkStr{
    input: 
    	partStrs = partStrs,
        output_prefix = output_prefix,
        docker_image = gotc_docker,
	preemptible_tries = preemptible_tries
  }
     
  File allStrs = CombineChunkStr.allStrs
  
  scatter (sampleid in inputSamples) {
    call RunMultiTestCorrection {
        input: 
          allStrs = allStrs,
          sampleid = sampleid,
          docker_image = gotc_docker,
	  preemptible_tries = preemptible_tries
    }
  }

  Array[File] sampleStrs = RunMultiTestCorrection.sampleStr
  #combine all samples 
  call CombineSamples{
    input: 
      	 sampleStrs = sampleStrs,
         output_prefix = output_prefix,
         docker_image = gotc_docker,
	 preemptible_tries = preemptible_tries
  }
  
  #Outputs
  output {
      	  File mastertable = CombineSamples.mastertable
  }
}
    
#Task Definitions
task PartitionSampleNameMap {
  File sample_map
  Int preemptible_tries
  
  command <<<
    	grep -v ^SampleId ${sample_map} | cut -f1> sample_list
    	grep -v ^SampleId ${sample_map} | cut -f2 > medcov_paths
    	grep -v ^SampleId ${sample_map} | cut -f3> strct_paths
    	grep -v ^SampleId ${sample_map} | cut -f4> locusct_paths
  >>>

  output {
    File samplenames = "sample_list"
    File medcov = "medcov_paths"
    File strct = "strct_paths"
    File locusct = "locusct_paths"
  }

  runtime {
    memory: "1 GiB"
    disks: "local-disk 5 HDD"
    cpu: "1"
    zones: "us-central1-c us-central1-b"
    preemptible: preemptible_tries   
    #docker: "us.gcr.io/broad-gotc-prod/python:2.7"
    docker:"bharatij/stretch_pipeline:stretch"
  }
}
#estimate STRs and outlier
task RunEstimateSTR {
  Array[File] fMedCov
  Array[File] fLocusCt
  Array[File] fStrCt
  String chrom
  String output_prefix
  
  String docker_image
  Int preemptible_tries
  Int addtional_disk_size = 4
  String machine_mem_size = 32 # can be < 5GB for few hundred samples	

  Float outsize = 200 # can be < 2GB for few hundred samples
    
  Int disk_size = ceil(outsize) + addtional_disk_size
  #Calls estimate STR
  command {
   	    PATH=$PATH:/STRetch-master/tools/bin; \ 
            /STRetch-master/tools/bin/python /STRetch-master/scripts/ModifiedEstimateSTR.ForThausandsOfSamplesPerChromosome.V2.py \
                --locus_counts ${sep=' ' fLocusCt} \
                --STR_counts  ${sep=' ' fStrCt} \
                --median_cov ${sep=' ' fMedCov} \
                --model /STRetch-master/scripts/STRcov.model.csv \
                --chrom ${chrom} \
                --out ${output_prefix}
                
            	zcat  ${output_prefix}.${chrom}.part1.STRs.tsv.gz | head -n1 >  ${output_prefix}.${chrom}.STRs.tsv
   		zcat ${output_prefix}.${chrom}.part*.STRs.tsv.gz |  egrep -v "^chrom|XXX"   >>  ${output_prefix}.${chrom}.STRs.tsv
   		gzip ${output_prefix}.${chrom}.STRs.tsv
  }
  runtime {
   	docker:"bharatij/stretch_pipeline:stretch"
    	memory: machine_mem_size + " GB"
    	disks: "local-disk " + disk_size + " SSD"
    	cpu: "1"
    	zones: "us-central1-c us-central1-b"
    	preemptible: preemptible_tries   
    	continueOnReturnCode: [0,1]
  }
  
  output {
    	    File strtsvs = "${output_prefix}.${chrom}.STRs.tsv.gz" 
  }
}

#combine all chr and chunk into single file
task CombineChunkStr {

  Array[File] partStrs
  String output_prefix

  String docker_image
  Int preemptible_tries
  Int addtional_disk_size = 4
  String machine_mem_size = 32 #can be < 5GB for few hundred samples	

  Float outsize = 200 #can be < 2GB for few hundred samples
  Int disk_size = ceil(outsize) + addtional_disk_size
  
  #Calls locus cov
  command <<<
   	zcat ${sep=' ' partStrs} | grep ^chrom |sort|uniq >  ${output_prefix}.AllChr.STRs.tsv
   	zcat ${sep=' ' partStrs} |  grep -v ^chrom  >>  ${output_prefix}.AllChr.STRs.tsv
   	gzip ${output_prefix}.AllChr.STRs.tsv
   	
   >>> 
  
  runtime {
    	docker:"bharatij/stretch_pipeline:stretch"
    	memory: machine_mem_size + " GB"
    	disks: "local-disk " + disk_size + " SSD"
    	cpu: "1"
    	zones: "us-central1-c us-central1-b"
    	preemptible: preemptible_tries   
    	continueOnReturnCode: [0,1]
  }
  
  output {
    	  File allStrs = "${output_prefix}.AllChr.STRs.tsv.gz"
  }
}

#Mutiple correction
task RunMultiTestCorrection {
  
  File allStrs
  String sampleid
  
  String docker_image
  Int preemptible_tries
  Int addtional_disk_size = 4
  String machine_mem_size = 32 #can be < 5GB for few hundred samples	

  Float outsize = 200 #can be < 2GB for few hundred samples
    
  Int disk_size = ceil(outsize) + addtional_disk_size
  
  command <<<
   	echo "Processing: ${sampleid}" 
   	zcat ${allStrs} | head -n1 > /cromwell_root/${sampleid}.STRs.tsv
   	zcat ${allStrs} | grep -w ${sampleid}  >> /cromwell_root/${sampleid}.STRs.tsv	
	
      	PATH=$PATH:/STRetch-master/tools/bin; \ 
        /STRetch-master/tools/bin/python /STRetch-master/scripts/MultipleTestCorrection.py --pvalfile /cromwell_root/${sampleid}.STRs.tsv
        
        gzip /cromwell_root/${sampleid}.FDR.STRs.tsv
  >>>

  runtime {
  	docker:"bharatij/stretch_pipeline:stretch"
    	memory: machine_mem_size + " GB"
    	disks: "local-disk " + disk_size + " SSD"
    	cpu: "1"
    	zones: "us-central1-c us-central1-b"
    	preemptible: preemptible_tries   
    	continueOnReturnCode: [0,1]
  }
  output {
    	    File sampleStr = "${sampleid}.FDR.STRs.tsv.gz"
  }
}

#combine all samples, all STRs single file
task CombineSamples {
  #Command parameters
  Array[File] sampleStrs
  String output_prefix
  # Runtime parameters
  String docker_image
  Int preemptible_tries
  Int addtional_disk_size = 4
  String machine_mem_size = 32 #can be < 5GB for few hundred samples	

  Float outsize = 200 #can be < 2GB for few hundred samples
    
  Int disk_size = ceil(outsize) + addtional_disk_size
  #Calls locus cov
  command <<<
   	zcat ${sep=' ' sampleStrs} | grep  ^chrom | sort | uniq >  ${output_prefix}.AllSample.AllChr.FDR.STRs.tsv
   	zcat ${sep=' ' sampleStrs} |  grep -v ^chrom  >>  ${output_prefix}.AllSample.AllChr.FDR.STRs.tsv
   	gzip ${output_prefix}.AllSample.AllChr.FDR.STRs.tsv
   	
   >>> 

   runtime {
    	docker:"bharatij/stretch_pipeline:stretch"
    	memory: machine_mem_size + " GB"
    	disks: "local-disk " + disk_size + " SSD"
    	cpu: "1"
    	zones: "us-central1-c us-central1-b"
    	preemptible: preemptible_tries   
    	continueOnReturnCode: [0,1]
  }
  
  output {
    	  File mastertable = "${output_prefix}.AllSample.AllChr.FDR.STRs.tsv.gz"
  }
}
