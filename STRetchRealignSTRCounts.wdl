#Workflow : STRetchPipeline-Step2 : Read realignment and generating STR read counts
workflow STRetchPipeline {
  String pipeline_version = "1.0"
  meta {
	author: "Bharati Jadhav"
    	email: "bharati.jadhav@mssm.edu"
    	description: "STRetch read realignment and generating STR read counts step"	
  }
    File cram_fasta 
    File cram_fasta_fai
    File cramFile
    File craiFile
    File TRF_BED
    String output_prefix = basename(cramFile, ".cram")
    # STRetch STRdecoy reference genome fasta and index files
    File ref_fasta 
    File ref_fasta_fai
    File ref_fasta_amb
    File ref_fasta_ann
    File ref_fasta_bwt
    File ref_fasta_pac
    File ref_fasta_sa 
    File ref_fasta_genome 
    File STR_BED  # STRdecoys.sorted.bed"
    
    Int preemptible_tries = 3
   
    call RunMosDepth { 
    	input: 
          cramFile = cramFile,
          craiFile = craiFile,
          cram_fasta = cram_fasta,
          cram_fasta_fai = cram_fasta_fai,
          output_prefix = output_prefix,
          
      	  preemptible_tries = preemptible_tries
    }
    
    call Runbwa { 
        input: 
          cramFile = cramFile,
          craiFile = craiFile,
          cram_fasta = cram_fasta,
          cram_fasta_fai = cram_fasta_fai,
          ref_fasta = ref_fasta,
          ref_fasta_fai = ref_fasta_fai,
          ref_fasta_amb = ref_fasta_amb,
          ref_fasta_ann = ref_fasta_ann,
          ref_fasta_bwt = ref_fasta_bwt,
          ref_fasta_pac = ref_fasta_pac,
          ref_fasta_sa = ref_fasta_sa,
          TRF_BED = TRF_BED,
          output_prefix = output_prefix,
          
          preemptible_tries = preemptible_tries
    }
    
    call RunSTRcov {
        input: 
          bamOut = Runbwa.bamOut,
          baiOut = Runbwa.baiOut,
          ref_fasta_genome = ref_fasta_genome,
          STR_BED = STR_BED,
          output_prefix = output_prefix,
          
          preemptible_tries = preemptible_tries
    }
     
    call RunLocuscov {
      input: 
          bamOut = Runbwa.bamOut,
          baiOut = Runbwa.baiOut,
          TRF_BED = TRF_BED,
          output_prefix = output_prefix,
          
          preemptible_tries = preemptible_tries
    }
	
    output {
      File medCov = RunMosDepth.medCov
      File strCov = RunSTRcov.strCov
      File locusCov = RunLocuscov.locusCov
    }   
}

#Task Definitions
task RunMosDepth {
    File cramFile
    File craiFile
    File cram_fasta 
    File cram_fasta_fai
    String output_prefix
    Int addtional_disk_size = 1 
    String machine_mem_size = 5	
    
    Int preemptible_tries
    Float distTxt_size = 0.01
    Float ref_size = size(cram_fasta, "GB") 
    Float refidx_size = size(cram_fasta_fai, "GB") 
    Int disk_size = ceil(size(cramFile, "GB") + size(craiFile, "GB") + distTxt_size + ref_size + refidx_size) + addtional_disk_size
    
    command {
  	  /STRetch-master/tools/bin/mosdepth -n -t 8 -f ${cram_fasta} ${output_prefix} ${cramFile} 
	    	  
     	  /STRetch-master/tools/bin/python /STRetch-master/scripts/mosdepth_median.py \
      	    				  --out  ${output_prefix}.median_cov \
      					  ${output_prefix}.mosdepth.global.dist.txt
   }

  runtime {
    docker:"bharatij/stretch_pipeline:stretch"
    memory: machine_mem_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    zones: "us-central1-c us-central1-b"
    preemptible: preemptible_tries   
    continueOnReturnCode: [0,1]
  }
  output {
    File medCov = "${output_prefix}.median_cov"
  }
}

task Runbwa {
    File cramFile
    File craiFile
    File cram_fasta 
    File cram_fasta_fai
    File TRF_BED
    String output_prefix
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_amb
    File ref_fasta_ann
    File ref_fasta_bwt
    File ref_fasta_pac
    File ref_fasta_sa 	
   
    Int addtional_disk_size = 10
    String machine_mem_size = 32
    
    Int preemptible_tries
      
    Float out_size = size(cramFile, "GB") * 3
    Float ref_size = size(ref_fasta, "GB") * 3
    Int disk_size = ceil(out_size + ref_size) + addtional_disk_size
   
   command {
   	  set -o pipefail
      	  java -Xmx16g -Dsamjdk.reference_fasta=${cram_fasta}  \
        	       -jar /STRetch-master/tools/bin/bazam.jar \
                       -pad 5 \
                       -n 6 \
                       -L ${TRF_BED} \
                       -bam ${cramFile} | \
          /STRetch-master/tools/bin/bwa mem -p -M -t 7 -R "@RG\tID:${output_prefix}\tPL:illumina\tPU:NA\tLB:NA\tSM:${output_prefix}"  ${ref_fasta} - | \
          /STRetch-master/tools/bin/samtools view -bSuh - | \
          /STRetch-master/tools/bin/samtools sort -o ${output_prefix}.STRdecoy.bam -T ${output_prefix}.STRdecoy
				                
          /STRetch-master/tools/bin/samtools index ${output_prefix}.STRdecoy.bam

   }  
   
    runtime {
    	docker:"bharatij/stretch_pipeline:stretch"
    	memory: machine_mem_size + " GB"
    	disks: "local-disk " + disk_size + " HDD"
    	cpu: "1"
    	zones: "us-central1-c us-central1-b"
    	preemptible: preemptible_tries
    	continueOnReturnCode: [0,1]
   }
   
    output {
    	File bamOut = "${output_prefix}.STRdecoy.bam"
    	File baiOut = "${output_prefix}.STRdecoy.bam.bai"
    }
}

task RunSTRcov {

    File bamOut
    File baiOut
    String output_prefix
    File ref_fasta_genome 
    File STR_BED
    
    Int addtional_disk_size = 1 
    String machine_mem_size = 4	
    
    Int preemptible_tries
    Float strCov_size = 0.01
    Int disk_size = ceil(size(bamOut, "GB") + size(baiOut, "GB") + size(STR_BED, "GB") + strCov_size) + addtional_disk_size
  
  
    command {
     		/STRetch-master/tools/bin/bedtools coverage -counts -sorted -g ${ref_fasta_genome} -a ${STR_BED} -b ${bamOut} > ${output_prefix}.STRdecoy.STR_counts
           }

    runtime {
    	docker:"bharatij/stretch_pipeline:stretch"
    	memory: machine_mem_size + " GB"
    	disks: "local-disk " + disk_size + " HDD"
    	cpu: "1"
    	zones: "us-central1-c us-central1-b"
    	preemptible: preemptible_tries   
    	continueOnReturnCode: [0,1]
  }
 
  output {
    	  File strCov = "${output_prefix}.STRdecoy.STR_counts"
  }
}

task RunLocuscov {
    File bamOut
    File baiOut
    String output_prefix
    File TRF_BED = "TRF.bed"
    
    Int addtional_disk_size = 1 
    String machine_mem_size = 4	
    
    Int preemptible_tries
    Float strLocus_size = 0.01
    Int disk_size = ceil(size(bamOut, "GB") + size(baiOut, "GB") + strLocus_size) + addtional_disk_size
  
  
    command {
  	  STRPATH=$PATH; PATH=/STRetch-master/tools/bin:$PATH; \
              
  	  /STRetch-master/tools/bin/python  /STRetch-master/scripts/identify_locus.py \
          			           --bam ${bamOut} \
                			   --bed ${TRF_BED} \
                                           --output ${output_prefix}.STRdecoy.locus_counts \
             				   ;PATH=$STRPATH
             
          rm ${bamOut}
          rm ${baiOut}
  }

  runtime {
    	docker:"bharatij/stretch_pipeline:stretch"
    	memory: machine_mem_size + " GB"
    	disks: "local-disk " + disk_size + " HDD"
    	cpu: "1"
    	zones: "us-central1-c us-central1-b"
    	preemptible: preemptible_tries   
    	continueOnReturnCode: [0,1]
  }
  
  output {
    	    File locusCov = "${output_prefix}.STRdecoy.locus_counts"
  }
}
