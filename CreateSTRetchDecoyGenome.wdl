#WORKFLOW: Generate STRetch Decoy Gennome: STRetchPipeline-Step1
workflow CreateSTRetchDecoyGenome {
  String pipeline_version = "1.0"
  meta {
	author: "Bharati Jadhav"
    	email: "bharati.jadhav@mssm.edu"
    	description: "Generate STRetch reference genome with decoy chromosomes"	
  }
    
    File RefFasta
    File ref_fasta_index
    String fastaName = basename(RefFasta, ".fasta")
    String gotc_docker = "bharatij/stretch_pipeline:stretch"
    Int preemptible_tries = 3
    #Generate STR decoy chromosomes
    call RunSTRDecoyFasta { 
    	input: 
          RefFasta = RefFasta,
          ref_fasta_index = ref_fasta_index,
          fastaName = fastaName,
          docker_image = gotc_docker,
      	  preemptible_tries = preemptible_tries
    }
    output {
      File STRdecoys = RunSTRDecoyFasta.STRdecoys
      File decoyFasta = RunSTRDecoyFasta.decoyFasta
      File sortedFasta = RunSTRDecoyFasta.sortedFasta      
      File bwtFasta = RunSTRDecoyFasta.bwtFasta
      File pacFasta = RunSTRDecoyFasta.pacFasta
      File annFasta = RunSTRDecoyFasta.annFasta
      File ambFasta = RunSTRDecoyFasta.ambFasta
      File saFasta = RunSTRDecoyFasta.saFasta
      File faiFasta = RunSTRDecoyFasta.faiFasta
      File genomeFasta = RunSTRDecoyFasta.genomeFasta
      File strBed = RunSTRDecoyFasta.strBed   
    }
}

#tasks
task RunSTRDecoyFasta {
    File RefFasta
    File ref_fasta_index
    String fastaName
   
    Int addtional_disk_size = 15
    String machine_mem_size = 10
    String docker_image
    Int preemptible_tries
    Int disk_size = ceil(size(RefFasta, "GB") + size(ref_fasta_index, "GB")) + addtional_disk_size
   
    command <<<
   	/STRetch-master/tools/bin/python /STRetch-master/scripts/decoy_STR.py --length 2000 > STRdecoys.fasta                
        cat ${RefFasta} STRdecoys.fasta > ${fastaName}.STRdecoys.fasta                
        /STRetch-master/tools/bin/python /STRetch-master/scripts/sort_fasta.py --infile ${fastaName}.STRdecoys.fasta --outfile ${fastaName}.STRdecoys.sorted.fasta                 
        /STRetch-master/tools/bin/bwa index ${fastaName}.STRdecoys.sorted.fasta                
        /STRetch-master/tools/bin/samtools faidx ${fastaName}.STRdecoys.sorted.fasta
   	grep STR ${fastaName}.STRdecoys.sorted.fasta.fai | awk -v OFS='\t' '{ print $1, 0, $2 -1 }' | /STRetch-master/tools/bin/bedtools sort -i - > STRdecoys.sorted.bed  		              
        cut -f1,2 ${fastaName}.STRdecoys.sorted.fasta.fai > ${fastaName}.STRdecoys.sorted.fasta.genome
    >>>  
      
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
        File STRdecoys = "STRdecoys.fasta"
        File decoyFasta = "${fastaName}.STRdecoys.fasta"
        File sortedFasta = "${fastaName}.STRdecoys.sorted.fasta"
        File bwtFasta = "${fastaName}.STRdecoys.sorted.fasta.bwt"
        File pacFasta = "${fastaName}.STRdecoys.sorted.fasta.pac"
        File annFasta = "${fastaName}.STRdecoys.sorted.fasta.ann"
        File ambFasta = "${fastaName}.STRdecoys.sorted.fasta.amb"
        File saFasta = "${fastaName}.STRdecoys.sorted.fasta.sa"
        File faiFasta = "${fastaName}.STRdecoys.sorted.fasta.fai"
        File genomeFasta = "${fastaName}.STRdecoys.sorted.fasta.genome"
        File strBed = "STRdecoys.sorted.bed"
    }
}
