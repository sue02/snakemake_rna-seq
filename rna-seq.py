import os
#shell.prefix(config)
shell.prefix("set -o pipefail; ")
configfile: "config.yaml"
#print(config['htseq']['counts'])

data="/data"
out_dir="/snakemake/out/"

if not os.path.exists(out_dir+'counts'):
	os.mkdir(out_dir+'counts')

(SAMPLE) = glob_wildcards(data+"/{sample}_1.fastq").sample

#print (SAMPLE)

rule all: 
    input:    expand("singleton/{sample}.fastq", sample=SAMPLE),
              "genome/chrName.txt",
              expand(out_dir+"{sample}/{sample}_sortedbyread.bam", sample=SAMPLE),
              expand(out_dir+"{sample}/{sample}Log.final.out", sample=SAMPLE),
              expand(out_dir+"counts/{sample}.counts", sample=SAMPLE)
	

    

rule subsample:
	input:
		in1=data+"/{sample}_1.fastq", 
		in2=data+"/{sample}_2.fastq"
	output:
		out1="sub_data/{sample}_sub_1.fastq",
		out2="sub_data/{sample}_sub_2.fastq",
		trim1="trimmed/{sample}_trimmed_1.fastq",
		trim2="trimmed/{sample}_trimmed_2.fastq",
		final1="input1/{sample}_1.fastq",
		final2="input1/{sample}_2.fastq",
		single="singleton/{sample}.fastq"
	shell:"""
		#echo {input} {output} 
		ml BBMap
		ml fastx
		reformat.sh in1={input.in1} in2={input.in2} out1={output.out1} out2={output.out2} samplerate=0.001
		fastq_quality_trimmer -Q33 -t 30 -l 50 -i {output.out1} -o {output.trim1} 
		fastq_quality_trimmer -Q33 -t 30 -l 50 -i {output.out2} -o {output.trim2} 
		repair.sh in1={output.trim1} in2={output.trim2} out1={output.final1} out2={output.final2} outsingle={output.single} -ow=t

		"""

rule index_genome:
		input:
			genome="genome/chr1.fa"
		output:"genome/chrName.txt"
		shell:"""
			ml STAR
			STAR --runThreadN 10 --runMode genomeGenerate --genomeDir genome/ --genomeFastaFiles {input.genome} --limitGenomeGenerateRAM 43524399488 
			echo "indexing done"
			"""

rule mapping:
		input:
			in11="input1/{sample}_1.fastq",
			in22="input1/{sample}_2.fastq"
		params:
        		outdir = out_dir + "{sample}"
			
		output: 
        		out =  out_dir + "{sample}/{sample}_sortedbyread.bam",
			out1 = out_dir + "{sample}/{sample}Log.final.out"
		threads:10
		shell:"""
			#echo {input.in11} {output}
			module load STAR
			module load samtools
			STAR  --runThreadN {threads} --runMode alignReads --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical --genomeDir genome --readFilesIn {input.in11} {input.in22} --twopassMode Basic --outFileNamePrefix {params.outdir}/{wildcards.sample}
			samtools view -bS {params.outdir}/{wildcards.sample}Aligned.out.sam > {params.outdir}/{wildcards.sample}.bam
 			samtools sort -n {params.outdir}/{wildcards.sample}.bam -o {params.outdir}/{wildcards.sample}_sortedbyread.bam
  			samtools sort {params.outdir}/{wildcards.sample}.bam -o {params.outdir}/{wildcards.sample}_sortedbycoord.bam
  			samtools index {params.outdir}/{wildcards.sample}_sortedbycoord.bam	
			"""

rule HTSeq:
	input:	
		bam=out_dir +"{sample}/{sample}_sortedbyread.bam",
		counts_dir=config['htseq']['counts'],
		gtf=config['gtf']['annotation'] 
	output:
		out=out_dir+"counts/{sample}.counts"
	shell: """
		ml htseq
    		mkdir {input.counts_dir}
		htseq-count -m union -s reverse -t gene -i gene_id -f bam {input.bam} {input.gtf} > {input.counts_dir}/{wildcards.sample}.counts

		"""
