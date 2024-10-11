import os
tr="input_folder/regions/vamos.motif.hg38.v2.1.e0.1.tsv"
tr_pathogenic="/confidential/tGenVar/Lion/TandemTwist/pathogenic.bed"
run_pathogenic = True
tr_trgt="input_folder/regions/vamos.motif.hg38.v2.1.e0.1_trgt_format.tsv"
ref="input_folder/ref/hg38.fna"
tools_paths={\
"TandemTwist":"/confidential/tGenVar/vntr/output_maryam/tools/TandemTwist/TandemTwister",\
"TRGT":"/confidential/tGenVar/vntr/output_maryam/tools/trgt/trgt-v1.1.1-x86_64-unknown-linux-gnu/trgt",\
"vamos":"/confidential/tGenVar/vntr/output_maryam/tools/vamos/src/vamos"}


#/confidential/tGenVar/Lion/TandemTwist/TandemTwister
#/confidential/tGenVar/Lion/TandemTwist/TandemTwister

(datasets, techs, samples) = glob_wildcards("input_folder/samples/{dataset}/{tech}/{sample}.bam.bai")

# datasets can be trios or populations such hgsvc, ...

print(samples)
print(datasets)


datasets_to_samples = {dataset:[] for dataset in set(datasets)}
for i in range(len(datasets)):
    if samples[i] in datasets_to_samples[datasets[i]]:
        continue
    datasets_to_samples[datasets[i]].append(samples[i])


datasets_to_tech = {dataset:[] for dataset in set(datasets)}
for i in range(len(datasets)):
    if techs[i] not in datasets_to_tech[datasets[i]]:
        datasets_to_tech[ datasets[i] ].append(techs[i])


tools = list(tools_paths.keys())
techs = list(set(techs))

def get_outputs(datasets_to_samples, datasets_to_tech):
    tools_list = tools.copy()
    tools_minus_trgt = list(set(tools_list)-set(["TRGT"]))
    extended_outputs = []
    
    for dataset, samples in datasets_to_samples.items():
        
        print (f"dataset = {dataset}, sample = {samples}")
        techs_list = datasets_to_tech[dataset].copy()
        techs_ccs = list(set(techs_list) & set(["CCS"]))
        techs_minus_ccs = list(set(techs_list) - set(["CCS"]))
        
        extended_outputs.extend( expand('output/{dataset}/{tool}/{tech}/{sample}.vcf.gz', \
            dataset=dataset, tool = tools_minus_trgt, tech = techs_list, sample = samples))
        extended_outputs.extend( expand('output/{dataset}/{tool}/{tech}/{sample}.vcf.gz.tbi', \
            dataset=dataset, tool = tools_minus_trgt, tech = techs_list, sample = samples)), \

        extended_outputs.extend( expand('output/{dataset}/TRGT/{tech}/{sample}.vcf.gz', \
            dataset=dataset, tech = techs_ccs, sample = samples))
        extended_outputs.extend( expand('output/{dataset}/TRGT/{tech}/{sample}.vcf.gz.tbi', \
            dataset=dataset, tech = techs_ccs, sample = samples))
        
        print(f"output: {extended_outputs}")

    return extended_outputs

# sample sexes are too be defined here. female:0, male:1
def get_sex(sample):
    sex_map = {
        'HG002_head_1000': 1,
        'HG003_head_1000': 1,
        'HG004_head_1000': 0,
        'HG002': 1,
        'HG003': 1,
        'HG004': 0,
        'HG00731': 1,
        'HG00732': 0,
        'HG00733': 0,
        'NA19238': 0,
        'NA19239': 1,
        'NA19240': 0
    }
    if sample not in sex_map:
        # raise ValueError(f"Sample {sample} not found in sex_map")
        sex_map[sample] = 0 # change it later
    return sex_map[sample]

print("************************************************")
print(get_outputs(datasets_to_samples, datasets_to_tech))

rule all:
    input:
        get_outputs(datasets_to_samples, datasets_to_tech)


rule reformat_repeats_to_trgt:
    input:
        tr=tr,
    output:
        tr_trgt=tr_trgt
    log:
        'logs/reformat_repeats_to_trgt.log'
    shell:
        '''
        (time python scripts/reformat_repeats_for_trgt.py -i {input} -o {output}) > {log} 2>&1
        '''

rule TandemTwist:
    input:
        bam="input_folder/samples/{datasets}/{tech}/{sample}.bam",
        bai="input_folder/samples/{datasets}/{tech}/{sample}.bam.bai",
        tr=tr,
        tr_pathogenic=tr_pathogenic,
        ref=ref
    output:
        vcf='output/{datasets}/TandemTwist/{tech}/{sample}.vcf.gz'
    log:
        'logs/TR_type_{datasets}_TandemTwist_{tech}_{sample}.log'
    benchmark:
        'benchmarks/TR_type_{datasets}_TandemTwist_{tech}_{sample}.benchmark'
    threads: lambda wildcards: 1 if wildcards.tech == "asm" else 32, 
    params:
        tool = tools_paths["TandemTwist"],
        outdir = lambda wildcards, output: os.path.dirname(output.vcf),
        sex = lambda wildcards: get_sex(wildcards.sample),
        input_type = lambda wildcards: "assembly" if wildcards.tech == "asm" else "reads",
        read_type = lambda wildcards: "" if wildcards.tech == "asm" else "-rt "+ wildcards.tech,
        pathogenic = run_pathogenic
    run:
        shell(f"(time {params.tool} -b {input.bam} -m {input.tr} -r {input.ref} \
        -o {params.outdir}/{wildcards.sample} -sn {wildcards.sample} -s {params.sex} -bt {params.input_type} {params.read_type} \
        -t {threads}) > {log} 2>&1")
        if parmas.pathogenic:
            shell(f"(time {params.tool} -b {input.bam} -m {input.tr} -r {input.ref} \
            -o {params.outdir}/{wildcards.sample}_pathogenic -sn {wildcards.sample} -s {params.sex} -bt {params.input_type} {params.read_type} \
            -t {threads}) > {log} 2>&1")
        

rule TRGT:
    input:
        bam="input_folder/samples/{dataset}/{tech}/{sample}.bam",
        bai="input_folder/samples/{dataset}/{tech}/{sample}.bam.bai",
        tr=tr_trgt,
        ref=ref
    output:
        vcf='output/{dataset}/TRGT/{tech}/{sample}.vcf.gz'
    log:
        'logs/TR_type_{dataset}_TRGT_{tech}_{sample}.log'
    benchmark:
        'benchmarks/TR_type_{dataset}_TRGT_{tech}_{sample}.benchmark'
    threads: lambda wildcards: 4 if wildcards.tech == "asm" else 32
    params:
        tool= tools_paths["TRGT"],
        outdir = lambda wildcards, output: os.path.dirname(output.vcf),
    shell:
        '''
        ({params.tool} genotype --genome {input.ref} --repeats {input.tr} --reads {input.bam} \
        --output-prefix {params.outdir}/{wildcards.sample} --threads {threads}) > {log} 2>&1
        '''


rule vamos:
    input:
        bam="input_folder/samples/{dataset}/{tech}/{sample}.bam",
        bai="input_folder/samples/{dataset}/{tech}/{sample}.bam.bai",
        tr=tr,
        ref=ref
    output:
        vcf=temp('output/{dataset}/vamos/{tech}/{sample}.vcf'),
        vcfgz='output/{dataset}/vamos/{tech}/{sample}.vcf.gz'
    log:
        'logs/TR_type_{dataset}_vamos_{tech}_{sample}.log'
    benchmark:
        'benchmarks/TR_type_{dataset}_vamos_{tech}_{sample}.benchmark'
    threads: lambda wildcards: 4 if wildcards.tech == "asm" else 32
    params:
        tool= tools_paths["vamos"],
        subcommand = lambda wildcards: "--contig" if wildcards.tech == "asm" else "--read"
    shell:
        '''
        (time {params.tool} {params.subcommand} -b {input.bam} -r {input.tr} -o {output.vcf} -s {wildcards.sample} -t {threads} \
                && bgzip -c {output.vcf} > {output.vcfgz} ) > {log} 2>&1
        '''

rule index_vcf_files:
    input:
        vcf='output/{dataset}/{tool}/{tech}/{sample}.vcf.gz'
    output:
        index='output/{dataset}/{tool}/{tech}/{sample}.vcf.gz.tbi'
    log:
        'logs/index_vcf_files_{dataset}_{tool}_{tech}_{sample}.log'
    benchmark:
        'benchmarks/index_vcf_files_{dataset}_{tool}_{tech}_{sample}.benchmark'
    shell:
        '''
        (time tabix -p vcf -c {input.vcf}) > {log} 2>&1
        '''

