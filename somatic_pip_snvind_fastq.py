import os
import logging
import json
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class Indexer:
    def __init__(self, reference_fasta:str, algotype:str):
        self.reference_fasta = reference_fasta
        self.algotype = algotype


    def index(self):
        
        logging.info("bwa - Indexing the reference genome")
        command  = f"bwa index -a {self.algotype} {self.reference_fasta}"

        try:
            os.system(command)
            logging.info("bwa - Successfully indexed the reference genome")
        except Exception as e: 
            logging.error(f"bwa - An error occurred while indexing the reference genome: {e}")


class Mapper:
    def __init__(self, reference_fasta: str,
                 sample_id: str, 
                 r1:str, 
                 r2:str, 
                 threads:int, 
                 output_path:str,
                 run_name: str, 
                 use_dragen:bool=False,
                 hash_table:str=None):
        
        self.reference_fasta = reference_fasta
        self.sample_id = sample_id
        self.r1 = r1
        self.r2 = r2
        self.threads = threads
        self.output_path = output_path
        self.run_name = run_name
        self.use_dragen = use_dragen
        self.hash_table = hash_table

    def map(self):


        sam_file = self.sample_id+".sam"
        
        if self.use_dragen:
          logging.info(f"dragen - Mapping {self.r1} and {self.r2} to the reference genome")
          command = (
              f"dragen-os -r {self.hash_table} "
              f"-1 {self.r1} -2 {self.r2} "
              f"--RGSM {self.sample_id} --RGID {self.sample_id} "
              f"--num-threads {self.threads} > {os.path.join(self.output_path, sam_file)}"
          )
          
          try:
            os.system(command)
            logging.info(f"dragen - Successfully mapped {self.r1} and {self.r2} to the reference genome")
          except Exception as e:
            logging.error(f"dragen - An error occurred while mapping {self.r1} and {self.r2}: {e}")

        
        else:
          logging.info(f"bwa - Mapping {self.r1} and {self.r2} to the reference genome")          
          command = (
              f"bwa mem -M -t {self.threads} "
              f"-R '@RG\\tID:{self.sample_id}\\tSM:{self.sample_id}\\tLB:Mylib\\tPU:Illumina' "
              f"{self.reference_fasta} {self.r1} {self.r2} > {os.path.join(self.output_path, sam_file)}"
          )

          try:
              os.system(command)
              logging.info(f"bwa - Successfully mapped {self.r1} and {self.r2} to the reference genome")
          except Exception as e:
              logging.error(f"bwa - An error occurred while mapping {self.r1} and {self.r2}: {e}")

          return os.path.join(self.output_path, sam_file)



class deduplicate:
    def __init__(self, 
                 sam_path: str,
                 sample_id: str,
                 threads: int,
                 use_mark_duplicates: bool=False,
                 remove_all_duplicates:bool=False,
                 remove_sequencing_duplicates:bool=False):
        
        self.sam_path = sam_path
        self.threads = threads
        self.sample_id = sample_id
        self.use_mark_duplicates = use_mark_duplicates
        self.remove_all_duplicates = remove_all_duplicates
        self.remove_sequencing_duplicates = remove_sequencing_duplicates

    def convert_markdedup(self):

        logging.info(f"Samtools - Creating a bam file for {self.sample_id} from sam file")
        bam_file = self.sample_id+".bam"
        
        out_path = os.path.dirname(self.sam_path)

        output_bam = os.path.join(out_path, bam_file)
        converter_cmd = f"samtools view -@ {self.threads} -bS -o {output_bam} {self.sam_path}"

        try: 
            os.system(converter_cmd)
            logging.info(f"Samtools - Successfully created bam file for {self.sample_id}")
        except Exception as e:
            logging.error(f"Samtools - An error occurred while creating bam file for {self.sample_id}")
        
        filtered_bam_file = "filtered_" + self.sample_id + ".bam"
        filtered_bam = os.path.join(out_path, filtered_bam_file)
        
        metrics_file = self.sample_id + "_dup_metrics.txt"
        metrics_path = os.path.join(out_path, metrics_file)

        if self.use_mark_duplicates:
            
            logging.info(f"MarkDuplicatesSpark - Marking and removing duplicates for bam file of {self.sample_id}")
            
            marker_cmd = f"gatk MarkDuplicatesSpark -I {output_bam} -O {filtered_bam} --metrics-file {metrics_path} --spark-master local[{self.threads}] --create-output-bam-index false"
            
            if self.remove_all_duplicates:
                
                marker_cmd += " --remove-all-duplicates"
            
            if self.remove_sequencing_duplicates:
                
                marker_cmd += " --remove-sequencing-duplicates"
                
            try:
                os.system(marker_cmd)
                logging.info(f"MarkDuplicatesSpark - Successfully marked and removed duplicates for {self.sample_id} > {filtered_bam}")
            except Exception as e:
                logging.error(f"An error occurred while marking and removing duplicates for bam file of {self.sample_id}")

        else:
            
            logging.info(f"Sambamba - Marking and removing duplicates for bam file of {self.sample_id}")
            marker_cmd = f"sambamba markdup -t {self.threads} {output_bam} {filtered_bam}"
            
            if self.remove_all_duplicates:
                
                marker_cmd += " --remove-duplicates"

            try:
                os.system(marker_cmd)
                logging.info(f"Sambamba - Successfully marked and removed duplicates for {self.sample_id} > {filtered_bam}")
            except Exception as e:
                logging.error(f"An error occurred while marking and removing duplicates for bam file of {self.sample_id}")

        os.system(f"rm {self.sam_path}")
        os.system(f"rm {output_bam}")

        sorted_bam_file = "sorted_" + filtered_bam_file
        sorted_bam = os.path.join(out_path, sorted_bam_file)

        sort_cmd = f"samtools sort -@ {self.threads} {filtered_bam} -o {sorted_bam}"

        os.system(sort_cmd)

        index_cmd = f"samtools index -@ {self.threads} {sorted_bam}"

        os.system(index_cmd)

        os.system(f"rm {filtered_bam}")
        
        return sorted_bam, out_path


class SomaticPipeline:

    """This class uses GATK somatic workflow to detect somatic mutations"""

    def __init__(self, 
                output_folder: str, 
                ref_fasta: str, 
                ref_dict: str, 
                interval_list: str, 
                scatter_count: int, 
                tumor_bam: str, 
                germline_resource: str, 
                index_image: str, 
                normal_bam: str = None, 
                panel_of_normals: str = None, 
                variants_for_contamination: str = None, 
                downsampling_stride: int = None, 
                max_reads_per_alignment_start: int = None, 
                max_suspicious_reads_per_alignment_start: int = None, 
                max_population_af : float = None,
                lrom : bool = False,
                interval_padding: int = 100):
        
        self.ref_fasta = ref_fasta
        self.ref_dict = ref_dict
        self.interval_list = interval_list
        self.scatter_count = scatter_count
        self.output_folder = output_folder
        self.normal_bam = normal_bam
        self.tumor_bam = tumor_bam
        self.panel_of_normals = panel_of_normals
        self.variants_for_contamination = variants_for_contamination
        self.germline_resource = germline_resource
        self.index_image = index_image
        self.downsampling_stride = downsampling_stride
        self.max_reads_per_alignment_start = max_reads_per_alignment_start
        self.max_suspicious_reads_per_alignment_start = max_suspicious_reads_per_alignment_start
        self.max_population_af = max_population_af
        self.lrom = lrom
        self.interval_padding = interval_padding

        self.scatter_folder = None
        self.variants_folder = None
        self.mutect_command_list = []
        self.getpileupsummariesn_command_list = []
        self.getpileupsummariest_command_list = []
        self.sortvcf_command_list = []
        self.normal_name_out = None
        self.tumor_name_out = None
        self.artifact_priors_path = None
        self.tumor_name = None
        self.normal_name = None
        self.merged_npu_path = None
        self.merged_tpu_path = None
        self.merged_ufvcf_path = None
        self.merged_out_folder = None
        self.merged_vcfstats_path = None
        self.segmentation_table = None
        self.contamination_table = None
        self.flagged_vcf_path = None
        self.filtering_stats_path = None
        self.filtered_vcf_path = None

    def split_intervals(self):

        """Uses gatk split intervals module to split interval list so that use multiple cores while running mutect2""" 

        logging.info("Splitting Interval List")

        self.scatter_folder = os.path.join(self.output_folder, "scattered_intervals")

        os.makedirs(self.scatter_folder, exist_ok=True)
        
        command = (
            f"gatk SplitIntervals -L {self.interval_list} "
            f"-R {self.ref_fasta} --sequence-dictionary {self.ref_dict} "
            f"--scatter-count {self.scatter_count} -O {self.scatter_folder}" 
        )

        try:
            result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"SplitIntervals - An error occurred while splitting the interval list: {stderr_output}")
            else:
                logging.info(f"SplitIntervals - Intervals have been splitted successfully: {stdout_output}")
                   
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while splitting the interval list: {e}")

    def variant_calling(self):

        """Uses mutect2 module to call somatic variants"""

        self.variants_folder = os.path.join(self.output_folder, "mutect2_outputs")

        os.makedirs(self.variants_folder, exist_ok=True)

        interval_lists = os.listdir(self.scatter_folder)

        interval_lists = sorted(interval_lists)

        logging.info("Getting sample names")

        self.tumor_name_out = os.path.join(self.output_folder, "tumor_name.txt")

        self.normal_name_out = os.path.join(self.output_folder, "normal_name.txt")

        try:
            result = subprocess.run(f"gatk GetSampleName -I {self.tumor_bam} -O {self.tumor_name_out}", shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"GetSampleName - An error occurred while getting the tumor sample name: {stderr_output}")
            else:
                logging.info(f"GetSampleName - Tumor sample name has been taken successfully {stdout_output}")

            if self.normal_bam:
                
                result = subprocess.run(f"gatk GetSampleName -I {self.normal_bam} -O {self.normal_name_out}", shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                print(result)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"GetSampleName - An error occurred while getting the normal sample name: {stderr_output}")
                else:
                    logging.info(f"GetSampleName - Normal sample name has been taken successfully {stdout_output}")

                   
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while getting the sample names: {e}")
        
        try:
                        
            with open (self.tumor_name_out, "r") as file:
                self.tumor_name = file.read().strip()
            
            if self.normal_bam:

                with open (self.normal_name_out, "r") as file: 
                    self.normal_name = file.read().strip()

            for item in interval_lists:
                number = item.split("-")[0]
                shard_folder = os.path.join(self.variants_folder, f"shard_{number}")
                
                os.makedirs(shard_folder)
                
                interval_list_path = os.path.join(self.scatter_folder, item)

                vcf_path = os.path.join(shard_folder, "output.vcf")
                f1r2_path = os.path.join(shard_folder, "f1r2.tar.gz")

                mutect_command = (
                    f"gatk Mutect2 -L {interval_list_path} -R {self.ref_fasta} --sequence-dictionary {self.ref_dict} "
                    f"-I {self.tumor_bam} -O {vcf_path} --f1r2-tar-gz {f1r2_path} "
                    f"--tumor-sample {self.tumor_name} --germline-resource {self.germline_resource} --interval-padding {self.interval_padding}"
                )

                if self.normal_bam:
                    mutect_command += f" -I {self.normal_bam} --normal-sample {self.normal_name}"
                if self.downsampling_stride:
                    mutect_command += f" --downsampling-stride {self.downsampling_stride}"
                if self.max_reads_per_alignment_start:
                    mutect_command += f" --max-reads-per-alignment-start {self.max_reads_per_alignment_start}"
                if self.max_suspicious_reads_per_alignment_start:
                    mutect_command += f" --max-suspicious-reads-per-alignment-start {self.max_suspicious_reads_per_alignment_start}"
                if self.max_population_af:
                    mutect_command += f" --max-population-af {self.max_population_af}"
                if self.panel_of_normals:
                    mutect_command += f" --panel-of-normals {self.panel_of_normals}"
            
                self.mutect_command_list.append(mutect_command)
                print(mutect_command)
                normal_pileups_table_path = os.path.join(shard_folder, "normal-pileups.table")
                tumor_pileups_table_path = os.path.join(shard_folder, "tumor-pileups.table")

                
                getpileupsummariesn_command = (
                    f"gatk GetPileupSummaries -R {self.ref_fasta} -I {self.normal_bam} "
                    f"--interval-set-rule INTERSECTION -L {interval_list_path} "
                    f"-V {self.variants_for_contamination} -L {self.variants_for_contamination} "
                    f"-O {normal_pileups_table_path}"
                )


                getpileupsummariest_command = (
                    f"gatk GetPileupSummaries -R {self.ref_fasta} -I {self.tumor_bam} "
                    f"--interval-set-rule INTERSECTION -L {interval_list_path} "
                    f"-V {self.variants_for_contamination} -L {self.variants_for_contamination} "
                    f"-O {tumor_pileups_table_path}"
                )
            
                self.getpileupsummariesn_command_list.append(getpileupsummariesn_command)
                self.getpileupsummariest_command_list.append(getpileupsummariest_command)


                sortvcf_command = (
                    f"gatk SortVcf -I {vcf_path} -O {vcf_path}"
                )
                
                self.sortvcf_command_list.append(sortvcf_command)

            logging.info("Somatic variant calling has been started")

            with ThreadPoolExecutor() as executor:
                futures = [executor.submit(self.run_mutect2_commands, command) for command in self.mutect_command_list]
        
                for future in futures:
                    try:
                        future.result() 
                    except Exception as e:
                        logging.error(f"Error occurred during running Mutect2 commands: {e}")
            
            logging.info("Sorting VCFs")

            with ThreadPoolExecutor() as executor:
                futures = [executor.submit(self.run_sortvcf_commands, command) for command in self.sortvcf_command_list]
                for future in futures:
                    try:
                        future.result() 
                    except Exception as e:
                        logging.error(f"Error occurred during running SortVcf commands: {e}")
            
            if self.normal_bam and self.variants_for_contamination:
                
                logging.info("Getting normal pileup summaries")
                
                with ThreadPoolExecutor() as executor:
                    futures = [executor.submit(self.run_getpileupsummariesn_commands, command) for command in self.getpileupsummariesn_command_list]
                    for future in futures:
                        try:
                            future.result() 
                        except Exception as e:
                            logging.error(f"Error occurred during running GetPileupSummaries (Normal) commands: {e}")
            
            if self.variants_for_contamination:
                
                logging.info("Getting tumor pileup summaries")
                
                with ThreadPoolExecutor() as executor:
                    futures = [executor.submit(self.run_getpileupsummariest_commands, command) for command in self.getpileupsummariest_command_list]
                    for future in futures:
                        try:
                            future.result() 
                        except Exception as e:
                            logging.error(f"Error occurred during running GetPileupSummaries (Tumor) commands: {e}")

        except Exception as e:
            logging.error(f"An internal error occurred while calling variants: {e}")

    def run_mutect2_commands(self, command):

        result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stderr_output = result.stderr.decode('utf-8').strip()
        stdout_output = result.stdout.decode('utf-8').strip()

        if result.returncode != 0:
            logging.error(f"Mutect2 - An error occurred while running command: {stderr_output}")
        else:
            logging.info(f"Mutect2 - Command completed successfully: {stdout_output}")

    def run_getpileupsummariesn_commands(self, command):
        result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stderr_output = result.stderr.decode('utf-8').strip()
        stdout_output = result.stdout.decode('utf-8').strip()

        if result.returncode != 0:
            logging.error(f"GetPileupSummaries (Normal) - An error occurred while running command")
        else:
            logging.info(f"GetPileupSummaries (Normal)- Command completed successfully: {stdout_output}")
        return
    
    def run_getpileupsummariest_commands(self, command):
        result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stderr_output = result.stderr.decode('utf-8').strip()
        stdout_output = result.stdout.decode('utf-8').strip()

        if result.returncode != 0:
            logging.error(f"GetPileupSummaries (Tumor) - An error occurred while running command")
        else:
            logging.info(f"GetPileupSummaries (Tumor) - Command completed successfully: {stdout_output}")
        return
    
    def run_sortvcf_commands(self, command):

        result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stderr_output = result.stderr.decode('utf-8').strip()
        stdout_output = result.stdout.decode('utf-8').strip()

        if result.returncode != 0:
            logging.error(f"SortVcf - An error occurred while running command: {stderr_output}")
        else:
            logging.info(f"SortVcf - Command completed successfully: {stdout_output}")
        return

    def learnreadorientation(self):
        
        if self.lrom:
            f1r2_list = []
            for root, dirs, files in os.walk(self.variants_folder):
                for file in files:
                    if file == "f1r2.tar.gz":
                        f1r2_list.append(os.path.join(root, file))

            line = ""
            for item in f1r2_list:
                line += "-I " + item + " "

            
            self.merged_out_folder = os.path.join(self.output_folder, "merged_outputs")
            os.makedirs(self.merged_out_folder, exist_ok=True)
            self.artifact_priors_path = os.path.join(self.merged_out_folder, "artifact-priors.tar.gz")

            command = f"gatk LearnReadOrientationModel {line.strip()} -O {self.artifact_priors_path}"                               

            logging.info("Starting to create artifact priors")

            try:
                result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"LearnReadOrientationModel - An error occurred while creating artifact priors: {stderr_output}")
                else:
                    logging.info(f"LearnReadOrientationModel - Artifact priors have been created successfully: {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somaticpip.py - An internal error occurred while creating the artifact priors: {e}")

    def mergenormalpileups(self):

        if self.normal_bam and self.variants_for_contamination:
            normal_pu_list = []
            for root,dirs,files in os.walk(self.variants_folder):
                for file in files:
                    if file == "normal-pileups.table":
                        normal_pu_list.append(os.path.join(root, file))

            line = ""
            for item in normal_pu_list:
                line += "-I " + item + " "

            self.merged_out_folder = os.path.join(self.output_folder, "merged_outputs")
            os.makedirs(self.merged_out_folder, exist_ok=True)
            self.merged_npu_path = os.path.join(self.merged_out_folder, f"{self.normal_name}_clean.tsv") 
        
            command = f"gatk GatherPileupSummaries --sequence-dictionary {self.ref_dict} {line.strip()} -O {self.merged_npu_path}"
            
            logging.info("Merging normal pileup summaries")

            try:
                result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"GatherPileupSummaries - An error occurred while gathering normal pileups: {stderr_output}")
                else:
                    logging.info(f"GatherPileupSummaries - Mormal pileups have been created successfully: {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somaticpip.py - An internal error occurred while gathering normal pileups: {e}")

    def mergetumorpileups(self):

        if self.variants_for_contamination:
            tumor_pu_list = []
            for root,dirs,files in os.walk(self.variants_folder):
                for file in files:
                    if file == "tumor-pileups.table":
                        tumor_pu_list.append(os.path.join(root, file))

            line = ""
            for item in tumor_pu_list:
                line += "-I " + item + " "

            self.merged_out_folder = os.path.join(self.output_folder, "merged_outputs")
            os.makedirs(self.merged_out_folder, exist_ok=True)
            self.merged_tpu_path = os.path.join(self.merged_out_folder, f"{self.tumor_name}_clean.tsv") 
        
            command = f"gatk GatherPileupSummaries --sequence-dictionary {self.ref_dict} {line.strip()} -O {self.merged_tpu_path}"
            
            logging.info("Merging tumor pileup summaries")

            try:
                result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"GatherPileupSummaries - An error occurred while gathering tumor pileups: {stderr_output}")
                else:
                    logging.info(f"GatherPileupSummaries - Tumor pileups have been created successfully: {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somaticpip.py - An internal error occurred while gathering tumor pileups: {e}")

    def mergevcfs(self):
        
        vcf_list = []
        for root,dirs,files in os.walk(self.variants_folder):
            for file in files:
                if file == "output.vcf":
                    vcf_list.append(os.path.join(root, file))

        line = ""
        for item in vcf_list:
            line += "-I " + item + " "

        self.merged_out_folder = os.path.join(self.output_folder, "merged_outputs")
        os.makedirs(self.merged_out_folder, exist_ok=True)
        self.merged_ufvcf_path = os.path.join(self.merged_out_folder, f"{self.tumor_name}_clean-unfiltered.vcf.gz") 
    
        command = f"gatk MergeVcfs {line.strip()} -O {self.merged_ufvcf_path}"

        logging.info("Merging VCFs")

        try:
            result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"MergeVcfs - An error occurred while merging VCFs: {stderr_output}")
            else:
                logging.info(f"MergeVcfs - Merged VCF (Unfiltered) created successfully: {stdout_output}")
                
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while merging VCFs: {e}")

    def mergestats(self):

        stats_list = []
        for root,dirs,files in os.walk(self.variants_folder):
            for file in files:
                if file == "output.vcf.stats":
                    stats_list.append(os.path.join(root, file))

        line = ""
        for item in stats_list:
            line += "-stats " + item + " "

            self.merged_out_folder = os.path.join(self.output_folder, "merged_outputs")
            os.makedirs(self.merged_out_folder, exist_ok=True)
            self.merged_vcfstats_path = os.path.join(self.merged_out_folder, "merged.stats") 
    
        command = f"gatk MergeMutectStats {line.strip()} -O {self.merged_vcfstats_path}"

        logging.info("Merging Stats")

        try:
            result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"MergeMutectStats- An error occurred while merging VCF stats: {stderr_output}")
            else:
                logging.info(f"MergeMutectStats - Merged stats created successfully: {stdout_output}")
                
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while merging VCF stats: {e}")


    def calculatecontamination(self):
        
        self.contamination_table = os.path.join(self.merged_out_folder, "contamination.table")
        self.segmentation_table = os.path.join(self.merged_out_folder, "segments.table")
        
        if self.variants_for_contamination:
            command = (
                f"gatk CalculateContamination -I {self.merged_tpu_path} "
                f"-O {self.contamination_table} --tumor-segmentation {self.segmentation_table}"
            )
            
            if self.normal_bam:
                command += f" -matched {self.merged_npu_path}"

            logging.info("Calculating contamination")

            try:
                result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"CalculateContamination- An error occurred while calculating contamination: {stderr_output}")
                else:
                    logging.info(f"CalculateContamination - Contamination data created successfully: {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somaticpip.py - An internal error occurred while calculating contamination: {e}")

    def filter(self):
        
        self.flagged_vcf_path = os.path.join(self.merged_out_folder, f"{self.tumor_name}_clean-flagged.vcf.gz")
        self.filtering_stats_path = os.path.join(self.merged_out_folder, "filtering.stats")
        command = (
            f"gatk FilterMutectCalls -V {self.merged_ufvcf_path} -R {self.ref_fasta} "
            f"-O {self.flagged_vcf_path} --stats {self.merged_vcfstats_path} "
            f"--filtering-stats {self.filtering_stats_path}"
        )

        if self.variants_for_contamination:
            command += f" --contamination-table {self.contamination_table} --tumor-segmentation {self.segmentation_table}"
        if self.lrom:
            command += f" --ob-priors {self.artifact_priors_path}"

        logging.info("Flagging germline variants")

        try:
            result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"FilterMutectCalls- An error occurred while flagging germline variants: {stderr_output}")
            else:
                logging.info(f"FilterMutectCalls - Germline flagged VCF is created successfully: {stdout_output}")
                
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while flagging germline variant calls: {e}")

    def filteralignment(self):
        
        self.filtered_vcf_path = os.path.join(self.merged_out_folder, f"{self.tumor_name}_clean-filtered.vcf.gz")
        self.further_filtered_vcf_path = os.path.join(self.merged_out_folder, f"{self.tumor_name}_clean-further-filtered.vcf.gz")
        
        filter_alignment_command = (
            f"gatk FilterAlignmentArtifacts -R {self.ref_fasta} -V {self.flagged_vcf_path} "
            f"-I {self.tumor_bam} --bwa-mem-index-image {self.index_image} "
            f"-O {self.filtered_vcf_path}"
        )

        logging.info("Filtering germline variants and flagging alignment artifacts")

        try:
            result = subprocess.run(filter_alignment_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"FilterAlignmentArtifacts- An error occurred while filtering germline variants and flagging alignment artifacts: {stderr_output}")
            else:
                logging.info(f"FilterAlignmentArtifacts - Germline filtered, alignment artifact flagged VCF is created successfully: {stdout_output}")
                
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while filtering germline variants and flagging alignment artifacts: {e}")

        further_filter_command = f"bcftools filter -i 'FILTER=\"PASS\"' -O z {self.filtered_vcf_path} -o {self.further_filtered_vcf_path}"

        logging.info("Filtering alignment artifacts")

        try:
            result = subprocess.run(further_filter_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"bcftools filter- An error occurred while filtering alignment artifacts: {stderr_output}")
            else:
                logging.info(f"bcftools filter - Filtered VCF is created successfully: {stdout_output}")
                
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while filtering alignment artifacts: {e}")

        index_command = f"tabix -p vcf {self.further_filtered_vcf_path}"

        logging.info("Indexing the final VCF")

        try:
            result = subprocess.run(index_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"tabix- An error occurred while indexing filtered VCF: {stderr_output}")
            else:
                logging.info(f"tabix - Indexed VCF is created successfully: {stdout_output}")
                
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while indexing filtered VCFs: {e}")


def process_samples(config_path):
    
    with open(config_path, 'r') as config_file:
        config = json.load(config_file)

    ref_fasta = config.get("ref_fasta")
    ref_dict = config.get("ref_dict")
    algotype = config.get("algotype")
    output_path = config.get("output_folder")
    run_name = config.get("run_name")
    threads = config.get("aligner_threads")
    remove_all_duplicates = config.get("remove_all_duplicates")
    remove_sequencing_duplicates = config.get("remove_sequencing_duplicates")
    use_mark_duplicates = config.get("use_gatk_mark_duplicates")
    use_dragen = config.get("use_dragen")
    index_fasta = config.get("index_fasta")
    interval_list = config.get("interval_list")
    scatter_count = config.get("scatter_count")
    tumor_samples = config.get("tumor_samples")
    germline_resource = config.get("germline_resource")
    index_image = config.get("index_image")
    normal_samples = config.get("normal_samples")
    panel_of_normals = config.get("panel_of_normals")
    variants_for_contamination = config.get("variants_for_contamination")
    downsampling_stride = config.get("downsampling_stride")
    max_reads_per_alignment_start = config.get("max_reads_per_alignment_start")
    max_suspicious_reads_per_alignment_start = config.get("max_suspicious_reads_per_alignment_start")
    max_population_af = config.get("max_population_af")
    lrom = config.get("lrom")
    interval_padding = config.get("interval_padding")
    
    if index_fasta:
    
        indexer = Indexer(reference_fasta=ref_fasta, algotype=algotype)
        indexer.index()

    for tumor_sample in tumor_samples: 
        
        tumor_sample_id = tumor_sample.get("sample_id")
        tumor_r1 = tumor_sample.get("r1")
        tumor_r2 = tumor_sample.get("r2")
        
        tumors_path = os.path.join(output_path, run_name, "tumors", tumor_sample_id)
        
        os.makedirs(tumors_path, exist_ok=True)

        mapper = Mapper(reference_fasta=ref_fasta, sample_id=tumor_sample_id, r1=tumor_r1, r2=tumor_r2, output_path=tumors_path, run_name=run_name, threads=threads, use_dragen=use_dragen)
        sam_file = mapper.map()

        deduplicator = deduplicate(sam_path=sam_file, sample_id=tumor_sample_id, threads=threads, remove_all_duplicates=remove_all_duplicates, remove_sequencing_duplicates=remove_sequencing_duplicates, use_mark_duplicates=use_mark_duplicates)
        
        tumor_bam, sample_path = deduplicator.convert_markdedup()
        
        index = next((i for i, sample in enumerate(tumor_samples) if sample["sample_id"] == tumor_sample_id))
        
        normal_bam = None

        if normal_samples:
            
            normal_sample = normal_samples[index]
            
            normal_sample_id = normal_sample.get("sample_id")
            print(normal_sample_id)
            normal_sample_r1 = normal_sample.get("r1")
            normal_sample_r2 = normal_sample.get("r2")
            
            normals_path = os.path.join(output_path, run_name, "normals", normal_sample_id)
            
            os.makedirs(normals_path, exist_ok=True)
                        
            mapper = Mapper(reference_fasta=ref_fasta, sample_id=normal_sample_id, r1=normal_sample_r1, r2=normal_sample_r2, output_path=normals_path, run_name=run_name, threads=threads, use_dragen=use_dragen)
            sam_file = mapper.map()

            deduplicator = deduplicate(sam_path=sam_file, sample_id=normal_sample_id, threads=threads, remove_all_duplicates=remove_all_duplicates, remove_sequencing_duplicates=remove_sequencing_duplicates, use_mark_duplicates=use_mark_duplicates)
            normal_bam, sample_path_normal = deduplicator.convert_markdedup()
            
        
        pipeline = SomaticPipeline(
            sample_path, ref_fasta, ref_dict, interval_list, scatter_count, 
            tumor_bam, germline_resource, index_image, normal_bam, 
            panel_of_normals, variants_for_contamination, downsampling_stride, 
            max_reads_per_alignment_start, max_suspicious_reads_per_alignment_start, 
            max_population_af, lrom, interval_padding)

        pipeline.split_intervals()
        pipeline.variant_calling()
        pipeline.learnreadorientation()
        pipeline.mergenormalpileups()
        pipeline.mergetumorpileups()
        pipeline.mergevcfs()
        pipeline.mergestats()
        pipeline.calculatecontamination()
        pipeline.filter()
        pipeline.filteralignment()


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Run Somatic Pipeline with configuration.')
    parser.add_argument('--config', type=str, help='Path to the JSON configuration file.')

    args = parser.parse_args()

    process_samples(args.config)
    


