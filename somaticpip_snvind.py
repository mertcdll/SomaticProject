import os
import logging
import subprocess
import json
import argparse
from concurrent.futures import ThreadPoolExecutor


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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
                lrom : bool = False):
        
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
                    f"--tumor-sample {self.tumor_name} --germline-resource {self.germline_resource} --interval-padding 100"
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

class Config:
    def __init__(self, config_file):
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        self.output_folder = config.get("output_folder")
        self.ref_fasta = config.get("ref_fasta")
        self.ref_dict = config.get("ref_dict")
        self.interval_list = config.get("interval_list")
        self.scatter_count = config.get("scatter_count")
        self.tumor_bam = config.get("tumor_bam")
        self.germline_resource = config.get("germline_resource")
        self.index_image = config.get("index_image")
        self.normal_bam = config.get("normal_bam")
        self.panel_of_normals = config.get("panel_of_normals")
        self.variants_for_contamination = config.get("variants_for_contamination")
        self.downsampling_stride = config.get("downsampling_stride")
        self.max_reads_per_alignment_start = config.get("max_reads_per_alignment_start")
        self.max_suspicious_reads_per_alignment_start = config.get("max_suspicious_reads_per_alignment_start")
        self.max_population_af = config.get("max_population_af")
        self.lrom = config.get("lrom", False)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Run Somatic Pipeline with configuration.')
    parser.add_argument('config_file', type=str, help='Path to the JSON configuration file.')

    args = parser.parse_args()

    config = Config(args.config_file)

    pipeline = SomaticPipeline(
        config.output_folder, config.ref_fasta, config.ref_dict, config.interval_list, config.scatter_count, 
        config.tumor_bam, config.germline_resource, config.index_image, config.normal_bam, 
        config.panel_of_normals, config.variants_for_contamination, config.downsampling_stride, 
        config.max_reads_per_alignment_start, config.max_suspicious_reads_per_alignment_start, 
        config.max_population_af, config.lrom )

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







