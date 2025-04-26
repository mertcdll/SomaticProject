import os
import logging
import subprocess
import json
import argparse
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional
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
        
        return sorted_bam
      
class SomaticPipelineCNV:
    def __init__(
            self,
            run_folder: str,
            tumor_output_folder: str,
            tumor_name : str,
            ref_fasta: str,
            ref_dict: str,
            interval_list: str,
            common_sites: str,
            tumor_bam: str,
            pon: str,
            normal_output_folder: str = None,
            normal_name: str = None,
            normal_bam: str = None,
            blacklist_intervals: str = None,
            minimum_base_quality: str = None,
            number_of_eigensamples: int = None,
            minimum_total_allele_count_case: int = None,
            minimum_total_allele_count_normal: int = None,
            genotyping_homozygous_log_ratio_threshold: float = None,
            genotyping_base_error_rate: float = None,
            maximum_number_of_segments_per_chromosome:int = None,
            kernel_variance_copy_ratio: float = None,
            kernel_variance_allele_fraction: float = None,
            kernel_scaling_allele_fraction: float = None,
            kernel_approximation_dimension: float = None,
            window_size: List[int] = None,
            number_of_changepoints_penalty_factor: float = None,
            minor_allele_fraction_prior_alpha: float = None,
            number_of_samples_copy_ratio: int = None,
            number_of_burn_in_samples_copy_ratio: int = None,
            number_of_samples_allele_fraction: int = None,
            number_of_burn_in_samples_allele_fraction: int = None,
            smoothing_credible_interval_threshold_copy_ratio: float = None,
            smoothing_credible_interval_threshold_allele_fraction: float = None,
            maximum_number_of_smoothing_iterations: int = None,
            number_of_smoothing_iterations_per_fit: int = None,
            neutral_segment_copy_ratio_lower_bound: float = None,
            neutral_segment_copy_ratio_upper_bound: float = None,
            outlier_neutral_segment_copy_ratio_z_score_threshold: float = None,
            calling_copy_ratio_z_score_threshold: float = None,
            minimum_contig_length: int = None,
            maximum_copy_ratio: float = None,
            point_size_copy_ratio: float = None,
            point_size_allele_fraction: float = None,
            padding: int = None,
            bin_length: int = None):
        
        self.run_folder = run_folder
        self.tumor_output_folder = tumor_output_folder,
        self.normal_output_folder = normal_output_folder
        self.tumor_name = tumor_name
        self.normal_name = normal_name
        self.ref_fasta = ref_fasta
        self.ref_dict = ref_dict
        self.interval_list = interval_list
        self.common_sites = common_sites
        self.tumor_bam = tumor_bam
        self.pon = pon
        self.normal_bam = normal_bam
        self.blacklist_intervals = blacklist_intervals
        self.minimum_base_quality = minimum_base_quality
        self.number_of_eigensamples = number_of_eigensamples
        self.minimum_total_allele_count_case = minimum_total_allele_count_case
        self.minimum_total_allele_count_normal = minimum_total_allele_count_normal
        self.genotyping_homozygous_log_ratio_threshold = genotyping_homozygous_log_ratio_threshold
        self.genotyping_base_error_rate = genotyping_base_error_rate
        self.maximum_number_of_segments_per_chromosome = maximum_number_of_segments_per_chromosome
        self.kernel_variance_copy_ratio = kernel_variance_copy_ratio
        self.kernel_variance_allele_fraction = kernel_variance_allele_fraction
        self.kernel_scaling_allele_fraction = kernel_scaling_allele_fraction
        self.kernel_approximation_dimension = kernel_approximation_dimension
        self.window_size = window_size
        self.number_of_changepoints_penalty_factor = number_of_changepoints_penalty_factor
        self.minor_allele_fraction_prior_alpha = minor_allele_fraction_prior_alpha
        self.number_of_samples_copy_ratio = number_of_samples_copy_ratio
        self.number_of_burn_in_samples_copy_ratio = number_of_burn_in_samples_copy_ratio
        self.number_of_samples_allele_fraction = number_of_samples_allele_fraction
        self.number_of_burn_in_samples_allele_fraction = number_of_burn_in_samples_allele_fraction
        self.smoothing_credible_interval_threshold_copy_ratio = smoothing_credible_interval_threshold_copy_ratio
        self.smoothing_credible_interval_threshold_allele_fraction = smoothing_credible_interval_threshold_allele_fraction
        self.maximum_number_of_smoothing_iterations = maximum_number_of_smoothing_iterations
        self.number_of_smoothing_iterations_per_fit = number_of_smoothing_iterations_per_fit
        self.neutral_segment_copy_ratio_lower_bound = neutral_segment_copy_ratio_lower_bound
        self.neutral_segment_copy_ratio_upper_bound = neutral_segment_copy_ratio_upper_bound
        self.outlier_neutral_segment_copy_ratio_z_score_threshold = outlier_neutral_segment_copy_ratio_z_score_threshold
        self.calling_copy_ratio_z_score_threshold = calling_copy_ratio_z_score_threshold
        self.minimum_contig_length = minimum_contig_length
        self.maximum_copy_ratio = maximum_copy_ratio
        self.point_size_copy_ratio = point_size_copy_ratio
        self.point_size_allele_fraction = point_size_allele_fraction
        self.padding = padding
        self.bin_length = bin_length


        self.tumor_counts_folder = None
        self.normal_counts_folder = None
        self.tumor_counts = None
        self.normal_counts = None
        self.tumor_alcounts = None
        self.normal_alcounts = None
        self.tumor_standardized = None
        self.tumor_denoised = None
        self.normal_denoised = None
        self.normal_standardized = None
        self.tsg_folder = None
        self.nsg_folder = None
        self.cr_seg_tumor = None
        self.called_cr_seg_tumor = None
        self.cr_seg_normal = None
        self.called_cr_seg_normal = None
        self.plots_folder = None
        self.normal_segments = None
        self.tumor_segments = None
        self.thets = None
        self.nhets = None
        self.preprocessed_intervals = None

    def preprocess_intervals(self):

        logging.info("Processing interval list")
        
        self.intervals_name = os.path.basename(self.interval_list).split("interval_list")[0]
        preprocessed_name = self.intervals_name + "preprocessed.interval_list"
        self.preprocessed_intervals = os.path.join(self.run_folder, preprocessed_name)

        command = (
            f" gatk PreprocessIntervals -L {self.interval_list} "
            f"--reference {self.ref_fasta} "
            f"--interval-merging-rule OVERLAPPING_ONLY "
            f"--output {self.preprocessed_intervals}"
        )

        if self.blacklist_intervals:
            command += f" --XL {self.blacklist_intervals}"
        if self.padding:
            command += f" --padding {self.padding}"
        if self.bin_length:
            command += f" --bin-length {self.bin_length}"


        try:
            result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"PreprocessIntervals - An error occurred while processing the interval list: {stderr_output}")
            else:
                logging.info(f"PreprocessIntervals - Interval list have been processed successfully: {stdout_output}")
                   
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while processing the interval list: {e}")
            
    def collect_counts(self):

        logging.info("Collecting tumor read counts")
        
        self.tumor_counts_folder = os.path.join(self.tumor_output_folder, "counts")
        
        os.makedirs(self.tumor_counts_folder, exist_ok=True)

        self.tumor_counts = os.path.join(self.tumor_counts_folder, f"{self.tumor_name}.counts.hdf5")

        tumor_command = (
            f"gatk CollectReadCounts -L {self.preprocessed_intervals} "
            f"--input {self.tumor_bam} "
            f"--reference {self.ref_fasta} "
            f"--format HDF5 --interval-merging-rule OVERLAPPING_ONLY "
            f"--output {self.tumor_counts}" 
        )
        print(tumor_command)
        try:
            result = subprocess.run(tumor_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"CollectReadCounts - An error occurred while counting the reads from tumor sample: {stderr_output}")
            else:
                logging.info(f"CollectReadCounts - Reads have been counted successfully (tumor): {stdout_output}")
                   
        except Exception as e:
            logging.error(f"somatic_pip_cnv.py - An internal error occurred while counting the reads from tumor sample: {e}")


        if self.normal_bam:
            
            logging.info("Collecting normal read counts")
            
            self.normal_counts_folder = os.path.join(self.normal_output_folder, "counts")

            self.normal_counts = os.path.join(self.normal_counts_folder, f"{self.normal_name}.counts.hdf5")


            normal_command = (
                f"gatk CollectReadCounts -L {self.preprocessed_intervals} "
                f"--input {self.normal_bam} "
                f"--reference {self.ref_fasta} "
                f"--format HDF5 --interval-merging-rule OVERLAPPING_ONLY "
                f"--output {self.normal_counts}" 
            )
            print(normal_command)
            try:
                result = subprocess.run(normal_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"CollectReadCounts - An error occurred while counting the reads from tumor sample: {stderr_output}")
                else:
                    logging.info(f"CollectReadCounts - Reads have been counted successfully (tumor): {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somatic_pip_cnv.py - An internal error occurred while counting the reads from tumor sample: {e}")


    def collect_allelic_counts(self):

        logging.info("Collecting tumor allelic counts")
        
        self.tumor_alcounts = os.path.join(self.tumor_counts_folder, f"{self.tumor_name}.allelicCounts.tsv")

        tumor_command = (
            f"gatk CollectAllelicCounts -L {self.common_sites} "
            f"--input {self.tumor_bam} "
            f"--reference {self.ref_fasta} "
            f"--output {self.tumor_alcounts}" 
        )
        
        if self.minimum_base_quality:
            tumor_command += f" --minimum-base-quality {self.minimum_base_quality}"

        print(tumor_command)
        try:
            result = subprocess.run(tumor_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"CollectAllelicCounts - An error occurred while getting allelic counts from tumor sample: {stderr_output}")
            else:
                logging.info(f"CollectAllelicCounts - Allelic counts created successfully (tumor): {stdout_output}")
                   
        except Exception as e:
            logging.error(f"somatic_pip_cnv.py - An internal error occurred while getting allelic counts from tumor sample: {e}")


        if self.normal_bam:
            
            logging.info("Collecting normal allelic counts")

            self.normal_alcounts = os.path.join(self.normal_counts_folder, f"{self.normal_name}.allelicCounts.tsv")

            normal_command = (
                f"gatk CollectAllelicCounts -L {self.common_sites} "
                f"--input {self.normal_bam} "
                f"--reference {self.ref_fasta} "
                f"--output {self.normal_alcounts}" 
            )
            
            if self.minimum_base_quality:
                normal_command += f" --minimum-base-quality {self.minimum_base_quality}"

            print(normal_command)
            try:
                result = subprocess.run(normal_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"CollectAllelicCounts - An error occurred while getting allelic counts from normal sample: {stderr_output}")
                else:
                    logging.info(f"CollectAllelicCounts - Allelic counts created successfully (normal): {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somatic_pip_cnv.py - An internal error occurred while getting allelic counts from normal sample: {e}")
                
    
    def denoise_counts(self):
        
        logging.info("Denoising tumor read counts")
        
        self.tumor_standardized = os.path.join(self.tumor_counts_folder, f"{self.tumor_name}.standardizedCR.tsv")
        self.tumor_denoised = os.path.join(self.tumor_counts_folder, f"{self.tumor_name}.denoisedCR.tsv")

        tumor_command = (
            f"gatk DenoiseReadCounts --input {self.tumor_counts} "
            f"--count-panel-of-normals {self.pon} "
            f"--standardized-copy-ratios {self.tumor_standardized} "
            f"--denoised-copy-ratios {self.tumor_denoised}"
        )
        print(tumor_command)
        try:
            result = subprocess.run(tumor_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"DenoiseReadCounts - An error occurred while denoising read counts of tumor sample: {stderr_output}")
            else:
                logging.info(f"DenoiseReadCounts - Denoised and standardized copy ratios created successfully (tumor): {stdout_output}")
                   
        except Exception as e:
            logging.error(f"somatic_pip_cnv.py - An internal error occurred while denoising read counts of tumor sample: {e}")


        if self.normal_bam:

            logging.info("Denoising normal read counts")
            
            self.normal_standardized = os.path.join(self.normal_counts_folder, f"{self.normal_name}.standardizedCR.tsv")
            self.normal_denoised = os.path.join(self.normal_counts_folder, f"{self.normal_name}.denoisedCR.tsv")

            normal_command = (
                f"gatk DenoiseReadCounts --input {self.normal_counts} "
                f"--count-panel-of-normals {self.pon} "
                f"--standardized-copy-ratios {self.normal_standardized} "
                f"--denoised-copy-ratios {self.normal_denoised}"
            )
            print(normal_command)
            try:
                result = subprocess.run(normal_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"DenoiseReadCounts - An error occurred while denoising read counts of normal sample: {stderr_output}")
                else:
                    logging.info(f"DenoiseReadCounts - Denoised and standardized copy ratios created successfully (normal): {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somatic_pip_cnv.py - An internal error occurred while denoising read counts of normal sample: {e}")
                
    def model_segments(self):
        
        logging.info("Modelling segments - Tumor")
        
        self.tsg_folder = os.path.join(self.tumor_output_folder, "tumor_segments")
        
        os.makedirs(self.tsg_folder, exist_ok=True)
        
        tumor_command = (
            f"gatk ModelSegments --denoised-copy-ratios {self.tumor_denoised} "
            f"--allelic-counts {self.tumor_alcounts} "
            f"--output {self.tsg_folder} "
            f"--output-prefix {self.tumor_name}"
        )
        
        if self.normal_bam:
            tumor_command += f" --normal-allelic-counts {self.normal_alcounts}"
        if self.minimum_total_allele_count_case:
            tumor_command += f" --minimum-total-allele-count-case {self.minimum_total_allele_count_case}"
        if self.minimum_total_allele_count_normal:
            tumor_command += f" --minimum-total-allele-count-normal {self.minimum_total_allele_count_normal}"
        if self.genotyping_homozygous_log_ratio_threshold:
            tumor_command += f" --genotyping-homozygous-log-ratio-threshold {self.genotyping_homozygous_log_ratio_threshold}"
        if self.genotyping_base_error_rate:
            tumor_command += f" --genotyping-base-error-rate {self.genotyping_base_error_rate}"
        if self.maximum_number_of_segments_per_chromosome:
            tumor_command += f" --maximum-number-of-segments-per-chromosome {self.maximum_number_of_segments_per_chromosome}"
        if self.kernel_variance_copy_ratio:
            tumor_command += f" --kernel-variance-copy-ratio {self.kernel_variance_copy_ratio}"
        if self.kernel_variance_allele_fraction:
            tumor_command += f" --kernel-variance-allele-fraction {self.kernel_variance_allele_fraction}"
        if self.kernel_scaling_allele_fraction:
            tumor_command += f" --kernel-scaling-allele-fraction {self.kernel_scaling_allele_fraction}"
        if self.kernel_approximation_dimension:
            tumor_command += f" --kernel-approximation-dimension {self.kernel_approximation_dimension}"
        if self.number_of_changepoints_penalty_factor:
            tumor_command += f" --number-of-changepoints-penalty-factor {self.number_of_changepoints_penalty_factor}"
        if self.minor_allele_fraction_prior_alpha:
            tumor_command += f" --minor-allele-fraction-prior-alpha {self.minor_allele_fraction_prior_alpha}"
        if self.number_of_samples_copy_ratio:
            tumor_command += f" --number-of-samples-copy-ratio {self.number_of_samples_copy_ratio}"
        if self.number_of_burn_in_samples_copy_ratio:
            tumor_command += f" --number-of-burn-in-samples-copy-ratio {self.number_of_burn_in_samples_copy_ratio}"
        if self.number_of_samples_allele_fraction:
            tumor_command += f" --number-of-samples-allele-fraction {self.number_of_samples_allele_fraction}"
        if self.number_of_burn_in_samples_allele_fraction:
            tumor_command += f" --number-of-burn-in-samples-allele-fraction {self.number_of_burn_in_samples_allele_fraction}"
        if self.smoothing_credible_interval_threshold_copy_ratio:
            tumor_command += f" --smoothing-credible-interval-threshold-copy-ratio {self.smoothing_credible_interval_threshold_copy_ratio}"
        if self.smoothing_credible_interval_threshold_allele_fraction:
            tumor_command += f" --smoothing-credible-interval-threshold-allele-fraction {self.smoothing_credible_interval_threshold_allele_fraction}"
        if self.maximum_number_of_smoothing_iterations:
            tumor_command += f" --maximum-number-of-smoothing-iterations {self.maximum_number_of_smoothing_iterations}"
        if self.number_of_smoothing_iterations_per_fit:
            tumor_command += f" --number-of-smoothing-iterations-per-fit {self.number_of_smoothing_iterations_per_fit}"

        if self.window_size:
            for size in self.window_size:
                tumor_command += f" --window-size {size}"
        print(tumor_command)
        try:
            result = subprocess.run(tumor_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"ModelSegments- An error occurred while modelling segments of tumor sample: {stderr_output}")
            else:
                logging.info(f"ModelSegments - Segments created successfully (tumor): {stdout_output}")
                   
        except Exception as e:
            logging.error(f"somatic_pip_cnv.py - An internal error occurred while modelling segments of tumor sample: {e}")
            
            
        if self.normal_bam:
            
            logging.info("Modelling segments - Normal")
            
            self.nsg_folder = os.path.join(self.normal_output_folder, "normal_segments")
            
            os.makedirs(self.nsg_folder, exist_ok=True)
            
            normal_command = (
                f"gatk ModelSegments --denoised-copy-ratios {self.normal_denoised} "
                f"--allelic-counts {self.normal_alcounts} "
                f"--output {self.nsg_folder} "
                f"--output-prefix {self.normal_name}"
            )
            
            if self.minimum_total_allele_count_case:
                normal_command += f" --minimum-total-allele-count-case {self.minimum_total_allele_count_case}"
            if self.minimum_total_allele_count_normal:
                normal_command += f" --minimum-total-allele-count-normal {self.minimum_total_allele_count_normal}"
            if self.genotyping_homozygous_log_ratio_threshold:
                normal_command += f" --genotyping-homozygous-log-ratio-threshold {self.genotyping_homozygous_log_ratio_threshold}"
            if self.genotyping_base_error_rate:
                normal_command += f" --genotyping-base-error-rate {self.genotyping_base_error_rate}"
            if self.maximum_number_of_segments_per_chromosome:
                normal_command += f" --maximum-number-of-segments-per-chromosome {self.maximum_number_of_segments_per_chromosome}"
            if self.kernel_variance_copy_ratio:
                normal_command += f" --kernel-variance-copy-ratio {self.kernel_variance_copy_ratio}"
            if self.kernel_variance_allele_fraction:
                normal_command += f" --kernel-variance-allele-fraction {self.kernel_variance_allele_fraction}"
            if self.kernel_scaling_allele_fraction:
                normal_command += f" --kernel-scaling-allele-fraction {self.kernel_scaling_allele_fraction}"
            if self.kernel_approximation_dimension:
                normal_command += f" --kernel-approximation-dimension {self.kernel_approximation_dimension}"
            if self.number_of_changepoints_penalty_factor:
                normal_command += f" --number-of-changepoints-penalty-factor {self.number_of_changepoints_penalty_factor}"
            if self.minor_allele_fraction_prior_alpha:
                normal_command += f" --minor-allele-fraction-prior-alpha {self.minor_allele_fraction_prior_alpha}"
            if self.number_of_samples_copy_ratio:
                normal_command += f" --number-of-samples-copy-ratio {self.number_of_samples_copy_ratio}"
            if self.number_of_burn_in_samples_copy_ratio:
                normal_command += f" --number-of-burn-in-samples-copy-ratio {self.number_of_burn_in_samples_copy_ratio}"
            if self.number_of_samples_allele_fraction:
                normal_command += f" --number-of-samples-allele-fraction {self.number_of_samples_allele_fraction}"
            if self.number_of_burn_in_samples_allele_fraction:
                normal_command += f" --number-of-burn-in-samples-allele-fraction {self.number_of_burn_in_samples_allele_fraction}"
            if self.smoothing_credible_interval_threshold_copy_ratio:
                normal_command += f" --smoothing-credible-interval-threshold-copy-ratio {self.smoothing_credible_interval_threshold_copy_ratio}"
            if self.smoothing_credible_interval_threshold_allele_fraction:
                normal_command += f" --smoothing-credible-interval-threshold-allele-fraction {self.smoothing_credible_interval_threshold_allele_fraction}"
            if self.maximum_number_of_smoothing_iterations:
                normal_command += f" --maximum-number-of-smoothing-iterations {self.maximum_number_of_smoothing_iterations}"
            if self.number_of_smoothing_iterations_per_fit:
                normal_command += f" --number-of-smoothing-iterations-per-fit {self.number_of_smoothing_iterations_per_fit}"

            if self.window_size:
                for size in self.window_size:
                    normal_command += f" --window-size {size}"
            print(normal_command)
            try:
                result = subprocess.run(normal_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"ModelSegments- An error occurred while modelling segments of normal sample: {stderr_output}")
                else:
                    logging.info(f"ModelSegments - Segments created successfully (normal): {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somatic_pip_cnv.py - An internal error occurred while modelling segments of normal sample: {e}")
                
                
    def call_crs(self):
        
        logging.info("Calling copy ratio segments - Tumor")
        
        self.cr_seg_tumor = os.path.join(self.tsg_folder, f"{self.tumor_name}.cr.seg")
        self.called_cr_seg_tumor = os.path.join(self.tsg_folder, f"{self.tumor_name}.called.cr.seg")
        
        tumor_command = (
            f"gatk CallCopyRatioSegments --input {self.cr_seg_tumor} "
            f"--output {self.called_cr_seg_tumor}"
        )
        
        
        if self.neutral_segment_copy_ratio_lower_bound:
            tumor_command += f" --neutral-segment-copy-ratio-lower-bound {self.neutral_segment_copy_ratio_lower_bound}"
        if self.neutral_segment_copy_ratio_upper_bound:
            tumor_command += f" --neutral-segment-copy-ratio-upper-bound {self.neutral_segment_copy_ratio_upper_bound}"
        if self.outlier_neutral_segment_copy_ratio_z_score_threshold:
            tumor_command += f" --outlier-neutral-segment-copy-ratio-z-score-threshold {self.outlier_neutral_segment_copy_ratio_z_score_threshold}"
        if self.calling_copy_ratio_z_score_threshold:
            tumor_command += f" --calling-copy-ratio-z-score-threshold {self.calling_copy_ratio_z_score_threshold}"
        print(tumor_command)
        try:
            result = subprocess.run(tumor_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"CallCopyRatioSegments- An error occurred while calling copy ratio segments of tumor sample: {stderr_output}")
            else:
                logging.info(f"CallCopyRatioSegments - Copy ratio segments called successfully (tumor): {stdout_output}")
                
        except Exception as e:
            logging.error(f"somatic_pip_cnv.py - An internal error occurred while calling copy ratio segments of tumor sample: {e}")            
            
        if self.normal_bam:
            
            logging.info("Calling copy ratio segments - Normal")
            
            self.cr_seg_normal = os.path.join(self.nsg_folder, f"{self.normal_name}.cr.seg")
            self.called_cr_seg_normal = os.path.join(self.nsg_folder, f"{self.normal_name}.called.cr.seg")
            
            normal_command = (
                f"gatk CallCopyRatioSegments --input {self.cr_seg_normal} "
                f"--output {self.called_cr_seg_normal}"
            )
            
            
            if self.neutral_segment_copy_ratio_lower_bound:
                normal_command += f" --neutral-segment-copy-ratio-lower-bound {self.neutral_segment_copy_ratio_lower_bound}"
            if self.neutral_segment_copy_ratio_upper_bound:
                normal_command += f" --neutral-segment-copy-ratio-upper-bound {self.neutral_segment_copy_ratio_upper_bound}"
            if self.outlier_neutral_segment_copy_ratio_z_score_threshold:
                normal_command += f" --outlier-neutral-segment-copy-ratio-z-score-threshold {self.outlier_neutral_segment_copy_ratio_z_score_threshold}"
            if self.calling_copy_ratio_z_score_threshold:
                normal_command += f" --calling-copy-ratio-z-score-threshold {self.calling_copy_ratio_z_score_threshold}"        
            print(normal_command)
            try:
                result = subprocess.run(normal_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"CallCopyRatioSegments- An error occurred while calling copy ratio segments of normal sample: {stderr_output}")
                else:
                    logging.info(f"CallCopyRatioSegments - Copy ratio segments called successfully (normal): {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somatic_pip_cnv.py - An internal error occurred while calling copy ratio segments of normal sample: {e}")
                
                
    def plot_dcr(self):
        
        logging.info("Plotting denoised copy ratios - Tumor")
        
        self.tumor_plots_folder = os.path.join(self.tumor_output_folder, "plots")
        
        os.makedirs(self.tumor_plots_folder, exist_ok=True)
                
        tumor_command = (
            f"gatk PlotDenoisedCopyRatios --standardized-copy-ratios {self.tumor_standardized} "
            f"--denoised-copy-ratios {self.tumor_denoised} "
            f"--sequence-dictionary {self.ref_dict} "
            f"--output {self.tumor_plots_folder} "
            f"--output-prefix {self.tumor_name}"
        )
        
        if self.minimum_contig_length:
            tumor_command += f" --minimum-contig-length {self.minimum_contig_length}"
        if self.maximum_copy_ratio:
            tumor_command += f" --maximum-copy-ratio {self.maximum_copy_ratio}"
        if self.point_size_copy_ratio:
            tumor_command += f" --point-size-copy-ratio {self.point_size_copy_ratio}"
        print(tumor_command)
        try:
            result = subprocess.run(tumor_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"PlotDenoisedCopyRatios- An error occurred while plotting denoised copy ratios of tumor sample: {stderr_output}")
            else:
                logging.info(f"PlotDenoisedCopyRatios - Denoised copy ratio plots created successfully (tumor): {stdout_output}")
                
        except Exception as e:
            logging.error(f"somatic_pip_cnv.py - An internal error occurred while plotting denoised copy ratios of tumor sample: {e}")              
        
        if self.normal_bam:
            
            self.normal_plots_folder = os.path.join(self.normal_output_folder, "plots")
            
            os.makedirs(self.normal_plots_folder, exist_ok=True)            
            
            logging.info("Plotting denoised copy ratios - Normal")
                    
            normal_command = (
                f"gatk PlotDenoisedCopyRatios --standardized-copy-ratios {self.normal_standardized} "
                f"--denoised-copy-ratios {self.normal_denoised} "
                f"--sequence-dictionary {self.ref_dict} "
                f"--output {self.normal_plots_folder} "
                f"--output-prefix {self.normal_name}"
            )
            
            if self.minimum_contig_length:
                normal_command += f" --minimum-contig-length {self.minimum_contig_length}"
            if self.maximum_copy_ratio:
                normal_command += f" --maximum-copy-ratio {self.maximum_copy_ratio}"
            if self.point_size_copy_ratio:
                normal_command += f" --point-size-copy-ratio {self.point_size_copy_ratio}"
            print(normal_command)
            try:
                result = subprocess.run(normal_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"PlotDenoisedCopyRatios- An error occurred while plotting denoised copy ratios of normal sample: {stderr_output}")
                else:
                    logging.info(f"PlotDenoisedCopyRatios - Denoised copy ratio plots created successfully (normal): {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somatic_pip_cnv.py - An internal error occurred while plotting denoised copy ratios of normal sample: {e}")
                
                
    def plot_msg(self):
        
        logging.info("Plotting modeled segments - Tumor")

        self.tumor_segments = os.path.join(self.tsg_folder, f"{self.tumor_name}.modelFinal.seg")
        self.thets =  os.path.join(self.tsg_folder, f"{self.tumor_name}.hets.tsv")
        
        tumor_command = (
            f"gatk PlotModeledSegments --denoised-copy-ratios {self.tumor_denoised} "
            f"--allelic-counts {self.thets} "
            f"--segments {self.tumor_segments} "
            f"--sequence-dictionary {self.ref_dict} "
            f"--output {self.tumor_plots_folder} "
            f"--output-prefix {self.tumor_name}"
        )
        
        if self.minimum_contig_length:
            tumor_command += f" --minimum-contig-length {self.minimum_contig_length}"
        if self.maximum_copy_ratio:
            tumor_command += f" --maximum-copy-ratio {self.maximum_copy_ratio}"
        if self.point_size_copy_ratio:
            tumor_command += f" --point-size-copy-ratio {self.point_size_copy_ratio}"
        print(tumor_command)
        try:
            result = subprocess.run(tumor_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"PlotModeledSegments- An error occurred while plotting modeled segments of tumor sample: {stderr_output}")
            else:
                logging.info(f"PlotModeledSegments - Modeled segments plots created successfully (tumor): {stdout_output}")
                
        except Exception as e:
            logging.error(f"somatic_pip_cnv.py - An internal error occurred while plotting modeled segments of tumor sample: {e}")              
        
        if self.normal_bam:
            
            logging.info("Plotting modeled segments - Normal")

            self.normal_segments = os.path.join(self.nsg_folder, f"{self.normal_name}.modelFinal.seg")
            self.nhets =  os.path.join(self.nsg_folder, f"{self.normal_name}.hets.tsv")
            
            normal_command = (
                f"gatk PlotModeledSegments --denoised-copy-ratios {self.normal_denoised} "
                f"--allelic-counts {self.nhets} "
                f"--segments {self.normal_segments} "
                f"--sequence-dictionary {self.ref_dict} "
                f"--output {self.normal_plots_folder} "
                f"--output-prefix {self.normal_name}"
            )
            
            if self.minimum_contig_length:
                normal_command += f" --minimum-contig-length {self.minimum_contig_length}"
            if self.maximum_copy_ratio:
                normal_command += f" --maximum-copy-ratio {self.maximum_copy_ratio}"
            if self.point_size_copy_ratio:
                normal_command += f" --point-size-copy-ratio {self.point_size_copy_ratio}"
            print(normal_command)
            try:
                result = subprocess.run(normal_command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                stderr_output = result.stderr.decode('utf-8').strip()
                stdout_output = result.stdout.decode('utf-8').strip()

                if result.returncode != 0:
                    logging.error(f"PlotModeledSegments- An error occurred while plotting modeled segments of normal sample: {stderr_output}")
                else:
                    logging.info(f"PlotModeledSegments - Modeled segments plots created successfully (normal): {stdout_output}")
                    
            except Exception as e:
                logging.error(f"somatic_pip_cnv.py - An internal error occurred while plotting modeled segments of normal sample: {e}")
                
                
def process_samples(config_path):
    
    with open(config_path, 'r') as f:
        config = json.load(f)

    algotype = config.get("output_folder")
    index_fasta = config.get("index_fasta")
    threads = config.get("aligner_threads")
    remove_all_duplicates = config.get("remove_all_duplicates")
    remove_sequencing_duplicates = config.get("remove_sequencing_duplicates")
    use_mark_duplicates = config.get("use_gatk_mark_duplicates")
    use_dragen = config.get("use_dragen")
    run_name = config.get("run_name")
    output_path = config.get("output_folder")
    ref_fasta = config.get("ref_fasta")
    ref_dict = config.get("ref_dict")
    interval_list = config.get("interval_list")
    tumor_samples = config.get("tumor_samples")
    normal_samples = config.get("normal_samples")
    common_sites = config.get("common_sites")
    pon = config.get("pon")
    blacklist_intervals = config.get("blacklist_intervals")
    minimum_base_quality = config.get("minimum_base_quality")
    number_of_eigensamples = config.get("number_of_eigensamples")
    minimum_total_allele_count_case = config.get("minimum_total_allele_count_case")
    minimum_total_allele_count_normal = config.get("minimum_total_allele_count_normal")
    genotyping_homozygous_log_ratio_threshold = config.get("genotyping_homozygous_log_ratio_threshold")
    genotyping_base_error_rate = config.get("genotyping_base_error_rate")
    maximum_number_of_segments_per_chromosome = config.get("maximum_number_of_segments_per_chromosome")
    kernel_variance_copy_ratio = config.get("kernel_variance_copy_ratio")
    kernel_variance_allele_fraction = config.get("kernel_variance_allele_fraction")
    kernel_scaling_allele_fraction = config.get("kernel_scaling_allele_fraction")
    kernel_approximation_dimension = config.get("kernel_approximation_dimension")
    window_size = config.get("window_size")
    number_of_changepoints_penalty_factor = config.get("number_of_changepoints_penalty_factor")
    minor_allele_fraction_prior_alpha = config.get("minor_allele_fraction_prior_alpha")
    number_of_samples_copy_ratio = config.get("number_of_samples_copy_ratio")
    number_of_burn_in_samples_copy_ratio = config.get("number_of_burn_in_samples_copy_ratio")
    number_of_samples_allele_fraction = config.get("number_of_samples_allele_fraction")
    number_of_burn_in_samples_allele_fraction = config.get("number_of_burn_in_samples_allele_fraction")
    smoothing_credible_interval_threshold_copy_ratio = config.get("smoothing_credible_interval_threshold_copy_ratio")
    smoothing_credible_interval_threshold_allele_fraction = config.get("smoothing_credible_interval_threshold_allele_fraction")
    maximum_number_of_smoothing_iterations = config.get("maximum_number_of_smoothing_iterations")
    number_of_smoothing_iterations_per_fit = config.get("number_of_smoothing_iterations_per_fit")
    neutral_segment_copy_ratio_lower_bound = config.get("neutral_segment_copy_ratio_lower_bound")
    neutral_segment_copy_ratio_upper_bound = config.get("neutral_segment_copy_ratio_upper_bound")
    outlier_neutral_segment_copy_ratio_z_score_threshold = config.get("outlier_neutral_segment_copy_ratio_z_score_threshold")
    calling_copy_ratio_z_score_threshold = config.get("calling_copy_ratio_z_score_threshold")
    minimum_contig_length = config.get("minimum_contig_length")
    maximum_copy_ratio = config.get("maximum_copy_ratio")
    point_size_copy_ratio = config.get("point_size_copy_ratio")
    point_size_allele_fraction = config.get("point_size_allele_fraction")
    
    if index_fasta:
    
        indexer = Indexer(reference_fasta=ref_fasta, algotype=algotype)
        indexer.index()

    for tumor_sample in tumor_samples: 
        
        tumor_sample_id = tumor_sample.get("sample_id")
        tumor_r1 = tumor_sample.get("r1")
        tumor_r2 = tumor_sample.get("r2")
        
        run_path = os.path.join(output_path, run_name)
        tumor_path = os.path.join(output_path, run_name, tumor_sample_id)
        
        os.makedirs(tumor_path, exist_ok=True)

        mapper = Mapper(reference_fasta=ref_fasta, sample_id=tumor_sample_id, r1=tumor_r1, r2=tumor_r2, output_path=tumor_path, run_name=run_name, threads=threads, use_dragen=use_dragen)
        sam_file = mapper.map()

        deduplicator = deduplicate(sam_path=sam_file, sample_id=tumor_sample_id, threads=threads, remove_all_duplicates=remove_all_duplicates, remove_sequencing_duplicates=remove_sequencing_duplicates, use_mark_duplicates=use_mark_duplicates)
        
        tumor_bam = deduplicator.convert_markdedup()
        
        index = next((i for i, sample in enumerate(tumor_samples) if sample["sample_id"] == tumor_sample_id))
        
        normal_bam = None

        if normal_samples:
            
            normal_sample = normal_samples[index]
            
            normal_sample_id = normal_sample.get("sample_id")
            normal_sample_r1 = normal_sample.get("r1")
            normal_sample_r2 = normal_sample.get("r2")
            
            normal_path = os.path.join(output_path, run_name, normal_sample_id)
            
                      
            os.makedirs(normal_path, exist_ok=True)
                        
            mapper = Mapper(reference_fasta=ref_fasta, sample_id=normal_sample_id, r1=normal_sample_r1, r2=normal_sample_r2, output_path=normal_path, run_name=run_name, threads=threads, use_dragen=use_dragen)
            sam_file = mapper.map()

            deduplicator = deduplicate(sam_path=sam_file, sample_id=normal_sample_id, threads=threads, remove_all_duplicates=remove_all_duplicates, remove_sequencing_duplicates=remove_sequencing_duplicates, use_mark_duplicates=use_mark_duplicates)
            normal_bam = deduplicator.convert_markdedup()
            
            
        pipeline = SomaticPipelineCNV(
            run_path,
            tumor_path,
            tumor_sample_id, 
            ref_fasta,
            ref_dict,
            interval_list,
            common_sites,
            tumor_bam,
            pon,
            normal_path,
            normal_sample_id,
            normal_bam,
            blacklist_intervals,
            minimum_base_quality,
            number_of_eigensamples,
            minimum_total_allele_count_case,
            minimum_total_allele_count_normal,
            genotyping_homozygous_log_ratio_threshold,
            genotyping_base_error_rate,
            maximum_number_of_segments_per_chromosome,
            kernel_variance_copy_ratio,
            kernel_variance_allele_fraction,
            kernel_scaling_allele_fraction,
            kernel_approximation_dimension,
            window_size,
            number_of_changepoints_penalty_factor,
            minor_allele_fraction_prior_alpha,
            number_of_samples_copy_ratio,
            number_of_burn_in_samples_copy_ratio,
            number_of_samples_allele_fraction,
            number_of_burn_in_samples_allele_fraction,
            smoothing_credible_interval_threshold_copy_ratio,
            smoothing_credible_interval_threshold_allele_fraction,
            maximum_number_of_smoothing_iterations,
            number_of_smoothing_iterations_per_fit,
            neutral_segment_copy_ratio_lower_bound,
            neutral_segment_copy_ratio_upper_bound,
            outlier_neutral_segment_copy_ratio_z_score_threshold,
            calling_copy_ratio_z_score_threshold,
            minimum_contig_length,
            maximum_copy_ratio,
            point_size_copy_ratio,
            point_size_allele_fraction
        )
        
        intervals_processed = False
        
        if not intervals_processed:
            pipeline.preprocess_intervals()
            intervals_processed = True
            
        pipeline.collect_counts()
        pipeline.collect_allelic_counts()
        pipeline.denoise_counts()
        pipeline.model_segments()
        pipeline.call_crs()
        pipeline.plot_dcr()
        pipeline.plot_msg()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run somatic CNV pipeline with configuration.')
    parser.add_argument('--config', type=str, help='Path to the JSON configuration file.')

    args = parser.parse_args()

    process_samples(args.config)

