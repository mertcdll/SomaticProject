import os
import logging
import subprocess
import json
import argparse
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class SomaticPipelineCNV:
    def __init__(
            self,
            output_folder: str,
            ref_fasta: str,
            ref_dict: str,
            processed_interval_list: str,
            common_sites: str,
            tumor_bam: str,
            pon: str,
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
            point_size_allele_fraction: float = None):
        
        self.output_folder = output_folder
        self.ref_fasta = ref_fasta
        self.ref_dict = ref_dict
        self.processed_interval_list = processed_interval_list
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


        self.tumor_name = None
        self.counts_folder = None
        self.tumor_counts = None
        self.normal_name = None
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
        
    def collect_counts(self):

        logging.info("Collecting tumor read counts")

        self.tumor_name = os.path.splitext(os.path.basename(self.tumor_bam))[0]  
        
        self.counts_folder = os.path.join(self.output_folder, "counts")
        
        os.makedirs(self.counts_folder, exist_ok=True)

        self.tumor_counts = os.path.join(self.counts_folder, f"{self.tumor_name}.counts.hdf5")

        tumor_command = (
            f"gatk CollectReadCounts -L {self.processed_interval_list} "
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

            self.normal_name = os.path.splitext(os.path.basename(self.normal_bam))[0]

            self.normal_counts = os.path.join(self.counts_folder, f"{self.normal_name}.counts.hdf5")


            normal_command = (
                f"gatk CollectReadCounts -L {self.processed_interval_list} "
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
        
        self.tumor_alcounts = os.path.join(self.counts_folder, f"{self.tumor_name}.allelicCounts.tsv")

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

            self.normal_alcounts = os.path.join(self.counts_folder, f"{self.normal_name}.allelicCounts.tsv")

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
        
        self.tumor_standardized = os.path.join(self.counts_folder, f"{self.tumor_name}.standardizedCR.tsv")
        self.tumor_denoised = os.path.join(self.counts_folder, f"{self.tumor_name}.denoisedCR.tsv")

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
            
            self.normal_standardized = os.path.join(self.counts_folder, f"{self.normal_name}.standardizedCR.tsv")
            self.normal_denoised = os.path.join(self.counts_folder, f"{self.normal_name}.denoisedCR.tsv")

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
        
        self.tsg_folder = os.path.join(self.output_folder, "tumor_segments")
        
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
            
            self.nsg_folder = os.path.join(self.output_folder, "normal_segments")
            
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
        
        self.plots_folder = os.path.join(self.output_folder, "plots")
        
        os.makedirs(self.plots_folder, exist_ok=True)
                
        tumor_command = (
            f"gatk PlotDenoisedCopyRatios --standardized-copy-ratios {self.tumor_standardized} "
            f"--denoised-copy-ratios {self.tumor_denoised} "
            f"--sequence-dictionary {self.ref_dict} "
            f"--output {self.plots_folder} "
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
            
            logging.info("Plotting denoised copy ratios - Normal")
                    
            normal_command = (
                f"gatk PlotDenoisedCopyRatios --standardized-copy-ratios {self.normal_standardized} "
                f"--denoised-copy-ratios {self.normal_denoised} "
                f"--sequence-dictionary {self.ref_dict} "
                f"--output {self.plots_folder} "
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
            f"--output {self.plots_folder} "
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
                f"--output {self.plots_folder} "
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
                
                
    def run_workflow(self):
        self.collect_counts()
        self.collect_allelic_counts()
        self.denoise_counts()
        self.model_segments()
        self.call_crs()
        self.plot_dcr()
        self.plot_msg()
        
    
class Config:

    def __init__(self, config_file):
        with open(config_file, 'r') as f:
            config = json.load(f)

        self.output_folder = config.get("output_folder")
        self.ref_fasta = config.get("ref_fasta")
        self.ref_dict = config.get("ref_dict")
        self.processed_interval_list = config.get("processed_interval_list")
        self.common_sites = config.get("common_sites")
        self.tumor_bam = config.get("tumor_bam")
        self.pon = config.get("pon")
        self.normal_bam = config.get("normal_bam")
        self.blacklist_intervals = config.get("blacklist_intervals")
        self.minimum_base_quality = config.get("minimum_base_quality")
        self.number_of_eigensamples = config.get("number_of_eigensamples")
        self.minimum_total_allele_count_case = config.get("minimum_total_allele_count_case")
        self.minimum_total_allele_count_normal = config.get("minimum_total_allele_count_normal")
        self.genotyping_homozygous_log_ratio_threshold = config.get("genotyping_homozygous_log_ratio_threshold")
        self.genotyping_base_error_rate = config.get("genotyping_base_error_rate")
        self.maximum_number_of_segments_per_chromosome = config.get("maximum_number_of_segments_per_chromosome")
        self.kernel_variance_copy_ratio = config.get("kernel_variance_copy_ratio")
        self.kernel_variance_allele_fraction = config.get("kernel_variance_allele_fraction")
        self.kernel_scaling_allele_fraction = config.get("kernel_scaling_allele_fraction")
        self.kernel_approximation_dimension = config.get("kernel_approximation_dimension")
        self.window_size = config.get("window_size")
        self.number_of_changepoints_penalty_factor = config.get("number_of_changepoints_penalty_factor")
        self.minor_allele_fraction_prior_alpha = config.get("minor_allele_fraction_prior_alpha")
        self.number_of_samples_copy_ratio = config.get("number_of_samples_copy_ratio")
        self.number_of_burn_in_samples_copy_ratio = config.get("number_of_burn_in_samples_copy_ratio")
        self.number_of_samples_allele_fraction = config.get("number_of_samples_allele_fraction")
        self.number_of_burn_in_samples_allele_fraction = config.get("number_of_burn_in_samples_allele_fraction")
        self.smoothing_credible_interval_threshold_copy_ratio = config.get("smoothing_credible_interval_threshold_copy_ratio")
        self.smoothing_credible_interval_threshold_allele_fraction = config.get("smoothing_credible_interval_threshold_allele_fraction")
        self.maximum_number_of_smoothing_iterations = config.get("maximum_number_of_smoothing_iterations")
        self.number_of_smoothing_iterations_per_fit = config.get("number_of_smoothing_iterations_per_fit")
        self.neutral_segment_copy_ratio_lower_bound = config.get("neutral_segment_copy_ratio_lower_bound")
        self.neutral_segment_copy_ratio_upper_bound = config.get("neutral_segment_copy_ratio_upper_bound")
        self.outlier_neutral_segment_copy_ratio_z_score_threshold = config.get("outlier_neutral_segment_copy_ratio_z_score_threshold")
        self.calling_copy_ratio_z_score_threshold = config.get("calling_copy_ratio_z_score_threshold")
        self.minimum_contig_length = config.get("minimum_contig_length")
        self.maximum_copy_ratio = config.get("maximum_copy_ratio")
        self.point_size_copy_ratio = config.get("point_size_copy_ratio")
        self.point_size_allele_fraction = config.get("point_size_allele_fraction")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run somatic CNV pipeline with configuration.')
    parser.add_argument('config_file', type=str, help='Path to the JSON configuration file.')

    args = parser.parse_args()

    config = Config(args.config_file)

    pipeline = SomaticPipelineCNV(
        config.output_folder,
        config.ref_fasta,
        config.ref_dict,
        config.processed_interval_list,
        config.common_sites,
        config.tumor_bam,
        config.pon,
        config.normal_bam,
        config.blacklist_intervals,
        config.minimum_base_quality,
        config.number_of_eigensamples,
        config.minimum_total_allele_count_case,
        config.minimum_total_allele_count_normal,
        config.genotyping_homozygous_log_ratio_threshold,
        config.genotyping_base_error_rate,
        config.maximum_number_of_segments_per_chromosome,
        config.kernel_variance_copy_ratio,
        config.kernel_variance_allele_fraction,
        config.kernel_scaling_allele_fraction,
        config.kernel_approximation_dimension,
        config.window_size,
        config.number_of_changepoints_penalty_factor,
        config.minor_allele_fraction_prior_alpha,
        config.number_of_samples_copy_ratio,
        config.number_of_burn_in_samples_copy_ratio,
        config.number_of_samples_allele_fraction,
        config.number_of_burn_in_samples_allele_fraction,
        config.smoothing_credible_interval_threshold_copy_ratio,
        config.smoothing_credible_interval_threshold_allele_fraction,
        config.maximum_number_of_smoothing_iterations,
        config.number_of_smoothing_iterations_per_fit,
        config.neutral_segment_copy_ratio_lower_bound,
        config.neutral_segment_copy_ratio_upper_bound,
        config.outlier_neutral_segment_copy_ratio_z_score_threshold,
        config.calling_copy_ratio_z_score_threshold,
        config.minimum_contig_length,
        config.maximum_copy_ratio,
        config.point_size_copy_ratio,
        config.point_size_allele_fraction
    )

    pipeline.run_workflow()
