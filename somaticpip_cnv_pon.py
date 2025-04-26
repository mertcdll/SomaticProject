import os
import logging
import subprocess
import json
import argparse
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class PanelOfNormals:
    def __init__(
            self,
            output_folder: str,
            ref_fasta: str,
            ref_dict: str,
            interval_list: str,
            pon_id: str,
            normal_bams: List[str],
            gc_correction: bool = None,
            bin_length: int = None,
            blacklist_intervals: str = None,
            padding: int = None,
            do_impute_zeros: bool = None,
            number_of_eigensamples: int = None,
            feature_query_lookahead: int = None,
            minimum_interval_median_percentile: float = None,
            maximum_zeros_in_sample_percentage: float = None,
            maximum_zeros_in_interval_percentage: float = None,
            extreme_sample_median_percentile: float = None,
            extreme_outlier_truncation_percentile: float = None,
            maximum_chunk_size: int = None):

        self.output_folder = output_folder
        self.ref_fasta = ref_fasta
        self.ref_dict = ref_dict
        self.interval_list = interval_list
        self.pon_id = pon_id
        self.normal_bams = normal_bams
        self.gc_correction = gc_correction
        self.bin_length = bin_length
        self.blacklist_intervals = blacklist_intervals
        self.padding = padding
        self.do_impute_zeros = do_impute_zeros
        self.number_of_eigensamples = number_of_eigensamples
        self.feature_query_lookahead = feature_query_lookahead
        self.minimum_interval_median_percentile = minimum_interval_median_percentile
        self.maximum_zeros_in_sample_percentage = maximum_zeros_in_sample_percentage
        self.maximum_zeros_in_interval_percentage = maximum_zeros_in_interval_percentage
        self.extreme_sample_median_percentile = extreme_sample_median_percentile
        self.extreme_outlier_truncation_percentile = extreme_outlier_truncation_percentile
        self.maximum_chunk_size = maximum_chunk_size
        
        self.intervals_name = None
        self.preprocessed_intervals = None
        self.annotated_intervals = None
        self.counts_folder = None
        self.count_commands = []

    def preprocess_intervals(self):

        logging.info("Processing interval list")
        
        self.intervals_name = os.path.basename(self.interval_list).split("interval_list")[0]
        preprocessed_name = self.intervals_name + "preprocessed.interval_list"
        self.preprocessed_intervals = os.path.join(self.output_folder, preprocessed_name)

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


    def annotate_intervals(self):

        logging.info("Annotating invervals")
        
        annotated_name = self.intervals_name + "preprocessed.annotated.tsv"
        self.annotated_intervals = os.path.join(self.output_folder, annotated_name)
        
        command = (
            f"gatk AnnotateIntervals -L {self.preprocessed_intervals} "
            f"--reference {self.ref_fasta} "
            f"--interval-merging-rule OVERLAPPING_ONLY "
            f"--output {self.annotated_intervals}"
        )

        if self.feature_query_lookahead:
            command += f" --feature-query-lookahead {self.feature_query_lookahead}"


        try:
            result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"AnnotateIntervals - An error occurred while annotating the intervals: {stderr_output}")
            else:
                logging.info(f"AnnotateIntervals - Intervals have been annotated successfully: {stdout_output}")
                   
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while annotating the intervals: {e}")

    
    def collect_counts(self):

        logging.info("Collecting Read Counts")

        self.counts_folder = os.path.join(self.output_folder, "read_counts")

        os.makedirs(self.counts_folder, exist_ok=True)
        
        for bam in self.normal_bams:
            bam_name = os.path.basename(bam)
            sample_name = os.path.splitext(bam_name)[0]
            counts_name = sample_name + ".counts.hdf5"
            counts_path = os.path.join(self.counts_folder, counts_name)


            count_command = (
                f"gatk CollectReadCounts -L {self.preprocessed_intervals} "
                f"--input {bam} --reference {self.ref_fasta} "
                f"--format HDF5 --interval-merging-rule OVERLAPPING_ONLY "
                f"--output {counts_path}"
            )

            self.count_commands.append(count_command)

        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(self.run_collectreadcounts_commands, command) for command in self.count_commands]
            for future in futures:
                try:
                    future.result() 
                except Exception as e:
                    logging.error(f"Error occurred during running CollectReadCounts commands: {e}")

    def run_collectreadcounts_commands(self, command):
        
        result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stderr_output = result.stderr.decode('utf-8').strip()
        stdout_output = result.stdout.decode('utf-8').strip()

        if result.returncode != 0:
            logging.error(f"CollectReadCounts - An error occurred while running command: {stderr_output}")
        else:
            logging.info(f"CollectReadCounts - Command completed successfully: {stdout_output}")

    
    def create_pon(self):

        logging.info("Creating panel of normals")

        count_files = os.listdir(self.counts_folder)

        count_paths = [os.path.join(self.counts_folder, count_file) for count_file in count_files]

        pon_name = self.pon_id + ".pon.hdf5"

        self.pon_path = os.path.join(self.output_folder, pon_name)

        line = ""

        for item in count_paths:
            line += "--input " + item + " "

        command = (
            f"gatk CreateReadCountPanelOfNormals {line.strip()} "
            f"--output {self.pon_path}"
        )

        if self.minimum_interval_median_percentile:
            command += f" --minimum-interval-median-percentile {self.minimum_interval_median_percentile}"
        if self.maximum_zeros_in_sample_percentage:
            command += f" --maximum-zeros-in-sample-percentage {self.maximum_zeros_in_sample_percentage}"
        if self.maximum_zeros_in_interval_percentage:
            command += f" --maximum-zeros-in-interval-percentage {self.maximum_zeros_in_interval_percentage}"
        if self.do_impute_zeros:
            command += f" --do-impute-zeros {self.do_impute_zeros}"
        if self.extreme_outlier_truncation_percentile:
            command +=  f" --extreme-outlier-truncation-percentile {self.extreme_outlier_truncation_percentile}"
        if self.extreme_sample_median_percentile:
            command += f" --number-of-eigensamples {self.number_of_eigensamples}"
        if self.maximum_chunk_size:
            command +=  f" --maximum-chunk-size {self.maximum_chunk_size}"
        if self.gc_correction: 
            command += f" --annotated-intervals {self.annotated_intervals}"

        try:
            result = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_output = result.stderr.decode('utf-8').strip()
            stdout_output = result.stdout.decode('utf-8').strip()

            if result.returncode != 0:
                logging.error(f"CreateReadCountPanelOfNormals - An error occurred while creating the panel of normals: {stderr_output}")
            else:
                logging.info(f"CreateReadCountPanelOfNormals - Panel of normals have been created successfully: {stdout_output}")
                   
        except Exception as e:
            logging.error(f"somaticpip.py - An internal error occurred while creating the panel of normals: {e}")

        
        return self.pon_path
    
    def run_workflow(self):
        
        self.preprocess_intervals()
        if self.gc_correction:
            self.annotate_intervals()
        self.collect_counts()
        pon_file = self.create_pon()
        return pon_file

class Config:

    def __init__(self, config_file):
        with open(config_file, 'r') as f:
            config = json.load(f)

        
        self.output_folder = config.get("output_folder")
        self.ref_fasta = config.get("ref_fasta")
        self.ref_dict = config.get("ref_dict")
        self.interval_list = config.get("interval_list")
        self.pon_id = config.get("pon_id")
        self.normal_bams = config.get("normal_bams")
        self.gc_correction = config.get("gc_correction")
        self.bin_length = config.get("bin_length")
        self.blacklist_intervals = config.get("blacklist_intervals")
        self.padding = config.get("padding")
        self.do_impute_zeros = config.get("do_impute_zeros")
        self.number_of_eigensamples = config.get("number_of_eigensamples")
        self.feature_query_lookahead = config.get("feature_query_lookahead")
        self.minimum_interval_median_percentile = config.get("minimum_interval_median_percentile")
        self.maximum_zeros_in_sample_percentage = config.get("maximum_zeros_in_sample_percentage")
        self.maximum_zeros_in_interval_percentage = config.get("maximum_zeros_in_interval_percentage")
        self.extreme_sample_median_percentile = config.get("extreme_sample_median_percentile")
        self.extreme_outlier_truncation_percentile = config.get("extreme_outlier_truncation_percentile")
        self.maximum_chunk_size = config.get("maximum_chunk_size")


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Run panel of normals with configuration.')
    parser.add_argument('config_file', type=str, help='Path to the JSON configuration file.')

    args = parser.parse_args()

    config = Config(args.config_file)

    pon = PanelOfNormals(
        config.output_folder,
        config.ref_fasta,
        config.ref_dict,
        config.interval_list,
        config.pon_id,
        config.normal_bams,
        config.gc_correction,
        config.bin_length,
        config.blacklist_intervals,
        config.padding,
        config.do_impute_zeros,
        config.number_of_eigensamples,
        config.feature_query_lookahead,
        config.minimum_interval_median_percentile,
        config.maximum_zeros_in_sample_percentage,
        config.maximum_zeros_in_interval_percentage,
        config.extreme_sample_median_percentile,
        config.extreme_outlier_truncation_percentile,
        config.maximum_chunk_size
    )

    pon.run_workflow()