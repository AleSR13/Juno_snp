import imp
import pandas as pd
import pathlib
import argparse
import sys
import zipfile
import subprocess

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GenomeApi as DatasetsGenomeApi

from ncbi.datasets.package import dataset


def join_full_dfs(left, right):
    merged = pd.merge(left, right, on='#ID', how='outer')
    return merged


def find_referenceseeker_res(input_dir):
    input_dir = pathlib.Path(input_dir)
    return list(input_dir.glob('**/referenceseeker_*.tab'))


def append_all_referenceseeker_res(file_list):
    res_referenceseeker = []
    for file_ in file_list:
        sample_name = file_.stem.replace('referenceseeker_', '')
        rows_to_skip = subprocess.check_output(
            f'grep -n "#ID" {file_} | cut -d : -f 1', shell=True
        )
        rows_to_skip = int(rows_to_skip) - 1
        sample_res = pd.read_csv(
            file_, skiprows=rows_to_skip, sep='\t', index_col=None, header=0
        )
        sample_res = sample_res[['#ID', 'ANI', 'Mash Distance', 'Con. DNA']]
        sample_res['Sample'] = sample_name
        res_referenceseeker.append(sample_res)
        res = pd.concat(res_referenceseeker)
    return res


def score_candidates(file_list):
    df_all = append_all_referenceseeker_res(file_list)
    candidates = df_all.groupby('#ID')[['ANI', 'Mash Distance', 'Con. DNA']]\
        .agg(['mean', 'size'])
    candidates.columns = ['_'.join(head) for head in candidates.columns]
    imp_cols = ['ANI_size', 'ANI_mean', 'Con. DNA_mean', 'Mash Distance_mean']
    candidates = candidates.sort_values(
            by=imp_cols, 
            ascending=[False, False, False, True]
        )
    return candidates[imp_cols]


def get_best_hit(candidates_df):
    best_refseq_hit = candidates_df.index[0]
    print(f'The best scored reference genome for this dataset is: {best_refseq_hit}')
    return best_refseq_hit


def download_from_ncbi(accession_nr, output_dir):
    output_dir = pathlib.Path(output_dir)
    accessions = [accession_nr]
    zipfile_name = output_dir.joinpath('ncbi_dataset.zip')
    ref_genome = output_dir.joinpath('ref_genome.fasta')
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        try:
            print("Begin download of genome data package ...")
            genome_ds_download = genome_api.download_assembly_package(
                accessions,
                include_annotation_type=["RNA_FASTA", "PROT_FASTA"],
                _preload_content=False,
            )
            with open(zipfile_name, "wb") as f:
                f.write(genome_ds_download.data)
            print(f"Download completed -- see {zipfile_name}")
        except DatasetsApiException as e:
            sys.exit(
                f"Exception when calling download_assembly_package: {e}\n"
            )
    with zipfile.ZipFile(zipfile_name, 'r') as zip_ref:
        zip_ref.extractall(output_dir)
    for item in output_dir.joinpath('ncbi_dataset', 'data').glob('**/*.fna'):
        item.rename(ref_genome)
    output_dir.joinpath('README.md').unlink()
    zipfile_name.unlink()
    return ref_genome


def get_user_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-o', '--output',
        metavar='DIR',
        default='output',
        type=pathlib.Path
    )
    parser.add_argument(
        '-id', '--input-dir',
        metavar='DIR',
        required=not ('-if' in sys.argv or '--input-files' in sys.argv),
        type=pathlib.Path,
        default=None
    )
    parser.add_argument(
        '-if', '--input-files',
        metavar='FILES',
        nargs='+',
        required=not ('-id' in sys.argv or '--input-dir' in sys.argv),
        type=pathlib.Path,
        default=None
    )
    args = parser.parse_args()
    if args.input_files is None:
        args.input_files = find_referenceseeker_res(args.input_dir)
    return args


def main():
    args = get_user_args()
    res = score_candidates(args.input_files)
    res.to_csv(args.output.joinpath('scores_refseq_candidates.csv'))
    best_hit = get_best_hit(res)
    download_from_ncbi(best_hit, output_dir=args.output)


if __name__ == '__main__':
    main()