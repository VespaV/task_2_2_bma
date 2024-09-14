import pandas as pd
import subprocess


class FindCompleteHomologous:
    def __init__(self, bed_file, genome, min_length = 0, target_fasta='target_sequences.fasta', blast_result_file='blast_result.txt', output_file='results_with_coordinates.tsv'):
        self.bed_file = bed_file
        self.genome = genome
        self.min_length = min_length
        self.target_fasta = target_fasta
        self.blast_result_file = blast_result_file
        self.output_file = output_file

    def extract_complete_homologues(self):
        try:
            self.get_fasta_from_bed()
            print(f'Сформирован файл Fasta из файла BED {self.target_fasta}')
        except Exception as e:
            print(e)
        try:
            print(f'Производится выравнивание BLAST последовательностей {self.target_fasta} с геномом {self.genome}')
            self.blast_amplicons()
            print(f'Выравнивание завешено. Результаты записаны в {self.blast_result_file}')
        except Exception as e:
            print(e)
        try:
            print('Производится поиск нецелевых регионов')
            self.blast_amplicons()
            blast_results = self.create_df_from_blast_results()
            filtered_not_target = self.check_not_target(blast_results)
            filter_length_df = self.filter_length(filtered_not_target)
            final_df = self.join_columns(filter_length_df)
            final_df.to_csv(self.output_file, index=False, sep='\t', header=True)
            print(f'Файл с координатами целевых и нецелевых регионов сформирован: {self.output_file}')
        except Exception as e:
            print(e)

    def get_fasta_from_bed(self):
        cmd = f"bedtools getfasta -fi {self.genome} -bed {self.bed_file} > {self.target_fasta}"
        subprocess.run(cmd, shell=True, check=True)

    def blast_amplicons(self):
        cmd = f"blastn -query {self.target_fasta} -db {self.genome } -out {self.blast_result_file} -perc_identity 100 -outfmt 6"
        subprocess.run(cmd, shell=True, check=True)

    def create_df_from_blast_results(self):
        columns = [
            "amplicon", "chr", "% ident", "length", "mismatch", "gapopen",
            "align_start", "align_end", "genome_start", "genome_end", "evalue", "bitscore"
        ]
        blast_results = pd.read_csv(self.blast_result_file, sep="\t", names=columns)

        return blast_results

    @staticmethod
    def extract_coordinates(amplicon):
        chrom, coords = amplicon.split(":")
        start, end = map(int, coords.split("-"))
        return chrom, start, end

    def check_not_target(self, blast_results):
        filtered_results = []
        for _, row in blast_results.iterrows():
            chrom, target_start, target_end = self.extract_coordinates(row['amplicon'])

            if row['chr'] != chrom and not (row['genome_start'] >= target_start and row['genome_end'] <= target_end):
                filtered_results.append(row)

        return pd.DataFrame(filtered_results)

    def filter_length(self, blast_results):
        return blast_results[(abs(blast_results['genome_start'] - blast_results['genome_end'] + 1)) > self.min_length]

    @staticmethod
    def join_columns(df):
        df = df.copy()
        df['100% homolog'] = df.apply(lambda row: f"{row['chr']}:{row['genome_start']}-{row['genome_end']}", axis=1)
        df = df[['amplicon', '100% homolog']].reset_index(drop=True)
        return df


if __name__ == "__main__":
    prossesor = FindCompleteHomologous(bed_file='IAD143293_241_Designed.bed', genome='hg19.fa', min_length=50)
    prossesor.extract_complete_homologues()
