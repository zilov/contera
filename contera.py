#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 25.05.2021
# @author: Danil
# @contact: Zilov

import argparse
import os
from collections import defaultdict
from inspect import getsourcefile
#from datetime import datetime


def run_blast_nt(assembly_fasta, blast_db, threads, output_file, debug=False):
    command = f"blastn -query {assembly_fasta} -db {blast_db} -max_target_seqs 5 " \
              f"-outfmt '6 qseqid sseqid pident qstart qend length evalue sscinames staxids' " \
              f"-evalue 1e-5 -num_threads {threads} -out {output_file}"
    if os.path.exists(output_file):
        print("Already blasted, continue!\n")
        return True
    else:
        if debug:
            print(f"\nCommand to find assembly content\n\n{command}\n")
            return False
        else:
            print(f"\nRunning {command}\n")
            os.system(command)
            print("Blast of assembly content is done!\n")
            return True


def run_blast_adapters(assembly_fasta, adapters_db, threads, output_file, debug=False):
    command = f"blastn -query {assembly_fasta} -db {adapters_db} -outfmt 6 -max_target_seqs 1000 " \
              f"-task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -num_threads {threads} -dust yes " \
              f"-soft_masking true -evalue 0.00001 -searchsp 1750000000000 -out {output_file}"

    if os.path.exists(output_file):
        print("Adapters is already blasted, continue!\n")
        return True
    else:
        if debug:
            print(f"\nCommand to find adapters\n\n{command}\n")
            return False
        else:
            print(f"\nRunning {command}\n")
            os.system(command)
            print("Blast of adapters is done!\n")
            return True


def read_blast_matches(blast_out_file):
    contig2bactlength = defaultdict(list)
    with open(blast_out_file) as fh:
        for line in fh:
            line = line.strip().split("\t")
            contig = line[0]
            contig2bactlength[contig].append(line)
    return contig2bactlength


def fasta_reader(fasta_file):
    fasta = {}
    header = None
    with open(fasta_file) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if line.startswith(">"):
                if header:
                    fasta[header]="".join(seq)
                header = line.split()[0]
                seq = []
            else:
                seq.append(line)
    if header:
        fasta[header] = "".join(seq)
    return fasta


def fasta_length(fasta_directory):
    total_length = 0
    for seq in fasta_directory.values():
        total_length += len(seq)
    return total_length


def bacteria_percentage_in_contig(contig_to_bactlength):
    bact_percentage_in_contig = defaultdict(dict)
    for contig_blast in contig_to_bactlength:
        value = contig_to_bactlength[contig_blast]
        value.sort(key=lambda x: int(x[5]), reverse=True)
        bact_seq = defaultdict(list)
        i = 0
        for blast_match in value:
            start_stop = [int(blast_match[3]), int(blast_match[4])]
            bacteria = " ".join(blast_match[-2].split()[0:2])
            if not bact_seq[bacteria]:
                bact_seq[bacteria].append(start_stop)
                # print(f"written first {bacteria}")
            else:
                # print(f"whatchin inside of {bacteria}")
                # print(f"wanna put inside {start_stop}")
                i = 1
                for seq_in in bact_seq[bacteria]:
                    # print(f"will compare it to {seq_in}")
                    start_in = seq_in[0]
                    stop_in = seq_in[1]
                    start_out = start_stop[0]
                    stop_out = start_stop[1]
                    # print(bact_seq)
                    if start_out <= start_in and stop_out >= stop_in:
                        seq_in = start_stop
                    elif (start_in <= start_out < stop_in) and stop_out > stop_in:
                        seq_in[1] = stop_out
                    elif start_out < start_in and (start_in < stop_out <= stop_in):
                        seq_in[0] = start_out
                    elif start_out < start_in and stop_out < start_in:
                        if i == len(bact_seq[bacteria]):
                            bact_seq[bacteria].append(start_stop)
                        else:
                            i += 1
                            continue
                    elif start_out > stop_in and stop_out > stop_in:
                        if i == len(bact_seq[bacteria]):
                            bact_seq[bacteria].append(start_stop)
                        else:
                            i += 1
                            continue
                    elif start_out == start_in and stop_out == stop_in:
                        continue
                    elif start_out > start_in and stop_out < stop_in:
                        continue
        bact_percentage_in_contig[contig_blast] = bact_seq
    return bact_percentage_in_contig


def bacteria_sequence_in_contig(bact_percentage_in_contig):
    bact2seqparts = defaultdict(dict)
    for contig in bact_percentage_in_contig:
        bacteria2len = bact_percentage_in_contig[contig]
        bact_seq = defaultdict(list)
        for bacteria in bacteria2len:
            lengths = bacteria2len[bacteria]
            lengths.sort(key=lambda x: x[0], reverse=False)
            for blast_match in lengths:
                start_stop = [blast_match[0], blast_match[1]]
                if not bact_seq[bacteria]:
                    bact_seq[bacteria].append(start_stop)
                    # print(f"written first {bacteria} to {contig}")
                else:
                    # print(f"whatchin inside of {bacteria}")
                    # print(f"wanna put inside {start_stop}")
                    i = 1
                    for seq_in in bact_seq[bacteria]:
                        # print(f"will compare it to {seq_in}")
                        start_in = seq_in[0]
                        stop_in = seq_in[1]
                        start_out = start_stop[0]
                        stop_out = start_stop[1]
                        # print(bact_seq)
                        if start_out <= start_in and stop_out >= stop_in:
                            seq_in = start_stop
                        elif (start_in <= start_out < stop_in) and stop_out > stop_in:
                            seq_in[1] = stop_out
                        elif start_out < start_in and (start_in < stop_out <= stop_in):
                            seq_in[0] = start_out
                        elif start_out < start_in and stop_out < start_in:
                            if i == len(bact_seq[bacteria]):
                                bact_seq[bacteria].append(start_stop)
                            else:
                                i += 1
                                continue
                        elif start_out > stop_in and stop_out > stop_in:
                            if i == len(bact_seq[bacteria]):
                                bact_seq[bacteria].append(start_stop)
                            else:
                                i += 1
                                continue
                        elif start_out == start_in and stop_out == stop_in:
                            continue
                        elif start_out > start_in and stop_out < stop_in:
                            continue
        bact2seqparts[contig] = bact_seq
    return bact2seqparts


def contig_inside(bacteria_sequence_in_contig, scaffolds_directory, file_to_write_contig_content):
    with open(file_to_write_contig_content, "w") as fw:
        fw.write("# report is created by Contera tool, please cite us https://github.com/zilov/contera!\n")
        fw.write("# contig\tcontig_length\torganism\tstart\tstop\tpercentage_of_contig\n")
        contig_content = defaultdict(dict)
        for contig in bacteria_sequence_in_contig:
            bac_content = defaultdict(list)
            bacterias = bacteria_sequence_in_contig[contig]
            contig_length = len(scaffolds_directory[f">{contig}"])
            for bact in bacterias:
                bact_seqs = bacterias[bact]
                bact_length = 0
                for seq in bact_seqs:
                    start = seq[0]
                    stop = seq[1]
                    bact_length += stop - start
                    percentage = bact_length / contig_length * 100
                    #print(f"{contig}\t{bact}\t{start}\t{stop}")
                    fw.write(f"{contig}\t{contig_length}\t{bact}\t{start}\t{stop}\t{percentage}\n")
                match_percent = bact_length / contig_length * 100
                bac_content[bact] = [bact_length, match_percent]
            contig_content[contig] = bac_content
    return contig_content


def write_most_common_percentage(assembly_fasta, scaffolds_directory, contig_content, file_to_write_assembly_content):
    total_length = fasta_length(scaffolds_directory)
    most_common_bact = defaultdict(int)
    for contig in contig_content:
        bact_content = contig_content[contig]
        for bacteria in bact_content:
            length = bact_content[bacteria]
            if length:
                most_common_bact[bacteria] += length[0]

    bacteria_top = []
    for bacteria, length in most_common_bact.items():
        percentage = length / total_length * 100
        bacteria_top.append((bacteria, length, percentage))
    bacteria_top = sorted(bacteria_top, key=lambda bact: bact[1], reverse=True)

    with open(file_to_write_assembly_content, "w") as fw:
        fw.write("# report is created by Contera tool, please cite us https://github.com/zilov/contera!\n")
        fw.write("# bacteria\tlength\tpercentage_in_genome\n")
        for bact in bacteria_top:
            fw.write(f"{bact[0]}\t{bact[1]}\t{bact[2]}\n")


def adapters_blast_matches(blast_adapters_matches_file):
    adapters_to_del = {}
    with open(blast_adapters_matches_file) as fh:
        for line in fh:
            line = line.strip().split()
            ids = line[1].split(":")[1]
            contig = line[0]
            start = line[6]
            stop = line[7]
            adapters_to_del[contig] = start + ":" + stop
    return adapters_to_del


def remove_adapters_and_short_contigs(assembly_fasta, adapters_to_del_dir):
    fasta = fasta_reader(assembly_fasta)
    number_of_adapters = 0
    if adapters_to_del_dir:
        number_of_adapters = len(adapters_to_del_dir.keys())
        for contig, position in adapters_to_del_dir.items():
            start = int(position.split(":")[0]) - 1
            stop = int(position.split(":")[1])
            for key, value in fasta.items():
                if key.startswith(">" + contig):
                    line_to_del = value[start:stop]
                    value = value.replace(line_to_del, "")
                    fasta[key] = value
    else:
        print("No adapters to remove!\n")

    short_contig = []
    for key, value in fasta.items():
        if len(value) < 200:
            short_contig.append(key)

    number_of_short_contigs = len(short_contig)
    if number_of_short_contigs > 0:
        for i in short_contig:
            if i in fasta.keys():
                fasta.pop(i)
    else:
        print("No short contigs! DONE!")

    print(f"{number_of_short_contigs} short contigs and {number_of_adapters} adapters will be deleted!")
    return fasta


def write_filtered_fasta(filtered_fasta_dir, filtered_fasta_file):
    with open(filtered_fasta_file, "w") as fw:
        for key, value in filtered_fasta_dir.items():
            fw.write(f">{key}\n{value}\n")


def run_contig_content_check(assembly_fasta, blast_db, threads, blast_content_file, debug,
                             contigs_content_file, assembly_content_file):

    blasted = run_blast_nt(assembly_fasta, blast_db, threads, blast_content_file, debug)
    if blasted:
        blast_matches_dir = read_blast_matches(blast_content_file)
        bact_percentage = bacteria_percentage_in_contig(blast_matches_dir)
        bact_sequences = bacteria_sequence_in_contig(bact_percentage)
        assembly_fasta_dir = fasta_reader(assembly_fasta)
        contig_content_dir = contig_inside(bact_sequences, assembly_fasta_dir, contigs_content_file)
        write_most_common_percentage(assembly_fasta, assembly_fasta_dir, contig_content_dir, assembly_content_file)
        print("Successfully done! Please cite Contera https://github.com/zilov/contera!\n")
    else:
        print("Blast file was no created! Close all tasks!\n")


def run_adapter_removement(assembly_fasta, adapters_db, threads, adapters_blast_file, debug, assembly_filtered):
    blasted = run_blast_adapters(assembly_fasta, adapters_db, threads, adapters_blast_file, debug)
    if blasted:
        adapters_to_del = adapters_blast_matches(adapters_blast_file)
        filtered_fasta_dir = remove_adapters_and_short_contigs(assembly_fasta, adapters_to_del)
        write_filtered_fasta(filtered_fasta_dir, assembly_filtered)
        print("Successfully done! Please cite Contera https://github.com/zilov/contera!\n")
    else:
        print("Blast file was no created! Close all tasks!")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Contera -  a tool for asseembly contamination detection and erasing')
    parser.add_argument('-m', '--mode', help="mode to use [default = adapter]",
                        choices=["aio", "adapter"], default="adapter")
    parser.add_argument('-a', '--assembly', help="path to genome asssembly in FASTA format", required=True)
    parser.add_argument('-p', '--prefix', help="prefix for output files", default="0")
    parser.add_argument('-o', '--outdir', help='output directory', required=True)
    parser.add_argument('-db', '--blast_db', help='path to blast nucleotide (nt) database', required=False)
    parser.add_argument('-t', '--threads', help='number of threads [default == 8]', default="8")
    parser.add_argument('-d', '--debug', help='do not run, print commands only', action='store_true', default=False)

    args = vars(parser.parse_args())

    assembly_fasta = os.path.abspath(args["assembly"])
    prefix = args["prefix"]
    outdir = os.path.abspath(args["outdir"])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    threads = args["threads"]
    if args["blast_db"]:
        blast_db = os.path.abspath(args["blast_db"])
    mode = args["mode"]
    debug = args["debug"]

    execution_folder = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))
    adapters_db = os.path.join(execution_folder, "adapters_db/adaptor_fasta.fna")
    if prefix == "0":
        prefix = os.path.splitext(os.path.basename(assembly_fasta))[0]

    # define output files
    content_blast_file = os.path.join(outdir, f"{prefix}_content_blast.outfmt6")
    adapters_blast_file = os.path.join(outdir, f"{prefix}_adapters_blast.outfmt6")
    contera_contig_content = os.path.join(outdir, f"{prefix}_contig_content.tsv")
    contera_assembly_content = os.path.join(outdir, f"{prefix}_assembly_content.tsv")
    assembly_filtered = os.path.join(outdir, "assembly_filtered.fasta")

    # Check that all files required for each mode are in and run pipeline
    if mode == "adapter":
        run_adapter_removement(assembly_fasta, adapters_db, threads, adapters_blast_file, debug, assembly_filtered)
    elif mode == "aio":
        run_contig_content_check(assembly_fasta, blast_db, threads, content_blast_file, debug,
                                 contera_contig_content, contera_assembly_content )
        run_adapter_removement(assembly_fasta, adapters_db, threads, adapters_blast_file, debug, assembly_filtered)
