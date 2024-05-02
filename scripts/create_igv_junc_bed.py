#!/usr/bin/env python

"""Created on Mon Jun 19 15:04:37 MDT 2023

Synopsis:

create_igv_junc.py 1099_SJ.out.tab

Description:

A script to convert the STAR v2 aligner *.SJ.out.tab file to a IGV
splice junction bed file.

"""

def main():
    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser(
        description='A script to convert the STAR v2 aligner *.SJ.out.tab file to a IGV splice junction bed file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file',
                        help='STAR v2 aligner *.SJ.out.tab file')
    args = parser.parse_args()

    df = pd.read_table(args.file, names=('chrom', 'start', 'end', 'strand',
                                    'motif', 'annotated', 'unique',
                                    'multi', 'overhang'))

    # From the STAR maunal (https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
    #-------------------------------------------------------------------------------------------------------------------
    # SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format. Note that
    # STAR defines the junction start/end as intronic bases, while many other software define them as
    # exonic bases. The columns have the following meaning:
    # column 1: chromosome
    # column 2: first base of the intron (1-based)
    # column 3: last base of the intron (1-based)
    # column 4: strand (0: undefined, 1: +, 2: -)
    # column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:
    # AT/AC, 6: GT/AT
    # column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
    # column 7: number of uniquely mapping reads crossing the junction
    # column 8: number of multi-mapping reads crossing the junction
    # column 9: maximum spliced alignment overhang

    # To prepare the data for IGV Splice Junction format we need to do some mapping
    # below we're mapping columns 4, 5 & 6 from SJ.out.tab output above

    # column 2: (start) is 1-based in SJ.out.tab and needs to be 0-based for BED
    df['start'] = df['start'].apply(lambda x: x - 1)

    # column 4: (strand) strand (0: undefined, 1: +, 2: -)
    df['strand'] = df['strand'].map({0: 'undefined', 1: '+', 2: '-'})

    # column 5: (motif) intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
    df['motif'] = df['motif'].map({0: 'non-canonical', 1: 'GT/AG', 2: 'CT/AC',
                                   3: 'GC/AG', 4: 'CT/GC', 5: 'AT/AC',
                                   6: 'GT/AT'})

    # column 6: (annotated_junction) 0: unannotated, 1: annotated (only if splice junctions database is used)
    df['annotated'] = df['annotated'].map({0: 'unannotated', 1: 'annotated'})

    # Writing output for IGV splice juction BED/GFF3 format
    # The format is a hybrid BED/GFF3 format specific to IGV splice junction tracks:
    # Here are links for details on the format we're going to write.
    # https://github.com/igvteam/igv.js/wiki/Splice-Junctions
    # https://software.broadinstitute.org/software/igv/BED
    # http://genome.ucsc.edu/FAQ/FAQformat#format1

    # Example of an IGV splice juction track row:
    # chr15  92883778  92885514  motif=GT/AG;uniquely_mapped=95;multi_mapped=10;maximum_spliced_alignment_overhang=38;annotated_junction=true95+

    # Columns of the IGV splice junction bed file
    # column 1: chromosome
    # column 2: start
    # column 3: end
    # column 4: name (GFF3 formatted attributes from the STAR SJ.out.tab file described above
    # column 5: score (uses the 'column 7: number of uniquely mapping reads crossing the junction' from STAR above)
    # column 6: strand

    # Column 4 attributes examples
    #--------------------
    # motif=GT/AG;
    # uniquely_mapped=95;
    # multi_mapped=10;
    # maximum_spliced_alignment_overhang=38;
    # annotated_junction=true

    # Description from UCSC Genome Browser of their standard BED format which 
    # we follow except for the 4th column, the modified name column.

    # 1) chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
    # 2) chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    # 3) chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature, however, the number in position format will be represented. For example, the first 100 bases of chromosome 1 are defined as chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 0-99 in our software (not 0-100), but will represent the position notation chr1:1-100. Read more here.
    #    chromStart and chromEnd can be identical, creating a feature of length 0, commonly used for insertions. For example, use chromStart=0, chromEnd=0 to represent an insertion before the first nucleotide of a chromosome.
    # 4) name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
    # 5) score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
    # 6) strand - Defines the strand. Either "." (=no strand) or "+" or "-".

    print(df)

    format_str='motif={};uniquely_mapped={};maximum_spliced_alignment_overhang={};annotated_junction={}'
    for (idx, row) in df.iterrows():
        # ('chrom', 'start', 'end', 'strand',
        #  'motif', 'annotated', 'unique',
        #  'multi', 'overhang')
        attributes = format_str.format(*[str(x) for x in row.loc[['motif', 'unique', 'overhang', 'annotated']].values])
        # format_str.format(row.loc[['motif', 'unique', 'overhang', 'annotated']].values)

        row_data = (row['chrom'],
                    str(row['start']),
                    str(row['end']),
                    attributes,
                    str(row['unique']),
                    row['strand'])
    
        print('\t'.join(row_data))

if __name__ == "__main__":
    main()

