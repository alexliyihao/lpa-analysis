"""Collecting the locus information of Coassin's output

Example::

    lc = locus_collector.LocusCollector(
        input_path = "/some/parent/path/of/bam/output" #choose this or next line
        bam_list = "/paths/to/a/file/recording/bam/path/line/by/line.txt")

    locus_table = lc.generate_locus_table()

"""

import numpy as np
import pandas as pd
import gc
import os
import glob


class LocusCollector:
    """
    Collecting the locus information of Coassin's output

    Initialize:

        lc = LocusCollector(input_path: Optional[str],
                            bam_list: Optional[str])
    Generate Table:

        locus_table = lc.generate_locus_table()
    """

    def __init__(self, input_path: str = None, bam_list: str = None):
        if bam_list is None:
            self.annotated_variant_iter = glob.iglob(
                os.path.join(
                    input_path,
                    "**",
                    "variantsAnnotate",
                    "variantsAnnotate.txt"),
                recursive=True)
        else:
            with open(bam_list, "r") as file:
                self.annotated_variant_iter = (os.path.join(
                    line.rstrip(),
                    "variantsAnnotate",
                    "variantsAnnotate.txt")
                    for line in file)

    def read_locus(self, path: str):
        """read a coassin output, tidy the output"""
        df = pd.read_csv(path, sep="\t")
        df = df[['Pos', "Ref", "Variant", "mylocus","wt","mut"]]
        return df

    def generate_locus_table(self):
        """Generate the locus table in one line

        Returns:

            pd.DataFrame, all the locus result
        """
        locus_iter = (self.read_locus(output_files)
                      for output_files in self.annotated_variant_iter)
        # for memory efficiency
        locus_table = pd.concat(locus_iter, axis=0).drop_duplicates()
        gc.collect()
        # If the mylocus gives Exon421 or 422, give "Exon", otherwise "Intron"
        locus_table["coding"] = np.where(
            locus_table["mylocus"].isin(["Exon421", "Exon422"]),
            "Exon",
            "Intron")
        # generate "pos-ref/var" format
        locus_table["pos-ref/var"] = locus_table["Pos"].astype(str) + "-" + \
                                     locus_table["Ref"].astype(str) + "/" + \
                                     locus_table["Variant"].astype(str)
        locus_table = locus_table.sort_values(by="Pos").reset_index(drop=True)
        locus_table = locus_table.set_index("pos-ref/var")
        return locus_table


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="locus_collector.py",
        description="Collecting the locus information of Coassin's output",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-I', "--input_path", type=str, default=None,
                             help="path to the folder storing Coassin's output folders")
    input_group.add_argument('-L', "--bam_list", type=str, default=None,
                             help="path to a file recording Coassin's output folders path by row")
    parser.add_argument('-O', "--output_path", type=str, required=True,
                        help="path saving the output and intermediate results")
    Args = parser.parse_args()

    lc = LocusCollector(input_path=Args.input_path,
                        bam_list=Args.bam_list)
    locus_table = lc.generate_locus_table()
    locus_table.to_csv(Args.output_path)
