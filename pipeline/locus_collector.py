import numpy as np
import pandas as pd
import gc
import glob

class LocusCollector():
    """
    Collecting the locus information of Coassin's output

    Usage:
    # glob.iglob will search all the subdirectory of parent_path for variantsAnnotate/variantsAnnotate.txt
    lc = LocusCollector(path = parent_path)
    locus_table = lc.generate_locus_table()
    """
    def __init__(self, path: str):
        self.path = path

    def read_locus(self, path: str):
        '''
        read a coassin output, tidy the output
        '''
        df = pd.read_csv(path, sep = "\t")
        df = df[['Pos', "Ref", "Variant", "mylocus"]]
        return(df)

    def generate_locus_table(self):
        annotated_variant_iter = glob.iglob(
            f"{self.path}/**/variantsAnnotate/variantsAnnotate.txt", recursive = True
            )
        locus_iter = (self.read_locus(output_files) for output_files in annotated_variant_iter)
        # for memory efficiency
        locus_table = pd.concat(locus_iter, axis = 0).drop_duplicates()
        gc.collect()
        #If the mylocus gives Exon421 or 422, give "Exon", otherwise "Intron"
        locus_table["coding"] = np.where(locus_table["mylocus"].isin(["Exon421","Exon422"]),
                                         "Exon",
                                         "Intron")
        # generate "pos-ref/var" format
        locus_table["pos-ref/var"] = locus_table["Pos"].astype(str) +  "-" + \
                                        locus_table["Ref"].astype(str) +  "/"+ \
                                        locus_table["Variant"].astype(str)
        locus_table = locus_table.sort_values(by = "Pos").reset_index(drop = True)
        locus_table = locus_table.set_index("pos-ref/var")
        return locus_table
